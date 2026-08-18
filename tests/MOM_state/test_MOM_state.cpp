// Sanity tests for the prognostic state and its initialization
// (src/types/MOM_state.cpp, src/initialization/MOM_state_initialization.cpp):
// staggering, the uniform thickness column, and rejection of unsupported
// configurations.
//
// The state lives on the domain's decomposition, so the infra layer must be
// up; hence the main() below.

#include <filesystem>
#include <utility>
#include <gtest/gtest.h>

#include "MOM_domain_infra.h"
#include "MOM_infra.h"
#include "MOM_logger.h"
#include "MOM_grid.h"
#include "MOM_grid_initialize.h"
#include "MOM_shared_initialization.h"
#include "MOM_state.h"
#include "MOM_state_initialization.h"
#include "MOM_vertical_grid.h"

namespace {

constexpr int NI = 44;
constexpr int NJ = 40;
constexpr int HALO = 2;
constexpr int NK = 2;
constexpr double MAX_DEPTH = 2000.0;   // [m]
constexpr double ANGSTROM = 1.0e-10;   // [m]

std::filesystem::path param_file(const std::string &name) {
  return std::filesystem::path(__FILE__).parent_path() / "MOM_param_files" / name;
}

MOM::Domain make_domain() {
  return MOM::Domain({.ni_global = NI, .nj_global = NJ,
                      .ni_halo = HALO, .nj_halo = HALO, .reentrant_x = false});
}

// A grid with a flat bottom at a chosen depth, so the expected column is
// exact. max_depth is the vertical extent the layers are distributed over and
// is independent of the bathymetry.
MOM::Grid flat_grid(const MOM::Domain &domain, const double depth) {
  const MOM::GridSpec spec = {.south_lat = 30.0, .len_lat = 20.0,
                              .west_lon = 0.0, .len_lon = 22.0,
                              .rad_earth = 6.378e6, .omega = 7.2921e-5};
  const MOM::TopoSpec topo = {.max_depth = depth, .min_depth = 0.0};
  MOM::GridFields fields = MOM::spherical_grid_fields(domain, spec);
  fields.CoriolisBu = MOM::planetary_rotation(domain, spec, fields.geoLatBu);
  fields.bathyT = MOM::named_topography(domain, "flat", spec, topo,
                                        fields.geoLonT, fields.geoLatT);
  fields.max_depth = MAX_DEPTH;
  MOM::initialize_masks(domain, topo, fields.bathyT, fields);
  MOM::set_derived_metrics(domain, fields);
  return MOM::Grid(std::move(fields));
}

} // namespace

// On a bottom deeper than MAXIMUM_DEPTH/nk everywhere, the uniform
// configuration gives each layer its even share except the bottom one, which
// takes the rest of the water column.
TEST(State, UniformThicknessOnADeepFlatBottom) {
  const MOM::Domain domain = make_domain();
  MOM::RuntimeParams params(param_file("MOM_input_test").string());

  const MOM::Grid grid = flat_grid(domain, MAX_DEPTH);
  const MOM::VerticalGrid vgrid(params);
  MOM::StateFields fields = MOM::initialize_state(domain, grid, vgrid, params);
  const MOM::State state(std::move(fields));

  // Each field sits at its C-grid point type and carries nk layers.
  EXPECT_TRUE(state.h().ixType().cellCentered());
  EXPECT_EQ(state.u().ixType(), amrex::IndexType(amrex::IntVect(1, 0, 0)));
  EXPECT_EQ(state.v().ixType(), amrex::IndexType(amrex::IntVect(0, 1, 0)));
  EXPECT_EQ(state.h().nComp(), 1);
  EXPECT_EQ(state.h().boxArray()[0].length(2), NK);

  // Both layers are MAX_DEPTH/nk thick over the whole domain.
  EXPECT_DOUBLE_EQ(state.h().min(0), MAX_DEPTH / NK);
  EXPECT_DOUBLE_EQ(state.h().max(0), MAX_DEPTH / NK);

  // The fluid starts at rest.
  EXPECT_DOUBLE_EQ(state.u().norm0(), 0.0);
  EXPECT_DOUBLE_EQ(state.v().norm0(), 0.0);
}

// Where the water column is shallower than the resting interface depth, the
// lower layers collapse to an Angstrom and the top layer takes what is left.
TEST(State, UniformThicknessCollapsesShallowColumns) {
  const MOM::Domain domain = make_domain();
  MOM::RuntimeParams params(param_file("MOM_input_test").string());

  const double depth = 300.0;  // shallower than MAX_DEPTH/nk = 1000
  const MOM::Grid grid = flat_grid(domain, depth);
  const MOM::VerticalGrid vgrid(params);
  MOM::StateFields fields = MOM::initialize_state(domain, grid, vgrid, params);
  const MOM::State state(std::move(fields));

  EXPECT_DOUBLE_EQ(state.h().min(0), ANGSTROM);
  EXPECT_DOUBLE_EQ(state.h().max(0), depth - ANGSTROM);
}

// Unsupported configurations abort rather than falling through to a default.
TEST(State, UnsupportedThicknessConfigIsFatal) {
  const MOM::Domain domain = make_domain();
  MOM::RuntimeParams params(param_file("MOM_input_bad_config").string());

  const MOM::Grid grid = flat_grid(domain, MAX_DEPTH);
  const MOM::VerticalGrid vgrid(params);
  EXPECT_THROW(MOM::initialize_state(domain, grid, vgrid, params),
               MOM::logger::FatalError);
}

/// @brief Test-binary entry point: initialize GTest, bring up the
/// infrastructure layer, and run all tests.
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  const MOM::Infra infra(argc, argv);
  return RUN_ALL_TESTS();
}
