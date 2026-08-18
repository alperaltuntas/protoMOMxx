// Unit tests for the VerticalGrid class (src/core/MOM_vertical_grid.cpp),
// which is constructed from an injected RuntimeParams object. The vertical
// grid is AMReX-free, so no infra (AMReX) bring-up is needed.

#include <filesystem>
#include <gtest/gtest.h>

#include "MOM_logger.h"
#include "MOM_vertical_grid.h"

// Parse a parameter file given relative to this test's directory.
MOM::RuntimeParams params_from(const std::string &rel_path) {
  return MOM::RuntimeParams((std::filesystem::path(__FILE__).parent_path() / rel_path).string());
}

// The gprime coordinate from the driving testcase's parameter file (reused
// here rather than duplicated): the surface interface takes GFS, the internal
// interfaces take GINT, and the layer densities build upward from RHO_0 via
// Rlay[k] = Rlay[k-1] + g_prime[k]*(Rho0/g_Earth).
TEST(MOMVerticalGridTest, ConstructsFromInjectedParams) {
  auto params = params_from("../double_gyre/MOM_input");
  const MOM::VerticalGrid vgrid(params);

  EXPECT_EQ(vgrid.nk(), 2);
  EXPECT_DOUBLE_EQ(vgrid.max_depth(), 2000.0);
  EXPECT_DOUBLE_EQ(vgrid.g_Earth(), 9.8);   // default
  EXPECT_DOUBLE_EQ(vgrid.Rho0(), 1035.0);   // default

  ASSERT_EQ(vgrid.g_prime().size(), 3u);
  EXPECT_DOUBLE_EQ(vgrid.g_prime()[0], 0.98);    // GFS
  EXPECT_DOUBLE_EQ(vgrid.g_prime()[1], 0.0098);  // GINT
  EXPECT_DOUBLE_EQ(vgrid.g_prime()[2], 98.0);    // bottom filler: 10*g_Earth

  ASSERT_EQ(vgrid.Rlay().size(), 2u);
  EXPECT_DOUBLE_EQ(vgrid.Rlay()[0], 1035.0);  // LIGHTEST_DENSITY default = RHO_0
  EXPECT_DOUBLE_EQ(vgrid.Rlay()[1], 1035.0 + 0.0098 * (1035.0 / 9.8));
}

// GFS defaults to G_EARTH and LIGHTEST_DENSITY defaults to RHO_0; with
// NK = 1 there are no internal interfaces to fill.
TEST(MOMVerticalGridTest, GprimeDefaultsAndSingleLayer) {
  auto params = params_from("MOM_param_files/MOM_input_gprime_defaults");
  const MOM::VerticalGrid vgrid(params);

  EXPECT_EQ(vgrid.nk(), 1);
  EXPECT_DOUBLE_EQ(vgrid.g_prime()[0], 9.8);  // GFS default = G_EARTH
  EXPECT_DOUBLE_EQ(vgrid.Rlay()[0], 1035.0);  // LIGHTEST_DENSITY default = RHO_0
}

// Unsupported COORD_CONFIG options and nonsensical values (NK = 0 standing in
// for the value checks) are fatal. Missing mandatory parameters throwing is
// the parameter table's contract, tested in test_MOM_file_parser -- it is
// deliberately not re-verified here.
TEST(MOMVerticalGridTest, InvalidConfigurationsAreFatal) {
  auto bad_coord = params_from("MOM_param_files/MOM_input_bad_coord");
  EXPECT_THROW(MOM::VerticalGrid{bad_coord}, MOM::logger::FatalError);

  auto bad_nk = params_from("MOM_param_files/MOM_input_bad_nk");
  EXPECT_THROW(MOM::VerticalGrid{bad_nk}, MOM::logger::FatalError);
}
