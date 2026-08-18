// Unit tests for make_grid (src/initialization/MOM_fixed_initialization.cpp):
// the parameter-driven computation of the grid fields, including defaults and
// the rejection of invalid extents. The grid values themselves are covered in
// tests/MOM_grid.
//
// make_grid requires the infra layer to be initialized, hence the
// main() below, which instantiates a MOM::Infra.

#include <filesystem>
#include <gtest/gtest.h>

#include "MOM_domains.h"
#include "MOM_fixed_initialization.h"
#include "MOM_grid.h"
#include "MOM_infra.h"
#include "MOM_logger.h"

namespace {

// Helper function to get the absolute path to the test data directory
std::filesystem::path get_test_data_dir() {
  return std::filesystem::path(__FILE__).parent_path() / "MOM_param_files";
}

// Construct the domain and the grid from one parameter file.
MOM::Grid grid_from_param_file(const std::string &file_name) {
  const auto path = get_test_data_dir() / file_name;
  EXPECT_TRUE(std::filesystem::exists(path)) << "Test file " << path << " does not exist";
  MOM::RuntimeParams params(path.string());
  const MOM::Domain domain = MOM::make_domain(params);
  return MOM::make_grid(domain, params);
}

} // namespace

// make_grid computes the grid fields from the parameter file, with
// defaults applied for the parameters the file omits (WESTLON, RAD_EARTH,
// ROTATION, OMEGA).
TEST(MOMFixedInitTest, MakeGridFromParamFile) {
  const MOM::Grid grid = grid_from_param_file("MOM_input_test");

  EXPECT_DOUBLE_EQ(grid.south_lat(), 30.0);
  EXPECT_DOUBLE_EQ(grid.len_lat(), 20.0);
  EXPECT_DOUBLE_EQ(grid.len_lon(), 22.0);
  EXPECT_DOUBLE_EQ(grid.west_lon(), 0.0);       // default
  EXPECT_DOUBLE_EQ(grid.rad_earth(), 6.378e6);  // default
}

// Invalid extents are rejected through logger::fatal:
TEST(MOMFixedInitTest, InvalidExtentIsFatal) {
  EXPECT_THROW(grid_from_param_file("MOM_input_bad_extent"), MOM::logger::FatalError);
}

/// @brief Test-binary entry point: initialize GTest, bring up the
/// infrastructure layer, and run all tests.
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  const MOM::Infra infra(argc, argv);
  return RUN_ALL_TESTS();
}
