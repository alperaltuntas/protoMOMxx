// Unit test for the Model initialization skeleton (src/core/MOM.cpp).
//
// A Model's only constructor argument is a RuntimeParams object, which is
// injected by the driver. Additionally, the Model implicitly depends on
// the infra (AMReX) layer being initialized. (Hence the main() below).

#include <filesystem>
#include <gtest/gtest.h>
#include <string>

#include "MOM.h"
#include "MOM_infra.h"

// Helper function to get the absolute path to the test data directory
std::filesystem::path get_test_data_dir() {
  return std::filesystem::path(__FILE__).parent_path() / "MOM_param_files";
}

// Check if the Model can be constructed from a RuntimeParams object.
TEST(MOMModelTest, ConstructsFromInjectedParams) {
  auto test_file_path = get_test_data_dir() / "MOM_input_test";
  ASSERT_TRUE(std::filesystem::exists(test_file_path)) << "Test file " << test_file_path << " does not exist";

  MOM::RuntimeParams params(test_file_path.string());
  const MOM::Model model(params);

  EXPECT_FALSE(model.config().split);
  EXPECT_TRUE(model.config().use_RK2);
  EXPECT_EQ(model.domain().ni_global(), 44);
  EXPECT_EQ(model.domain().nj_global(), 40);
  EXPECT_DOUBLE_EQ(model.grid().south_lat(), 30.0);  // grid wired from params
  EXPECT_EQ(model.vertical_grid().nk(), 2);
}

// Two Model instances can coexist in one process (multi-instance mode).
// Each is constructed from its own injected RuntimeParams on the shared
// infra runtime, and each ends up with its own independent domain.
// todo: Note that each RuntimeParams here is constructed without a parameter-doc
// writer; two instances documenting to the same MOM_parameter_doc files
// would collide, so per-instance doc file naming is expected to become the
// first real friction point when instances gain doc writers.
TEST(MOMModelTest, TwoInstancesCoexist) {
  auto params_path_a = get_test_data_dir() / "MOM_input_test";
  auto params_path_b = get_test_data_dir() / "MOM_input_test_smaller";
  ASSERT_TRUE(std::filesystem::exists(params_path_a));
  ASSERT_TRUE(std::filesystem::exists(params_path_b));

  MOM::RuntimeParams params_a(params_path_a.string());
  MOM::RuntimeParams params_b(params_path_b.string());

  const MOM::Model model_a(params_a);
  const MOM::Model model_b(params_b);

  // The two instances hold independent configurations and domains.
  EXPECT_FALSE(model_a.config().split);
  EXPECT_TRUE(model_b.config().split);
  EXPECT_EQ(model_a.domain().ni_global(), 44);
  EXPECT_EQ(model_a.domain().nj_global(), 40);
  EXPECT_EQ(model_b.domain().ni_global(), 14);
  EXPECT_EQ(model_b.domain().nj_global(), 10);
  EXPECT_TRUE(model_b.domain().reentrant_x());  // default in b, off in a
  EXPECT_FALSE(model_a.domain().reentrant_x());
  EXPECT_EQ(model_a.vertical_grid().nk(), 2);
  EXPECT_EQ(model_b.vertical_grid().nk(), 3);
}

/// @brief Test-binary entry point: initialize GTest, bring up the infrastructure
/// layer, and run all tests.
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  const MOM::Infra infra(argc, argv);
  return RUN_ALL_TESTS();
}
