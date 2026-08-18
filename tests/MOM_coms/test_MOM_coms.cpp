// Unit tests for the MOM_coms module (src/framework/MOM_coms.cpp): the
// order-invariant global sum. The tests check the properties the sum exists
// for -- exactness for a single value, and independence of the order the
// terms are added in -- rather than re-deriving the fixed-point arithmetic.
// (Hence the main() below, which instantiates a MOM::Infra instance.)

#include <algorithm>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "MOM_coms.h"
#include "MOM_domain_infra.h"
#include "MOM_infra.h"

// A single value round-trips exactly, for magnitudes spanning the range the
// six integers cover.
TEST(MOMComsTest, SingleValueRoundTrips) {
  for (const double value : {0.0, 1.0, -1.0, 1234.5678, -9.87654321e12, 3.5e-9, -2.25e-15}) {
    MOM::EFP efp;
    efp.add(value);
    EXPECT_EQ(efp.to_real(), value) << "value " << value;
  }
}

// The sum does not depend on the order the terms are added in, which is what
// makes it independent of the domain decomposition.
TEST(MOMComsTest, SumIsOrderInvariant) {
  std::vector<double> values;
  std::mt19937 gen(12345);
  std::uniform_real_distribution<double> dist(-1.0e6, 1.0e6);
  values.reserve(500);
  for (int i = 0; i < 500; ++i) {
    values.push_back(dist(gen));
  }

  MOM::EFP forward;
  for (const double v : values) forward.add(v);

  std::vector<double> shuffled = values;
  std::shuffle(shuffled.begin(), shuffled.end(), gen);
  MOM::EFP shuffled_sum;
  for (const double v : shuffled) shuffled_sum.add(v);

  EXPECT_EQ(forward.to_real(), shuffled_sum.to_real());
  // A naive sum is not order invariant at this width, so the test would be
  // vacuous if it were.
  double naive_forward = 0.0;
  for (const double v : values) naive_forward += v;
  double naive_shuffled = 0.0;
  for (const double v : shuffled) naive_shuffled += v;
  EXPECT_NE(naive_forward, naive_shuffled);
}

// The sum of a field is the same however the domain is split into boxes.
TEST(MOMComsTest, FieldSumIsLayoutIndependent) {
  const auto sum_with = [](const int n_boxes) {
    const MOM::Domain domain({.ni_global = 44,
                              .nj_global = 40,
                              .ni_halo = 4,
                              .nj_halo = 4,
                              .n_boxes = n_boxes});
    amrex::MultiFab field = domain.make_field(MOM::Stagger::Cell, 2, 1);
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
      const amrex::Box bx = mfi.validbox();
      const amrex::Array4<amrex::Real> a = field.array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        a(i, j, k) = std::sin(0.1 * i) * std::cos(0.07 * j) * (1.0 + k) * 1.0e7;
      });
    }
    return MOM::reproducing_sum(field);
  };

  const double one_box = sum_with(1);
  EXPECT_EQ(sum_with(2), one_box);
  EXPECT_EQ(sum_with(4), one_box);
  EXPECT_EQ(sum_with(8), one_box);
}

/// @brief Test-binary entry point: initialize GTest, bring up the
/// infrastructure layer, and run all tests.
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  const MOM::Infra infra(argc, argv);
  return RUN_ALL_TESTS();
}
