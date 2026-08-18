#pragma once
/// @file MOM_coms.h
/// @brief Order-invariant global sums. The analogue of MOM6's MOM_coms
///        (src/framework/MOM_coms.F90), which carries the same facility.

#include <array>
#include <cstdint>
#include <vector>

#include <AMReX_MultiFab.H>

namespace MOM {

/// @class EFP
/// @brief A real number held as a fixed-point value in six 64-bit integers.
///
/// Adding two of these is integer addition, which is exact and associative,
/// so a sum of many terms does not depend on the order they are added in and
/// therefore does not depend on the domain decomposition. The technique is
/// Hallberg & Adcroft, 2014, Parallel Computing 40(5-6),
/// doi:10.1016/j.parco.2014.04.007, and the layout of the six integers here
/// is MOM6's, so a value round-trips to the same double MOM6 would report.
///
/// The digits are only canonicalized on conversion back to a real number, so
/// intermediate representations may differ from MOM6's while the value does
/// not.
class EFP {
public:
  /// @brief The number of 64-bit integers a value is spread over.
  static constexpr int n_digits = 6;

  /// @brief Add a real number to this value, without carrying between digits.
  /// @param value The number to add.
  void add(double value);

  /// @brief Carry between digits so that each holds less than one unit of the
  /// digit above it. Does not change the value.
  void carry();

  /// @brief Carry, then give every digit the sign of the overall value. Does
  /// not change the value, and yields the same digits for equal values.
  void regularize();

  /// @brief The real number this value represents.
  /// @return The value as a double.
  double to_real() const;

  /// @brief Add another value to this one.
  /// @param other The value to add.
  /// @return This value.
  EFP &operator+=(const EFP &other);

  /// @brief The difference of two values.
  /// @param other The value to subtract.
  /// @return The difference.
  EFP operator-(const EFP &other) const;

  /// @brief The six integers holding the value.
  /// @return The digits, most significant first.
  std::array<std::int64_t, n_digits> &digits() { return v_; }
  /// @brief The six integers holding the value.
  /// @return The digits, most significant first.
  const std::array<std::int64_t, n_digits> &digits() const { return v_; }

private:
  std::array<std::int64_t, n_digits> v_{};  ///< The digits, most significant first.
};

/// @brief Sum a field over its valid region and over all ranks, in a way that
/// does not depend on the decomposition. The analogue of MOM6's
/// reproducing_sum.
/// @param mf The field to sum. Every level of every valid box is included.
/// @param layer_sums If given, filled with the sum of each k level.
/// @param efp_sum If given, set to the total in fixed-point form, so that
///        differences between calls can be taken without losing digits.
/// @return The total.
double reproducing_sum(const amrex::MultiFab &mf,
                       std::vector<double> *layer_sums = nullptr,
                       EFP *efp_sum = nullptr);

} // namespace MOM
