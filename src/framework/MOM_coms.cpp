#include "MOM_coms.h"

#include <cmath>
#include <cstdlib>

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include "MOM_logger.h"

namespace MOM {

namespace {

// The value of one unit of digit 2 (the digit holding the ones place is 3).
// MOM6 calls this prec; the six digits step down from it by this factor.
constexpr std::int64_t PREC = std::int64_t(1) << 46;
constexpr double R_PREC = 70368744177664.0;  // 2^46 as a real.

// The real value of one unit of each digit, and its inverse. Every entry is a
// power of two, so scaling by them is exact.
constexpr std::array<double, EFP::n_digits> PR = {
    R_PREC * R_PREC, R_PREC, 1.0, 1.0 / R_PREC, 1.0 / (R_PREC * R_PREC),
    1.0 / (R_PREC * R_PREC * R_PREC)};
constexpr std::array<double, EFP::n_digits> I_PR = {
    1.0 / (R_PREC * R_PREC), 1.0 / R_PREC, 1.0, R_PREC, R_PREC * R_PREC,
    R_PREC * R_PREC * R_PREC};

// The largest magnitude with a fixed-point representation: the top digit is
// the only one allowed to exceed PREC, and it is still a signed 64-bit value.
const double MAX_EFP_FLOAT = PR[0] * (9223372036854775808.0 - 1.0);

} // namespace

void EFP::add(const double value) {
  if (std::isnan(value)) {
    logger::fatal("reproducing_sum: NaN in the field being summed.");
  }
  if (std::fabs(value) > MAX_EFP_FLOAT) {
    logger::fatal("reproducing_sum: ", value, " is too large to represent.");
  }

  const std::int64_t sgn = (value < 0.0) ? -1 : 1;
  double rs = std::fabs(value);
  for (int i = 0; i < n_digits; ++i) {
    const std::int64_t ival = static_cast<std::int64_t>(rs * I_PR[i]);
    rs -= static_cast<double>(ival) * PR[i];
    v_[i] += sgn * ival;
  }
}

void EFP::carry() {
  // Integer division here where MOM6 goes through a real multiply by 1/PREC.
  // Both preserve the value, which is all that the final conversion needs.
  for (int i = n_digits - 1; i >= 1; --i) {
    if (std::llabs(v_[i]) >= PREC) {
      const std::int64_t num_carry = v_[i] / PREC;
      v_[i] -= num_carry * PREC;
      v_[i - 1] += num_carry;
    }
  }
}

void EFP::regularize() {
  carry();

  bool positive = true;
  for (int i = 0; i < n_digits; ++i) {
    if (v_[i] != 0) {
      positive = (v_[i] > 0);
      break;
    }
  }

  if (positive) {
    for (int i = n_digits - 1; i >= 1; --i) {
      if (v_[i] < 0) {
        v_[i] += PREC;
        v_[i - 1] -= 1;
      }
    }
  } else {
    for (int i = n_digits - 1; i >= 1; --i) {
      if (v_[i] > 0) {
        v_[i] -= PREC;
        v_[i - 1] += 1;
      }
    }
  }
}

double EFP::to_real() const {
  EFP tmp = *this;
  tmp.regularize();
  double r = 0.0;
  for (int i = 0; i < n_digits; ++i) {
    r += PR[i] * static_cast<double>(tmp.v_[i]);
  }
  return r;
}

EFP &EFP::operator+=(const EFP &other) {
  for (int i = n_digits - 1; i >= 1; --i) {
    v_[i] += other.v_[i];
    if (v_[i] > PREC) {
      v_[i] -= PREC;
      v_[i - 1] += 1;
    } else if (v_[i] < -PREC) {
      v_[i] += PREC;
      v_[i - 1] -= 1;
    }
  }
  v_[0] += other.v_[0];
  return *this;
}

EFP EFP::operator-(const EFP &other) const {
  EFP result;
  for (int i = 0; i < n_digits; ++i) {
    result.v_[i] = -other.v_[i];
  }
  result += *this;
  return result;
}

double reproducing_sum(const amrex::MultiFab &mf, std::vector<double> *layer_sums,
                       EFP *efp_sum) {
  const amrex::Box covered = mf.boxArray().minimalBox();
  const int k_lo = covered.smallEnd(2);
  const int nk = covered.length(2);

  std::vector<EFP> sums(nk);

  // The fixed-point conversion is serial per digit and the fields here are
  // small, so this stays on the host rather than going through ParallelFor.
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const amrex::Box &bx = mfi.validbox();
    const amrex::Array4<const amrex::Real> a = mf.const_array(mfi);
    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
      EFP &acc = sums[k - k_lo];
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          acc.add(a(i, j, k));
        }
      }
      acc.carry();
    }
  }

  // Integer addition, so the result is the same however the boxes were spread
  // over the ranks.
  std::vector<std::int64_t> flat(static_cast<std::size_t>(nk) * EFP::n_digits);
  for (int k = 0; k < nk; ++k) {
    for (int d = 0; d < EFP::n_digits; ++d) {
      flat[static_cast<std::size_t>(k) * EFP::n_digits + d] = sums[k].digits()[d];
    }
  }
  if (amrex::ParallelContext::NProcsSub() > 1) {
    amrex::ParallelAllReduce::Sum(flat.data(), static_cast<int>(flat.size()),
                                  amrex::ParallelContext::CommunicatorSub());
  }

  if (layer_sums != nullptr) {
    layer_sums->assign(static_cast<std::size_t>(nk), 0.0);
  }
  EFP total;
  double sum = 0.0;
  for (int k = 0; k < nk; ++k) {
    EFP &acc = sums[k];
    for (int d = 0; d < EFP::n_digits; ++d) {
      acc.digits()[d] = flat[static_cast<std::size_t>(k) * EFP::n_digits + d];
    }
    acc.regularize();
    const double val = acc.to_real();
    if (layer_sums != nullptr) {
      (*layer_sums)[static_cast<std::size_t>(k)] = val;
    }
    sum += val;
    total += acc;
  }
  if (efp_sum != nullptr) {
    *efp_sum = total;
  }
  return sum;
}

} // namespace MOM
