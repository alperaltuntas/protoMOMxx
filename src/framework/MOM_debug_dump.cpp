#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>

#include <AMReX_ParallelDescriptor.H>

#include "MOM_debug_dump.h"
#include "MOM_logger.h"

namespace MOM::debug {

void dump_field(const amrex::MultiFab &mf, const std::string &name, const int comp) {

  const int rank = amrex::ParallelDescriptor::MyProc();
  const std::string path = name + ".rank" + std::to_string(rank) + ".txt";

  std::FILE *out = std::fopen(path.c_str(), "w");
  if (out == nullptr) {
    return;
  }

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.validbox();
    const amrex::Array4<const amrex::Real> a = mf.const_array(mfi);
    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          std::fprintf(out, "%d %d %d %.17g %a\n", i, j, k,
                       static_cast<double>(a(i, j, k, comp)),
                       static_cast<double>(a(i, j, k, comp)));
        }
      }
    }
  }

  std::fclose(out);
}

namespace {
bool reporting_on = false;
} // namespace

bool reporting() { return reporting_on; }

void set_reporting(const bool on) { reporting_on = on; }

void report_field(const amrex::MultiFab &mf, const std::string &name, const int comp) {
  if (!reporting_on) return;

  const amrex::IndexType ix = mf.ixType();
  double amin = std::numeric_limits<double>::max();
  double amax = -std::numeric_limits<double>::max();
  double sum = 0.0;
  long npts = 0;
  long long bc = 0;

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    // Drop the low face row in each staggered direction: MOM6 reports a
    // velocity field over its non-symmetric index range.
    amrex::Box bx = mfi.validbox();
    if (ix.nodeCentered(0)) bx.growLo(0, -1);
    if (ix.nodeCentered(1)) bx.growLo(1, -1);

    const amrex::Array4<const amrex::Real> a = mf.const_array(mfi);
    for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          const double val = static_cast<double>(a(i, j, k, comp));
          amin = std::min(amin, val);
          amax = std::max(amax, val);
          sum += val;
          ++npts;
          std::uint64_t bits;
          const double av = std::abs(val);
          std::memcpy(&bits, &av, sizeof(bits));
          bc += __builtin_popcountll(bits);
        }
      }
    }
  }

  char buf[320];
  std::snprintf(buf, sizeof(buf),
                "%-24s mean=%23.16E min=%23.16E max=%23.16E c=%9lld",
                name.c_str(), sum / static_cast<double>(npts), amin, amax,
                bc % 1000000000LL);
  logger::note(buf);
}

} // namespace MOM::debug
