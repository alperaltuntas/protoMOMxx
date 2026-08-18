#include <cstdio>

#include <AMReX_ParallelDescriptor.H>

#include "MOM_debug_dump.h"

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

} // namespace MOM::debug
