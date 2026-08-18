#include <optional>
#include <string>
#include <vector>

#include "MOM_domains.h"
#include "MOM_logger.h"
#include "MOM_unsupported_params.h"

namespace MOM {

Domain make_domain(RuntimeParams &params) {

  params.doc_module("MOM_domains", "");

  bool reentrant_x = true;
  bool reentrant_y = false;
  bool tripolar_n = false;
  int ni_global = 0;
  int nj_global = 0;
  int ni_halo = 4;
  int nj_halo = 4;

  params.get("REENTRANT_X", reentrant_x,
             {.default_value = true,
              .desc = "If true, the domain is zonally reentrant."});

  params.get("REENTRANT_Y", reentrant_y,
             {.default_value = false,
              .desc = "If true, the domain is meridionally reentrant."});

  params.get("TRIPOLAR_N", tripolar_n,
             {.default_value = false,
              .desc = "Use tripolar connectivity at the northern edge of the domain. "
                      "With TRIPOLAR_N, NIGLOBAL must be even."});

  if (reentrant_y && tripolar_n) {
    logger::fatal("make_domain: REENTRANT_Y and TRIPOLAR_N cannot both be true.");
  }

  // defer: tripolar connectivity.
  if (tripolar_n) {
    logger::fatal("make_domain: TRIPOLAR_N is not implemented yet.");
  }

  params.get("NIGLOBAL", ni_global,
             {.desc = "The total number of thickness grid points in "
                      "the x-direction in the physical domain.",
              .fail_if_missing = true});

  params.get("NJGLOBAL", nj_global,
             {.desc = "The total number of thickness grid points in "
                      "the y-direction in the physical domain.",
              .fail_if_missing = true});

  if (ni_global < 1 || nj_global < 1) {
    logger::fatal("make_domain: NIGLOBAL and NJGLOBAL must be positive.");
  }

  params.get("NIHALO", ni_halo,
             {.default_value = 4,
              .desc = "The number of halo points on each side in the x-direction."});

  params.get("NJHALO", nj_halo,
             {.default_value = 4,
              .desc = "The number of halo points on each side in the y-direction."});

  if (ni_halo < 0 || nj_halo < 0) {
    logger::fatal("make_domain: NIHALO and NJHALO must be non-negative.");
  }

  // A protoMOMxx-only knob with no MOM6 counterpart: MOM6 fixes the
  // decomposition to one tile per PE and varies it through LAYOUT and the PE
  // count. Splitting a rank's share into several boxes is how the layout
  // independence of DESIGN.md section 5 gets exercised without MPI, which is
  // the only way to test halo handling in a serial build.
  int n_boxes = 0;
  params.get("NBOXES", n_boxes,
             {.default_value = 0,
              .desc = "The number of boxes to decompose the domain into. The default of 0 "
                      "gives one box per processor. Any other value is a debugging aid: the "
                      "answers must not depend on it.",
              .debugging_param = true});
  if (n_boxes < 0) {
    logger::fatal("make_domain: NBOXES must not be negative.");
  }

  unsupported_param(params, "GLOBAL_INDEXING", true,
                    "protoMOMxx always uses global index conventions; local "
                    "indexing is not supported.");
  unsupported_param(params, "LAYOUT", std::vector<int>{0, 0},
                    "the decomposition and its processor assignment are set "
                    "automatically.");
  unsupported_param(params, "IO_LAYOUT", std::vector<int>{1, 1},
                    "writing distributed files via an I/O processor layout is "
                    "not supported.");
  unsupported_param(params, "MASKTABLE", std::string{},
                    "masking out land-only processors via a table is not "
                    "supported.");
  unsupported_param(params, "AUTO_MASKTABLE", false,
                    "masking out land-only processors via a table is not "
                    "supported.");

  return Domain({.ni_global = ni_global,
                 .nj_global = nj_global,
                 .ni_halo = ni_halo,
                 .nj_halo = nj_halo,
                 .reentrant_x = reentrant_x,
                 .reentrant_y = reentrant_y,
                 .tripolar_n = tripolar_n,
                 .n_boxes = (n_boxes > 0) ? std::optional<int>(n_boxes) : std::nullopt});
}

} // namespace MOM
