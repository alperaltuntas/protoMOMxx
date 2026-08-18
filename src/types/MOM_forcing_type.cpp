#include "MOM_forcing_type.h"

namespace MOM {

MechForcing::MechForcing(const Domain &domain)
  : taux_(domain.make_field(Stagger::XFace, 1, 1)),
    tauy_(domain.make_field(Stagger::YFace, 1, 1)) {
  taux_.setVal(0.0);
  tauy_.setVal(0.0);
}

} // namespace MOM
