#include <string_view>

#include "MOM_state.h"

#include "MOM_logger.h"

namespace MOM {

namespace {

void check_field(const amrex::MultiFab &field, const std::string_view name) {
  if (field.empty()) {
    logger::fatal("State: the ", name, " field is not created.");
  }
}

} // namespace

State::State(StateFields &&fields)
  : fields_(std::move(fields)) {

  check_field(fields_.h, "h");
  check_field(fields_.u, "u");
  check_field(fields_.v, "v");
}

} // namespace MOM
