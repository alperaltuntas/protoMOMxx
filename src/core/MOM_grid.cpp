#include <string_view>

#include "MOM_grid.h"

#include "MOM_logger.h"

namespace MOM {

namespace {

// Check that a grid field is created. Each field is created and computed by
// the same setup function (src/initialization), so a field that is not
// created is one whose setup step was skipped.
void check_field(const amrex::MultiFab &field, const std::string_view name) {
  if (field.empty()) {
    logger::fatal("Grid: the ", name, " field is not created.");
  }
}

} // namespace

Grid::Grid(GridFields &&fields)
  : fields_(std::move(fields)) {

  check_field(fields_.geoLatT, "geoLatT");
  check_field(fields_.geoLonT, "geoLonT");
  check_field(fields_.dxT, "dxT");
  check_field(fields_.dyT, "dyT");
  check_field(fields_.areaT, "areaT");

  check_field(fields_.geoLatCu, "geoLatCu");
  check_field(fields_.geoLonCu, "geoLonCu");
  check_field(fields_.dxCu, "dxCu");
  check_field(fields_.dyCu, "dyCu");

  check_field(fields_.geoLatCv, "geoLatCv");
  check_field(fields_.geoLonCv, "geoLonCv");
  check_field(fields_.dxCv, "dxCv");
  check_field(fields_.dyCv, "dyCv");

  check_field(fields_.geoLatBu, "geoLatBu");
  check_field(fields_.geoLonBu, "geoLonBu");
  check_field(fields_.dxBu, "dxBu");
  check_field(fields_.dyBu, "dyBu");

  check_field(fields_.CoriolisBu, "CoriolisBu");

  check_field(fields_.bathyT, "bathyT");
  check_field(fields_.mask2dT, "mask2dT");
  check_field(fields_.mask2dCu, "mask2dCu");
  check_field(fields_.mask2dCv, "mask2dCv");
  check_field(fields_.mask2dBu, "mask2dBu");

  if (!(fields_.max_depth > 0.0)) {
    logger::fatal("Grid: max_depth is not set.");
  }
}

} // namespace MOM
