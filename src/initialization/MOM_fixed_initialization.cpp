#include <string>

#include "MOM_fixed_initialization.h"
#include "MOM_grid_initialize.h"
#include "MOM_logger.h"
#include "MOM_shared_initialization.h"

namespace MOM {

namespace {

// Read the extents of a simple spherical grid into the spec. The analogue of
// the parameter reads of MOM6's set_grid_metrics_spherical; the metric
// computation itself lives in MOM_grid_initialize.
void read_spherical_grid_params(RuntimeParams &params, GridSpec &spec) {

  params.get("SOUTHLAT", spec.south_lat,
             {.desc = "The southern latitude of the domain.",
              .units = "degrees_N",
              .fail_if_missing = true});

  params.get("LENLAT", spec.len_lat,
             {.desc = "The latitudinal length of the domain.",
              .units = "degrees_N",
              .fail_if_missing = true});

  params.get("WESTLON", spec.west_lon,
             {.default_value = 0.0,
              .desc = "The western longitude of the domain.",
              .units = "degrees_E"});

  params.get("LENLON", spec.len_lon,
             {.desc = "The longitudinal length of the domain.",
              .units = "degrees_E",
              .fail_if_missing = true});

  params.get("RAD_EARTH", spec.rad_earth,
             {.default_value = 6.378e6,
              .desc = "The radius of the Earth.",
              .units = "m"});

  if (!(spec.len_lat > 0.0) || !(spec.len_lon > 0.0)) {
    logger::fatal("initialize_fixed: LENLAT and LENLON must be positive.");
  }
  if (!(spec.rad_earth > 0.0)) {
    logger::fatal("initialize_fixed: RAD_EARTH must be positive.");
  }
}

// Read the planetary rotation configuration into the spec. The analogue of
// the parameter reads of MOM6's MOM_initialize_rotation
// (MOM_shared_initialization.F90).
void read_rotation_params(RuntimeParams &params, GridSpec &spec) {

  std::string rotation = "2omegasinlat";
  params.get("ROTATION", rotation,
             {.default_value = std::string("2omegasinlat"),
              .desc = "This specifies how the Coriolis parameter is specified:\n"
                      "\t 2omegasinlat - Use twice the planetary rotation rate\n"
                      "\t\t times the sine of latitude.\n"
                      "\t betaplane - Use a beta-plane or f-plane.\n"
                      "\t USER - call a user modified routine."});

  if (rotation == "2omegasinlat") {
    params.get("OMEGA", spec.omega,
               {.default_value = 7.2921e-5,
                .desc = "The rotation rate of the earth.",
                .units = "s-1"});
  } else if (rotation == "beta" || rotation == "betaplane") {
    // defer: the beta-plane/f-plane rotation (set_rotation_beta_plane).
    logger::fatal("initialize_fixed: ROTATION \"", rotation,
                  "\" is not implemented yet.");
  } else {
    logger::fatal("initialize_fixed: Unrecognized rotation setup \"", rotation, "\".");
  }
}

// Read the topography configuration into the spec. The analogue of the
// parameter reads of MOM6's MOM_initialize_topography
// (MOM_fixed_initialization.F90) and initialize_topography_named
// (MOM_shared_initialization.F90).
std::string read_topography_params(RuntimeParams &params, TopoSpec &spec) {

  std::string config;
  params.get("TOPO_CONFIG", config,
             {.desc = "This specifies how bathymetry is specified:\n"
                      "\t file - read bathymetric information from the file\n"
                      "\t\t specified by (TOPO_FILE).\n"
                      "\t flat - flat bottom set to MAXIMUM_DEPTH.\n"
                      "\t bowl - an analytically specified bowl-shaped basin\n"
                      "\t\t ranging between MAXIMUM_DEPTH and MINIMUM_DEPTH.\n"
                      "\t spoon - a similar shape to 'bowl', but with an vertical\n"
                      "\t\t wall at the southern face.\n"
                      "\t halfpipe - a zonally uniform channel with a half-sine\n"
                      "\t\t profile in the meridional direction.\n"
                      "\t USER - call a user modified routine.",
              .fail_if_missing = true});

  // MOM6 reads MAXIMUM_DEPTH unlogged here and logs it afterwards, once the
  // named configuration has confirmed it was set.
  params.get("MAXIMUM_DEPTH", spec.max_depth,
             {.desc = "The maximum depth of the ocean.",
              .units = "m",
              .fail_if_missing = true,
              .do_not_log = true});

  params.get("MINIMUM_DEPTH", spec.min_depth,
             {.default_value = 0.0,
              .desc = "The minimum depth of the ocean.",
              .units = "m"});

  if (config != "flat") {
    params.get("EDGE_DEPTH", spec.edge_depth,
               {.default_value = 100.0,
                .desc = "The depth at the edge of one of the named topographies.",
                .units = "m"});
    params.get("TOPOG_SLOPE_SCALE", spec.expdecay,
               {.default_value = 400000.0,
                .desc = "The exponential decay scale used in defining some of "
                        "the named topographies.",
                .units = "m"});
  }

  if (!(spec.max_depth > 0.0)) {
    logger::fatal("initialize_fixed: MAXIMUM_DEPTH must be positive.");
  }
  if (spec.min_depth < 0.0) {
    // MOM6 accepts a negative MINIMUM_DEPTH only with MASKING_DEPTH set, which
    // is not implemented here.
    logger::fatal("initialize_fixed: MINIMUM_DEPTH must not be negative.");
  }

  return config;
}

} // namespace

GridFields initialize_fixed(const Domain &domain, RuntimeParams &params) {

  params.doc_module("MOM_grid_init", "");

  std::string config;
  params.get("GRID_CONFIG", config,
             {.desc = "A character string that determines the method for defining the horizontal "
                      "grid. Current options are:\n"
                      "\t mosaic - read the grid from a mosaic (supergrid)\n"
                      "\t\t file set by GRID_FILE.\n"
                      "\t cartesian - use a (flat) Cartesian grid.\n"
                      "\t spherical - use a simple spherical grid.\n"
                      "\t mercator - use a Mercator spherical grid.",
              .fail_if_missing = true});

  GridSpec spec;
  GridFields fields;
  if (config == "spherical") {
    read_spherical_grid_params(params, spec);
    fields = spherical_grid_fields(domain, spec);
  } else if (config == "mosaic" || config == "cartesian" || config == "mercator") {
    // defer: the mosaic (file-based), cartesian, and mercator grid
    //        configurations.
    logger::fatal("initialize_fixed: GRID_CONFIG \"", config,
                  "\" is not implemented yet.");
  } else if (config == "file") {
    // Retired in MOM6 itself; carry its message.
    logger::fatal("initialize_fixed: GRID_CONFIG \"file\" is no longer a supported "
                  "option. Use a mosaic file (\"mosaic\") or one of the analytic "
                  "forms instead.");
  } else {
    logger::fatal("initialize_fixed: Unrecognized grid configuration \"", config, "\".");
  }

  // Topography and the land/sea masks, between the metrics and the rotation,
  // matching MOM6's MOM_initialize_fixed order.
  TopoSpec topo_spec;
  const std::string topo_config = read_topography_params(params, topo_spec);
  fields.max_depth = topo_spec.max_depth;
  fields.bathyT = named_topography(domain, topo_config, spec, topo_spec,
                                   fields.geoLonT, fields.geoLatT);

  // MOM6 logs MAXIMUM_DEPTH here, after the named configuration has used it
  // (log_param, following the unlogged read above). RuntimeParams has no
  // log-only entry point, so this is a second read of the same key, which
  // produces the same documentation line.
  params.get("MAXIMUM_DEPTH", topo_spec.max_depth,
             {.desc = "The maximum depth of the ocean.",
              .units = "m",
              .fail_if_missing = true});

  initialize_masks(domain, topo_spec, fields.bathyT, fields);
  set_derived_metrics(domain, fields);

  read_rotation_params(params, spec);
  fields.CoriolisBu = planetary_rotation(domain, spec, fields.geoLatBu);

  return fields;
}

} // namespace MOM
