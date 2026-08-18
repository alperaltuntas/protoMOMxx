# src/types

The passive model types: the objects the model carries and nearly every part of it reads, i.e., the grids, the prognostic state, and the forcing. The construction precursors (e.g., `GridFields`) live here too, since the types' constructors take them.

This directory is a layer below both `src/initialization` and `src/parameterizations`, because both need these types: the initialization code computes their contents, the dynamics reads them, and `src/core`'s stepper calls the dynamics. Keeping the types below everyone keeps the per-directory libraries one-way: framework <- types <- {initialization, parameterizations, diagnostics} <- core.

Two rules keep this directory what it is:

- A type belongs here if it is model-carried data that multiple layers read. Infrastructure wrappers (e.g., `Domain`) stay in `config_src/infra`, and services (parser, logger, IO) stay in `src/framework`.
- Nothing here reads runtime parameters or computes field values. That is `src/initialization`'s job; the constructors here validate what they are handed and take ownership.

MOM6 keeps most of these files in its `src/core` (e.g., `MOM_grid.F90`, `MOM_verticalGrid.F90`, `MOM_forcing_type.F90`). The move down is deliberate, so as to keep the library graph acyclic, and the MOM6 file names are kept.
