#pragma once
/// @file MOM_debug_dump.h
/// @brief An ASCII field dump for parity checking against legacy MOM6.
///        Temporary scaffolding: it exists because protoMOMxx has no I/O layer
///        yet, and it retires when diagnostics arrive through TIM.

#include <string>

#include <AMReX_MultiFab.H>

namespace MOM::debug {

/// @brief Write the valid region of a field to "<name>.rank<n>.txt", one line
/// per point, with the value in both decimal and hexadecimal floating point.
/// The hexadecimal form is exact, so a diff against the same dump from another
/// build is a bit-for-bit comparison.
/// @param mf The field to dump.
/// @param name The output file's base name.
/// @param comp The component to dump.
void dump_field(const amrex::MultiFab &mf, const std::string &name, int comp = 0);

/// @brief Log a field's mean, minimum, maximum and bit-count checksum over
/// the same points MOM6's chksum routines use, so the two logs can be diffed
/// stage by stage.
///
/// MOM6 reports a velocity field over its non-symmetric range: I = isc..iec
/// for a u point, which is AMReX face index 1..NI, so the western boundary
/// face is left out. The ranges here reproduce that, and the checksum is
/// MOM6's: the sum over the reported points of the population count of the
/// absolute value's bit pattern, modulo 1e9.
/// @param mf The field to report.
/// @param name The label to print.
/// @param comp The component to report.
void report_field(const amrex::MultiFab &mf, const std::string &name, int comp = 0);

/// @brief Whether the stage-by-stage reporting is on.
/// @return True if reporting is enabled.
bool reporting();

/// @brief Turn the stage-by-stage reporting on or off.
/// @param on Whether to report.
void set_reporting(bool on);

} // namespace MOM::debug
