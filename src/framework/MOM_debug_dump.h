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

} // namespace MOM::debug
