#pragma once
/// @file MOM_loop_boxes.h
/// @brief The iteration boxes the dynamics kernels use, named after the MOM6
///        loop bounds they correspond to.
///
/// MOM6 writes its loops over index ranges taken from the grid
/// (is..ie, Isq..Ieq, and so on). Under symmetric memory those ranges map
/// exactly onto AMReX boxes: a rank's computational domain is its valid box,
/// and MOM6's Isq..Ieq velocity range is the valid box's surrounding x faces
/// (MOM6's I = i-1 in AMReX face numbering). Naming the boxes here keeps the
/// correspondence visible at each kernel and keeps the halo reasoning in one
/// place.

#include <AMReX_Box.H>

namespace MOM::loops {

/// @brief The h points of the computational domain (MOM6's is..ie, js..je).
/// @param valid The valid (cell-centered) box of the current FAB.
/// @return The same box.
inline amrex::Box h_points(const amrex::Box &valid) { return valid; }

/// @brief The h points of the computational domain grown by n halo rings
/// (MOM6's Isq-n+1..Ieq+n, i.e. is-n..ie+n).
/// @param valid The valid (cell-centered) box of the current FAB.
/// @param n The number of halo rings.
/// @return The grown cell-centered box.
inline amrex::Box h_points_grown(const amrex::Box &valid, const int n) {
  return amrex::grow(valid, amrex::IntVect(n, n, 0));
}

/// @brief The u points of the computational domain (MOM6's Isq..Ieq, js..je).
/// @param valid The valid (cell-centered) box of the current FAB.
/// @param n Halo rings to add in each direction.
/// @return The x-face box.
inline amrex::Box u_points(const amrex::Box &valid, const int n = 0) {
  return amrex::surroundingNodes(amrex::grow(valid, amrex::IntVect(n, n, 0)), 0);
}

/// @brief The v points of the computational domain (MOM6's is..ie, Jsq..Jeq).
/// @param valid The valid (cell-centered) box of the current FAB.
/// @param n Halo rings to add in each direction.
/// @return The y-face box.
inline amrex::Box v_points(const amrex::Box &valid, const int n = 0) {
  return amrex::surroundingNodes(amrex::grow(valid, amrex::IntVect(n, n, 0)), 1);
}

/// @brief The q points of the computational domain (MOM6's Isq..Ieq, Jsq..Jeq).
/// @param valid The valid (cell-centered) box of the current FAB.
/// @param n Halo rings to add in each direction.
/// @return The nodal box.
inline amrex::Box q_points(const amrex::Box &valid, const int n = 0) {
  const amrex::Box grown = amrex::grow(valid, amrex::IntVect(n, n, 0));
  return amrex::surroundingNodes(amrex::surroundingNodes(grown, 0), 1);
}

/// @brief Collapse a box's vertical range to the single level k = 0, for
/// loops over 2-D fields driven by an MFIter over a 3-D field.
/// @param bx The box to collapse.
/// @return The same box with k running over 0 only.
inline amrex::Box flat(amrex::Box bx) {
  bx.setSmall(2, 0);
  bx.setBig(2, 0);
  return bx;
}

/// @brief Set a box's vertical range to the model's layers.
/// @param bx The box to set.
/// @param nk The number of layers.
/// @return The same box with k running over 0 .. nk-1.
inline amrex::Box layers(amrex::Box bx, const int nk) {
  bx.setSmall(2, 0);
  bx.setBig(2, nk - 1);
  return bx;
}

} // namespace MOM::loops
