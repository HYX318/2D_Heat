/**
 * @file jacobi_solver.cpp
 * @brief Jacobi solver implementation
 *
 * This file contains explicit template instantiations for the Jacobi solver.
 * The main implementation is in jacobi_solver.hpp.
 */

#include "jacobi_solver.hpp"

// Explicit template instantiations
// This allows the template code to be compiled into a library

// Serial version
template class JacobiSolver<false>;

// Parallel version
template class JacobiSolver<true>;
