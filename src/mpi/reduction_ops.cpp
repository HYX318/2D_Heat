#include "reduction_ops.hpp"
#include "mpi_context.hpp"
#include <limits>
#include <string>

// Constructor implementations

ReductionOps::ReductionOps() : comm_(MPI_COMM_WORLD) {
}

ReductionOps::ReductionOps(MPI_Comm comm) : comm_(comm) {
}

ReductionOps::ReductionOps(const MPIContext& ctx) : comm_(MPI_COMM_WORLD) {
    // The MPIContext uses MPI_COMM_WORLD internally
    // If we need to support custom communicators in MPIContext,
    // we would need to add a method to access it
}

// Barrier implementation

void ReductionOps::barrier() const {
    int result = MPI_Barrier(comm_);
    if (result != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Barrier failed with error code: " + std::to_string(result));
    }
}

// get_mpi_op implementation

MPI_Op ReductionOps::get_mpi_op(ReductionOp op) {
    switch (op) {
        case ReductionOp::SUM:
            return MPI_SUM;
        case ReductionOp::MAX:
            return MPI_MAX;
        case ReductionOp::MIN:
            return MPI_MIN;
        case ReductionOp::PROD:
            return MPI_PROD;
        case ReductionOp::LAND:
            return MPI_LAND;
        case ReductionOp::LOR:
            return MPI_LOR;
        case ReductionOp::LXOR:
            return MPI_LXOR;
        case ReductionOp::BAND:
            return MPI_BAND;
        case ReductionOp::BOR:
            return MPI_BOR;
        case ReductionOp::BXOR:
            return MPI_BXOR;
        default:
            throw std::invalid_argument("Unsupported reduction operation");
    }
}

// Template method implementations
// These are included in the header file for template instantiation,
// but we need explicit instantiations for supported types

// Explicit instantiations for reduce
template int ReductionOps::reduce<int>(int, ReductionOp, int) const;
template long ReductionOps::reduce<long>(long, ReductionOp, int) const;
template long long ReductionOps::reduce<long long>(long long, ReductionOp, int) const;
template unsigned ReductionOps::reduce<unsigned>(unsigned, ReductionOp, int) const;
template unsigned long ReductionOps::reduce<unsigned long>(unsigned long, ReductionOp, int) const;
template unsigned long long ReductionOps::reduce<unsigned long long>(unsigned long long, ReductionOp, int) const;
template float ReductionOps::reduce<float>(float, ReductionOp, int) const;
template double ReductionOps::reduce<double>(double, ReductionOp, int) const;
template long double ReductionOps::reduce<long double>(long double, ReductionOp, int) const;

// Explicit instantiations for allreduce
template int ReductionOps::allreduce<int>(int, ReductionOp) const;
template long ReductionOps::allreduce<long>(long, ReductionOp) const;
template long long ReductionOps::allreduce<long long>(long long, ReductionOp) const;
template unsigned ReductionOps::allreduce<unsigned>(unsigned, ReductionOp) const;
template unsigned long ReductionOps::allreduce<unsigned long>(unsigned long, ReductionOp) const;
template unsigned long long ReductionOps::allreduce<unsigned long long>(unsigned long long, ReductionOp) const;
template float ReductionOps::allreduce<float>(float, ReductionOp) const;
template double ReductionOps::allreduce<double>(double, ReductionOp) const;
template long double ReductionOps::allreduce<long double>(long double, ReductionOp) const;

// Explicit instantiations for scan
template int ReductionOps::scan<int>(int, ReductionOp) const;
template long ReductionOps::scan<long>(long, ReductionOp) const;
template long long ReductionOps::scan<long long>(long long, ReductionOp) const;
template unsigned ReductionOps::scan<unsigned>(unsigned, ReductionOp) const;
template unsigned long ReductionOps::scan<unsigned long>(unsigned long, ReductionOp) const;
template unsigned long long ReductionOps::scan<unsigned long long>(unsigned long long, ReductionOp) const;
template float ReductionOps::scan<float>(float, ReductionOp) const;
template double ReductionOps::scan<double>(double, ReductionOp) const;
template long double ReductionOps::scan<long double>(long double, ReductionOp) const;

// Explicit instantiations for exscan
template int ReductionOps::exscan<int>(int, ReductionOp) const;
template long ReductionOps::exscan<long>(long, ReductionOp) const;
template long long ReductionOps::exscan<long long>(long long, ReductionOp) const;
template unsigned ReductionOps::exscan<unsigned>(unsigned, ReductionOp) const;
template unsigned long ReductionOps::exscan<unsigned long>(unsigned long, ReductionOp) const;
template unsigned long long ReductionOps::exscan<unsigned long long>(unsigned long long, ReductionOp) const;
template float ReductionOps::exscan<float>(float, ReductionOp) const;
template double ReductionOps::exscan<double>(double, ReductionOp) const;
template long double ReductionOps::exscan<long double>(long double, ReductionOp) const;

// Explicit instantiations for reduce_array
template std::vector<int> ReductionOps::reduce_array<int>(const std::vector<int>&, ReductionOp, int) const;
template std::vector<long> ReductionOps::reduce_array<long>(const std::vector<long>&, ReductionOp, int) const;
template std::vector<long long> ReductionOps::reduce_array<long long>(const std::vector<long long>&, ReductionOp, int) const;
template std::vector<unsigned> ReductionOps::reduce_array<unsigned>(const std::vector<unsigned>&, ReductionOp, int) const;
template std::vector<unsigned long> ReductionOps::reduce_array<unsigned long>(const std::vector<unsigned long>&, ReductionOp, int) const;
template std::vector<unsigned long long> ReductionOps::reduce_array<unsigned long long>(const std::vector<unsigned long long>&, ReductionOp, int) const;
template std::vector<float> ReductionOps::reduce_array<float>(const std::vector<float>&, ReductionOp, int) const;
template std::vector<double> ReductionOps::reduce_array<double>(const std::vector<double>&, ReductionOp, int) const;
template std::vector<long double> ReductionOps::reduce_array<long double>(const std::vector<long double>&, ReductionOp, int) const;

// Explicit instantiations for allreduce_array
template std::vector<int> ReductionOps::allreduce_array<int>(const std::vector<int>&, ReductionOp) const;
template std::vector<long> ReductionOps::allreduce_array<long>(const std::vector<long>&, ReductionOp) const;
template std::vector<long long> ReductionOps::allreduce_array<long long>(const std::vector<long long>&, ReductionOp) const;
template std::vector<unsigned> ReductionOps::allreduce_array<unsigned>(const std::vector<unsigned>&, ReductionOp) const;
template std::vector<unsigned long> ReductionOps::allreduce_array<unsigned long>(const std::vector<unsigned long>&, ReductionOp) const;
template std::vector<unsigned long long> ReductionOps::allreduce_array<unsigned long long>(const std::vector<unsigned long long>&, ReductionOp) const;
template std::vector<float> ReductionOps::allreduce_array<float>(const std::vector<float>&, ReductionOp) const;
template std::vector<double> ReductionOps::allreduce_array<double>(const std::vector<double>&, ReductionOp) const;
template std::vector<long double> ReductionOps::allreduce_array<long double>(const std::vector<long double>&, ReductionOp) const;

// Explicit instantiations for maxloc
template MaxLocResult<int> ReductionOps::maxloc<int>(int) const;
template MaxLocResult<long> ReductionOps::maxloc<long>(long) const;
template MaxLocResult<long long> ReductionOps::maxloc<long long>(long long) const;
template MaxLocResult<unsigned> ReductionOps::maxloc<unsigned>(unsigned) const;
template MaxLocResult<unsigned long> ReductionOps::maxloc<unsigned long>(unsigned long) const;
template MaxLocResult<unsigned long long> ReductionOps::maxloc<unsigned long long>(unsigned long long) const;
template MaxLocResult<float> ReductionOps::maxloc<float>(float) const;
template MaxLocResult<double> ReductionOps::maxloc<double>(double) const;
template MaxLocResult<long double> ReductionOps::maxloc<long double>(long double) const;

// Explicit instantiations for minloc
template MinLocResult<int> ReductionOps::minloc<int>(int) const;
template MinLocResult<long> ReductionOps::minloc<long>(long) const;
template MinLocResult<long long> ReductionOps::minloc<long long>(long long) const;
template MinLocResult<unsigned> ReductionOps::minloc<unsigned>(unsigned) const;
template MinLocResult<unsigned long> ReductionOps::minloc<unsigned long>(unsigned long) const;
template MinLocResult<unsigned long long> ReductionOps::minloc<unsigned long long>(unsigned long long) const;
template MinLocResult<float> ReductionOps::minloc<float>(float) const;
template MinLocResult<double> ReductionOps::minloc<double>(double) const;
template MinLocResult<long double> ReductionOps::minloc<long double>(long double) const;

// Explicit instantiations for broadcast
template void ReductionOps::broadcast<int>(int&, int) const;
template void ReductionOps::broadcast<long>(long&, int) const;
template void ReductionOps::broadcast<long long>(long long&, int) const;
template void ReductionOps::broadcast<unsigned>(unsigned&, int) const;
template void ReductionOps::broadcast<unsigned long>(unsigned long&, int) const;
template void ReductionOps::broadcast<unsigned long long>(unsigned long long&, int) const;
template void ReductionOps::broadcast<float>(float&, int) const;
template void ReductionOps::broadcast<double>(double&, int) const;
template void ReductionOps::broadcast<long double>(long double&, int) const;

// Explicit instantiations for broadcast_array
template void ReductionOps::broadcast_array<int>(std::vector<int>&, int) const;
template void ReductionOps::broadcast_array<long>(std::vector<long>&, int) const;
template void ReductionOps::broadcast_array<long long>(std::vector<long long>&, int) const;
template void ReductionOps::broadcast_array<unsigned>(std::vector<unsigned>&, int) const;
template void ReductionOps::broadcast_array<unsigned long>(std::vector<unsigned long>&, int) const;
template void ReductionOps::broadcast_array<unsigned long long>(std::vector<unsigned long long>&, int) const;
template void ReductionOps::broadcast_array<float>(std::vector<float>&, int) const;
template void ReductionOps::broadcast_array<double>(std::vector<double>&, int) const;
template void ReductionOps::broadcast_array<long double>(std::vector<long double>&, int) const;
