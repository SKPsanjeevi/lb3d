#include "lbe.h"

!> Pointlike tracer particles, independent from MD modules and
!> therefore able to efficiently coexist with MD particles in a single
!> simulation
module lbe_tracer_globals_module
#ifdef TRACER
    use lbe_globals_module, only: rk

    implicit none
    public

    integer,parameter :: tracer_input_file_unit=20

    !> a type for tracer particles:
    type tracer_particle_type
       real(kind=rk),dimension(3) :: x !< position
       real(kind=rk),dimension(3) :: v !< velocity
       integer :: uid                  !< unique tracer id
       !> used to distinquish different groups of tracers
       integer :: kind
    end type tracer_particle_type

    !> contains what is needed for \c exchange()
    integer,save :: tracer_exch_mpitype

    integer,save :: sproc(3,4) !< rankc of cpu to send to for each dim and dir
    integer,save :: rproc(3,4) !< rankc of cpu to recv from in each dim and dir

    !> for building send mpitype
    integer,allocatable,save :: sidxs(:)

    !> for \c mpi_type_indexed() in \c exchange()
    integer,allocatable,save :: slens(:)

    !> tracers local on each processor:
    type(tracer_particle_type),allocatable,save :: T(:)

    !> particle recv buffer for \c exchange() :
    type(tracer_particle_type),allocatable,save :: tbuf(:)

    integer,save :: nfmax       !< size of \c pbuf(:)
    !> maximum number of owned particles, halo'ed particles start at
    !> \c P(npmax+1)
    integer,save :: npmax

    !> linked list of local atoms (last one -> \c npmax+1)
    integer,allocatable,save :: list(:)

    integer,save :: atompnt !< pointer to 1st atom in my list
    integer,save :: freepnt !< pointer to 1st free space in list (last one -> 0)
    !> \name indices referencing timers
    !> \{
    integer,save :: ti_tracer_comm,ti_tracer_dump,ti_tracer_move
    !> \}
    integer,save :: nlocal      !< # of atoms I currently own
#endif
end module lbe_tracer_globals_module
