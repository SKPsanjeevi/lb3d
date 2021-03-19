#include "lbe.h"

!> functionality to dynamically add and remove MD particles during a
!> simulation
module lbe_md_dynamic_module
#ifdef MD
    use lbe_globals_module, only: myrankc
    use lbe_log_module
    use lbe_md_globals_module
    use lbe_md_helper_module, only: count_particles_all,error_md,log_msg_md&
         &,trigger_list_update
    use lbe_md_memory_module, only: boost_npmax
    use lbe_parallel_module, only: check_allocate,comm_cart,nprocs
    use lbe_parms_module, only: nt
    use map_module, only: Mii_map
    implicit none
    include 'mpif.h'
    private

    !> determines how to obtain uid numbers for newly inserted
    !> particles (see \c get_free_puid_by_rank(), \c
    !> get_free_puid_by_time_step(), and \c set_free_puid_by_master()
    !> for explanations)
    character(len=32),save,public :: new_uids='rank'

    !> lowest particle uid which does not belong to a particle that
    !> was placed during initialization
    integer,save,public :: first_inserted_puid

    !> Assumed maximum number of slaves one master can spawn during a
    !> simulation (used for new_uids=='master')
    integer,parameter,public :: max_slaves_per_master = 10000

    public delete_particle,dissolve,insert,insert_particle,remove&
         &,set_dissolution_limits_z,setup_dynamic,to_be_deleted,uid_of_deletee

    !> Every particle that leaves this z-position range will dissolve
    !> (range is closed on the lower side and open on the upper one)
    real(kind=rk),save :: dissolution_limits_z(2)=&
         &(/-huge(1.0_rk),huge(1.0_rk)/)

    !> which bit in a particle uid triggers dissolution?
    integer,parameter :: deletion_bit=30

    !> next free particle uid to be assigned by the local process
    !> (used by \c get_free_puid_by_rank())
    integer,save :: next_free_puid

    !> local number of particles in \c insert_list
    integer,save :: n_insert=0

    !> buffer for particles to be inserted locally within the current
    !> time step
    type(md_particle_type),allocatable,save :: insert_list(:)

contains

    !> marks a particle for dissolution within one time step
    !>
    !> \param[in,out] p particle
    subroutine delete_particle(p)
        type(md_particle_type),intent(inout) :: p

        p%uid = ibset(p%uid,deletion_bit)
        call trigger_list_update
    end subroutine delete_particle

    !> Trigger deletion for all particles that exceed the z-position
    !> range given by \c dissolution_limits_z
    !>
    !> This works also on halo'ed particles, so dissolution does not
    !> need to be communicated (for example to update the rock state
    !> in case of \c interaction=='ladd').
    subroutine dissolve
        integer :: i,ii

        i = atompnt
        do ii = 1,nlocal+nother
           ! for the moment dissolve only particles that were inserted
           ! after the initialization
           if (P(i)%uid>=first_inserted_puid.and.&
                &(P(i)%x(3)<dissolution_limits_z(1)&
                &.or.P(i)%x(3)>=dissolution_limits_z(2))) then
              call delete_particle(P(i))
           end if

           if (ii<=nlocal) then
              i = list(i)
           else
              i = i+1
           endif
        end do
    end subroutine dissolve

    !> get next free particle uid according to \c new_uids=='rank'
    !>
    !> An interleaved scheme is used to distribute all possible
    !> integers on the different ranks. Each rank assignes his uid
    !> numbers in ascending order. This does not (yet) support
    !> checkpointing.
    integer function get_free_puid_by_rank()
        get_free_puid_by_rank = next_free_puid
        next_free_puid = next_free_puid+nprocs
    end function get_free_puid_by_rank

    !> get next free particle uid according to \c new_uids=='time_step'
    !>
    !> Particle uids are assigned based on the current time step. This
    !> allows the insertion of only one particle per time step within
    !> the whole run. The program will stop if more particles would be
    !> created.
    integer function get_free_puid_by_time_step()
        get_free_puid_by_time_step = first_inserted_puid+nt
    end function get_free_puid_by_time_step

    !> actually inserts the particles passed by \c insert_particle()
    !> into the local \c P
    !>
    !> \param[in] nunu determines how new particle uids are generated
    !> (optional, defaults to the the value in \c new_uids).
    !>
    !> \c nunu can be 'keep_uid' (in which case the uids stay as they
    !> were passed to \c insert_particle()) or any of the values
    !> acceptable for \c new_uids.
    subroutine insert(nunu)
        character(len=*),intent(in),optional :: nunu
        character(len=32) :: nu
        integer i,ierror,itmp,n_insert_global

        if (present(nunu)) then
           nu = nunu
        else
           nu = new_uids
        end if

        select case(nu)
        case ('keep_uid')
           ! just do nothing and keep the present uids

        case ('master')
           do i=1,n_insert
              call set_free_puid_by_master(insert_list(i))
           end do

        case ('rank')
           do i=1,n_insert
              insert_list(i)%uid = get_free_puid_by_rank()
           end do

        case ('time_step')
           call MPI_Reduce(n_insert,n_insert_global,1,MPI_INTEGER,MPI_SUM,0&
                &,comm_cart,ierror)
           if (myrankc==0.and.n_insert_global>1) &
                &call error_md('Inserting more than one particle per time step'&
                &//' is not supported with new_uids=''time_step''.')
           if (n_insert==1) insert_list(1)%uid = get_free_puid_by_time_step()

        case default
           call error_md('unknown policy for new particle uids: new_uids="'&
                &//nu//'"')
        end select

        ! can't do this before MPI_Reduce in case 'time_step'
        if (n_insert==0) return

        call trigger_list_update
        if (nlocal+n_insert>npmax) CALL boost_npmax(nlocal+n_insert)
        nlocal = nlocal+n_insert

        do i=1,n_insert
           itmp = atompnt
           atompnt = freepnt
           freepnt = list(freepnt)
           list(atompnt) = itmp

           P(atompnt) = insert_list(i)
        end do

        n_insert = 0
    end subroutine insert

    !> locally inserts a new particle within one time step
    !>
    !> \param[in] p particle
    !>
    !> In general, the actual uid of the new particle will be set
    !> automatically, so p%uid will be ignored.
    subroutine insert_particle(p)
        type(md_particle_type),intent(in) :: p
        type(md_particle_type),allocatable :: tmp(:)
        integer stat

        if (.not.allocated(insert_list)) then
           allocate (insert_list(1),stat=stat)
           call check_allocate(stat,'insert_particle(): insert_list')
        end if
        if (size(insert_list)<n_insert+1) then
           allocate(tmp(n_insert),stat=stat)
           call check_allocate(stat,'insert_particle(): tmp')
           tmp = insert_list(1:n_insert)
           deallocate(insert_list)
           allocate(insert_list(2*n_insert),stat=stat)
           call check_allocate(stat,'insert_particle(): insert_list')
           insert_list(1:n_insert) = tmp
           deallocate(tmp)
        end if

        n_insert = n_insert+1
        insert_list(n_insert) = p
    end subroutine insert_particle

    !> throw away particles to be deleted, update local list
    subroutine remove
        integer :: i,ii,iprev,itmp,ndelete

        ndelete = 0
        iprev = 0
        i = atompnt

        do ii = 1,nlocal
           if (to_be_deleted(P(i))) then
              ndelete = ndelete+1
              if (iprev==0) then
                 atompnt = list(i)
              else
                 list(iprev) = list(i)
              endif
              itmp = list(i)
              list(i) = freepnt
              freepnt = i
              i = itmp
           else
              iprev = i
              i = list(i)
           endif
        enddo
        nlocal = nlocal-ndelete
    end subroutine remove

    !> initializaes \c dissolution_limits_z
    !>
    !> \param[in] lo lower limit
    !>
    !> \param[in] hi upper limit
    !>
    !> Dissolution is triggered as soon as a particle center position
    !> gets lower than \c lo or equal to \c hi.
    subroutine set_dissolution_limits_z(lo,hi)
        real(kind=rk),intent(in) :: lo,hi

        dissolution_limits_z = (/lo,hi/)

        write (unit=msgstr,fmt='("Dissolution of non-periodic particles in '&
             &//'[*:",F9.3,"[ and [",F9.3,":*[")') dissolution_limits_z
        call log_msg_md(msgstr)
    end subroutine set_dissolution_limits_z

    !> sets a particle's uid to the next free particle uid according
    !> to \c new_uids=='master'
    !>
    !> \param[in] new_p the particle to insert
    !>
    !> The master particle of \c new_p is obtained from \c
    !> new_p%master. The whole range of possible uids is distributed
    !> on all possible master particles assuming that during one
    !> simulation no single master will spawn more than \c
    !> max_slaves_per_master slaves. For each master, his uids are
    !> assigned in ascending order.
    subroutine set_free_puid_by_master(new_p)
        type(md_particle_type),intent(inout) :: new_p
        integer mi

        mi = Mii_map(uid2i,new_p%master)

        new_p%uid = -P(mi)%master
        if (new_p%uid>=first_inserted_puid+P(mi)%uid*max_slaves_per_master) then
           write (unit=6,fmt='("master: uid=",I10,",master=",I10)') &
                &P(mi)%uid,P(mi)%master
           call error_md('Reached max_slaves_per_master for this simulation!')
        end if

        P(mi)%master = P(mi)%master-1
    end subroutine set_free_puid_by_master

    !> initialize this module
    subroutine setup_dynamic

        ! initialize first_inserted_puid which. This will stay
        ! constant all the time.
        call count_particles_all(first_inserted_puid)
        first_inserted_puid = first_inserted_puid+1

        write (unit=msgstr,fmt='("first_inserted_puid = ",I0)') &
             &first_inserted_puid
        call log_msg_md(msgstr)

        call setup_next_free_puid
    end subroutine setup_dynamic

    !> initialize \c next_free_puid which is used by \c new_uids=='rank'
    subroutine setup_next_free_puid
        integer :: gmaxu,i,ierror,ii,maxu

        maxu = -1
        i = atompnt
        do ii=1,nlocal
           maxu = max(P(i)%uid,maxu)
           i = list(i)
        end do

        call mpi_allreduce(maxu,gmaxu,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD&
             &,ierror)

        next_free_puid = gmaxu+1+myrankc
    end subroutine setup_next_free_puid

    !> returns whether a particle was marked for dissolution
    !>
    !> \param[in] p the particle
    pure logical function to_be_deleted(p)
        type(md_particle_type),intent(in) :: p

        to_be_deleted = btest(p%uid,deletion_bit)
    end function to_be_deleted

    !> returns the uid of a particle that was marked for dissolution
    !>
    !> \param[in] p the particle
    !>
    !> In the current implementation, one bit of the particle uid is
    !> used to mark particles for deletion, so if it is marked, \c uid
    !> does not contain the plain uid anymore.
    pure integer function uid_of_deletee(p)
        type(md_particle_type),intent(in) :: p

        uid_of_deletee = ibclr(p%uid,deletion_bit)
    end function uid_of_deletee

#endif
end module lbe_md_dynamic_module
