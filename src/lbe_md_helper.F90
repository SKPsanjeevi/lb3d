#include "lbe.h"

!> helper routines for md
module lbe_md_helper_module
#ifdef MD

    use lbe_helper_module, only: local_coordinates,present_and_true
    use lbe_globals_module,only: c,chunksize,cx,cy,cz,g,halo_extent,lp_sur,n_lp&
         &,nnonrest,opp_lp_sur,pi,tsize,myrankc
    use lbe_log_module
    use lbe_md_fluid_ladd_parms_module, only: initial_radii,R_max,R_orth,R_para
    use lbe_md_globals_module
    use lbe_parallel_module, only: ccoords,check_allocate,checkmpi,comm_cart&
         &,find_topology,nprocs,start,tnx,tny,tnz
    use lbe_parms_module, only: amass_r, amass_b, amass_s, folder, cpfolder, &
         gr_out_file,tau_b,tau_r,tau_s,g_accn_min,g_accn_max,g_accn_min_x&
         &,g_accn_max_x,g_accn_min_y,g_accn_max_y,g_accn, g_accn_x, g_accn_y&
         &,nx,ny,nz
    use lbe_io_helper_module, only: restore_fmt_p
    use lbe_types_module, only: lbe_site
    use lbe_collision_module, only:  md_calculate_sc_forces
    use lbe_parms_module, only: chk_uid
    implicit none
    include 'mpif.h'
    private

    public build_collect_mpitypes,build_ft_buffer_mpitype&
         &,build_particle_mpitype,count_particles,count_particles_all,error_md&
         &,fluid_velocity,fluid_velocity_and_viscosity,inside_rock&
         &,list_update_add_displacement,log_msg_md,log_msg_md_ws,log_msg_md_hdr&
         &,maximum_particle_radius,md_make_filename_particles,number_density&
         &,orientation,particle_half_axes,particle_inertia&
         &,particle_inverse_half_axes,particle_mass,particle_radii,quaternion&
         &,scatter_rock_state_halo,space_to_body_matrix,sum_dx_all&
         &,trigger_list_update,warning,check_quaternion

    character(len=4), parameter :: md_log_prefix = "[MD]"

contains

    !> consider an additional increase of particle displacement
    !> effective for determining whether a neighbor list update is
    !> required.
    !>
    !> \param[in] d additional displacement
    !>
    !> It can be called multiple times. The effect is cumulative. It
    !> must be called by all ranks with the same argument.
    subroutine list_update_add_displacement(d)
        real(kind=rk),intent(in) :: d

        list_update_displ_add = list_update_displ_add + d
    end subroutine list_update_add_displacement

!> Log message with MD prefix and surrounding whitespace.
subroutine log_msg_md_ws(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg
  logical, optional, intent(in) :: forall

  if ( present(forall) ) then
    call log_msg("", forall)
    call log_msg(md_log_prefix // " " // trim(msg), forall)
    call log_msg("", forall)
  else
    call log_msg("")
    call log_msg(md_log_prefix // " " // trim(msg))
    call log_msg("")
  end if

end subroutine log_msg_md_ws

!> Log message with MD prefix.
subroutine log_msg_md(msg, forall, error)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.
  logical, optional, intent(in) :: forall !< If supplied and set to true, all ranks will log.
  logical, optional, intent(in) :: error  !< If supplied and set to true, log an error.

  if ( present(forall) ) then
    if ( present(error) ) then
      call log_msg(md_log_prefix // " " // trim(msg), forall = forall, error = error)
    else
      call log_msg(md_log_prefix // " " // trim(msg), forall = forall)
    end if
  else
    if ( present(error) ) then
      call log_msg(md_log_prefix // " " // trim(msg), error = error)
    else
      call log_msg(md_log_prefix // " " // trim(msg))
    end if
  end if

end subroutine log_msg_md

!> Log header message with MD prefix.
subroutine log_msg_md_hdr(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.
  logical, optional, intent(in) :: forall !< If supplied and set to true, all ranks will log.

  if ( present(forall) ) then
    call log_msg_hdr(md_log_prefix // " " // trim(msg), forall)
  else
    call log_msg_hdr(md_log_prefix // " " // trim(msg))
  end if

end subroutine log_msg_md_hdr

!> Log error message with MD prefix and exit.
subroutine error_md(msg)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.

  integer :: mpierror

  call log_msg_md("ERROR: "//trim(msg), forall = .true., error = .false.)
  call log_msg_md("ERROR: "//trim(msg), forall = .true., error = .true. )
  call MPI_Abort(MPI_COMM_WORLD, -1, mpierror)
  stop

end subroutine error_md

  !> build two corresponding mpi data types for sending and receiving
  !> in a \c collect() operation that cover a matching selection of
  !> elements of \c md_particle_type and \c ft_buffer_type controlled
  !> by the optional arguments
  !>
  !> \param[out] pmt custom mpi type for \c md_particle_type object
  !>
  !> \param[out] fmt custom mpi type for \c ft_buffer_type object
  !>
  !> \param[in] f include \c f ?
  !>
  !> \param[in] t include \c t ?
  !>
  !> \param[in] ft_fluid include \c f_fluid and \c t_fluid ?
  !>
  !> \param[in] e_pot include \c e_pot ?
#ifdef PARTICLESTRESS
  !>
  !> \param[in] tau include \c tau ?
#endif
#ifdef LADD_SURR_RHOF
  !>
  !> \param[in] rhofn_acc include \c rhof_acc and \c n_acc ?
#endif
  subroutine build_collect_mpitypes(pmt,fmt,f,t,ft_fluid,e_pot&
#ifdef FORCECOMPONENT
         &,f_normal&
         &,f_tangent&
         &,f_n&
         &,f_t&
#endif
#ifdef PARTICLESTRESS
         &,tau&
#endif
#ifdef LADD_SURR_RHOF
         &,rhofn_acc&
#endif
!#ifdef MD_MAG
         &,ftmi&
!#endif
         &,freeze)
      integer,intent(out) :: pmt,fmt
      logical,intent(in),optional :: f,t,ft_fluid,e_pot
#ifdef FORCECOMPONENT
      logical,intent(in),optional :: f_normal
      logical,intent(in),optional :: f_tangent
      logical,intent(in),optional :: f_n
      logical,intent(in),optional :: f_t
#endif
#ifdef PARTICLESTRESS
      logical,intent(in),optional :: tau
#endif
#ifdef LADD_SURR_RHOF
      logical,intent(in),optional :: rhofn_acc
#endif
!#ifdef MD_MAG
      logical,intent(in),optional :: ftmi
!#endif
      logical,intent(in),optional :: freeze

      call build_particle_mpitype(pmt,f=f,t=t,ft_fluid=ft_fluid,e_pot=e_pot&
#ifdef FORCECOMPONENT
           &,f_normal=f_normal&
           &,f_tangent=f_tangent&
           &,f_n=f_n&
           &,f_t=f_t&
#endif
#ifdef PARTICLESTRESS
           &,tau=tau&
#endif
#ifdef LADD_SURR_RHOF
           &,rhofn_acc=rhofn_acc&
#endif
!#ifdef MD_MAG
           &,ftmi=ftmi&
!#endif
           &,freeze=freeze)
      call build_ft_buffer_mpitype(fmt,f=f,t=t,ft_fluid=ft_fluid,e_pot=e_pot&
#ifdef FORCECOMPONENT
           &,f_normal=f_normal&
           &,f_tangent=f_tangent&
           &,f_n=f_n&
           &,f_t=f_t&
#endif
#ifdef PARTICLESTRESS
           &,tau=tau&
#endif
#ifdef LADD_SURR_RHOF
           &,rhofn_acc=rhofn_acc&
#endif
!#ifdef MD_MAG
           &,ftmi=ftmi&
!#endif
           &)
  end subroutine build_collect_mpitypes

    !> build a custom mpi datatype representing a selection of the
    !> elements of \c ft_buffer_type controlled by the optional
    !> arguments
    !>
    !> \param[out] fmt custom mpi type
    !>
    !> \param[in] f include \c f ?
    !>
    !> \param[in] t include \c t ?
    !>
    !> \param[in] ft_fluid include \c f_fluid and \c t_fluid ?
    !>
    !> \param[in] e_pot include \c e_pot ?
#ifdef PARTICLESTRESS
    !>
    !> \param[in] tau include \c tau ?
#endif
#ifdef LADD_SURR_RHOF
    !>
    !> \param[in] rhofn_acc include \c rhof_acc and \c n_acc ?
#endif
    subroutine build_ft_buffer_mpitype(fmt,f,t,ft_fluid,e_pot&
#ifdef FORCECOMPONENT
         &,f_normal&
         &,f_tangent&
         &,f_n&
         &,f_t&
#endif
#ifdef PARTICLESTRESS
         &,tau&
#endif
#ifdef LADD_SURR_RHOF
         &,rhofn_acc&
#endif
!#ifdef MD_MAG
         &,ftmi&
!#endif
         &)
        integer,intent(out) :: fmt
        logical,intent(in),optional :: f,t,ft_fluid,e_pot
#ifdef FORCECOMPONENT	     
	logical,intent(in),optional :: f_normal,f_tangent,f_n,f_t
#endif	     
#ifdef PARTICLESTRESS
        logical,intent(in),optional :: tau
#endif
#ifdef LADD_SURR_RHOF
        logical,intent(in),optional :: rhofn_acc
#endif
!#ifdef MD_MAG
        logical,intent(in),optional :: ftmi
!#endif

        integer,parameter :: n_blocks_max=&
#ifdef FORCECOMPONENT	     
             &2+&	
#endif	     
#ifdef PARTICLESTRESS
             &1+&
#endif
#ifdef LADD_SURR_RHOF
             &2+&
#endif
!#ifdef MD_MAG
         &2+&
!#endif
             &8
        integer n_blocks,ierror
        integer lengths(n_blocks_max),types(n_blocks_max)
        integer(kind=MPI_ADDRESS_KIND) :: base,addrs(n_blocks_max),displs(n_blocks_max)
        type(ft_buffer_type) :: sample(2)

        n_blocks = 1
        lengths(n_blocks) = 1   ! start of particle in memory
        types(n_blocks) = MPI_LB
        call mpi_get_address(sample(1),addrs(n_blocks),ierror)

        if (present_and_true(f)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! f
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f(1),addrs(n_blocks),ierror)
        end if

        if (present_and_true(t)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! t
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%t(1),addrs(n_blocks),ierror)
        end if
#ifdef FORCECOMPONENT
        if (present_and_true(f_normal)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_normal 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_normal(1),addrs(n_blocks),ierror)
        end if
        if (present_and_true(f_tangent)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_tangent 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_tangent(1),addrs(n_blocks),ierror)
        end if
        if (present_and_true(f_n)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_n 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_n(1),addrs(n_blocks),ierror)
        end if
        if (present_and_true(f_t)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_t 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_t(1),addrs(n_blocks),ierror)
        end if
#endif	
        if (present_and_true(ft_fluid)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! f_fluid
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_fluid(1),addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! t_fluid
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%t_fluid(1),addrs(n_blocks),ierror)
        end if

#ifdef PARTICLESTRESS
        if (present_and_true(tau)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 9          ! tau
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%tau(1,1),addrs(n_blocks),ierror)
        end if

#endif
        if (present_and_true(e_pot)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1          ! e_pot
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%e_pot,addrs(n_blocks),ierror)
        end if

#ifdef LADD_SURR_RHOF
        if (present_and_true(rhofn_acc)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1          ! rhof_acc
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%rhof_acc,addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 1          ! n_acc
           types(n_blocks) = MPI_INTEGER
           call mpi_get_address(sample(1)%n_acc,addrs(n_blocks),ierror)
        end if
#endif
!#ifdef MD_MAG
!       if (present_and_true(ftmi)) then
!           n_blocks = n_blocks+1
!           lengths(n_blocks) = 3          ! fmi
!           types(n_blocks) = MPI_REAL8
!           call mpi_get_address(sample(1)%fmi(1),addrs(n_blocks),ierror)

!           n_blocks = n_blocks+1
 !          lengths(n_blocks) = 3          ! tmi
!           types(n_blocks) = MPI_REAL8
   !        call mpi_get_address(sample(1)%tmi(1),addrs(n_blocks),ierror)
  !      end if

!#endif
        n_blocks = n_blocks+1
        lengths(n_blocks) = 1   ! next particle
        types(n_blocks) = MPI_UB
        call mpi_get_address(sample(2),addrs(n_blocks),ierror)

        call mpi_get_address(sample(1),base,ierror) ! base address
        displs(1:n_blocks) = addrs(1:n_blocks) - base

        call mpi_type_create_struct(n_blocks,lengths,displs,types,fmt,ierror)
        call mpi_type_commit(fmt,ierror)
    end subroutine build_ft_buffer_mpitype

    !> build a custom mpi datatype representing a selection of the
    !> elements of \c md_particle_type controlled by the optional
    !> arguments
    !>
    !> \param[out] pmt custom mpi type
    !>
    !> \param[in] x include \c x ?
    !>
    !> \param[in] v include \c v ?
    !>
    !> \param[in] f include \c f ?
    !>
    !> \param[in] q include \c q ?
    !>
    !> \param[in] w include \c w ?
    !>
    !> \param[in] t include \c t ?
    !>
    !> \param[in] ft_fluid include \c f_fluid and \c t_fluid ?
    !>
    !> \param[in] ft_fluid_prev include \c f_fluid_prev and \c t_fluid_prev ?
    !>
    !> \param[in] vws_fluid include \c v_fluid and \c ws_fluid ?
    !>
    !> \param[in] vws_fluid_avg include \c v_fluid_avg and \c ws_fluid_avg ?
    !>
    !> \param[in] vws_fluid_acc include \c v_fluid_acc and \c ws_fluid_acc ?
    !>
    !> \param[in] uid include \c uid ?
    !>
    !> \param[in] vqw_new include \c vnew , \c qnew , and \c wnew ?
    !>
    !> \param[in] ws include \c ws ?
    !>
    !> \param[in] e_pot include \c e_pot ?
    !>
    !> \param[in] master include \c master ?
    !>
    !> \param[in] R include \c R_orth and \c R_para ?
#ifdef PARTICLESTRESS
    !>
    !> \param[in] tau include \c tau ?
#endif
#ifdef RWALK
    !>
    !> \param[in] v_r include \c v_r ?
    !>
    !> \param[in] sdx include \c sdx ?
#endif
#ifdef LADD_SURR_RHOF
    !>
    !> \param[in] rhof include \c rhof ?
    !>
    !> \param[in] rhofn_acc include \c rhof_acc and \c n_acc ?
#endif
    subroutine build_particle_mpitype(pmt,x,v,f,q,w,t,ft_fluid,ft_fluid_prev&
#ifdef FORCECOMPONENT    
         &,f_normal,f_tangent&
         &,f_n,f_t&
#endif	 
         &,vws_fluid,vws_fluid_avg,vws_fluid_acc,uid,vqw_new,ws,e_pot,master,R, mag&
#ifdef PARTICLESTRESS
         &,tau&
#endif
#ifdef RWALK
         &,v_r,sdx&
#endif
#ifdef LADD_SURR_RHOF
         &,rhof,rhofn_acc&
#endif
!#ifdef MD_MAG
         &,ftmi&
!#endif
         &,freeze)
        integer,intent(out) :: pmt
        logical,intent(in),optional :: x,v,f,q,w,t,ft_fluid,ft_fluid_prev&
             &,vws_fluid,vws_fluid_avg,vws_fluid_acc,uid,vqw_new,ws,e_pot&
#ifdef FORCECOMPONENT	     
	     &,f_normal,f_tangent&
             &,f_n,f_t&
#endif	     
             &,master,R, mag
#ifdef PARTICLESTRESS
        logical,intent(in),optional :: tau
#endif
#ifdef RWALK
        logical,intent(in),optional :: v_r,sdx
#endif
#ifdef LADD_SURR_RHOF
        logical,intent(in),optional :: rhof,rhofn_acc
#endif
!#ifdef MD_MAG
        logical,intent(in),optional :: ftmi
!#endif
        logical,intent(in),optional :: freeze  ! for friction
        integer,parameter :: n_blocks_max=&
#ifdef PARTICLESTRESS
             &1+&
#endif
#ifdef FORCECOMPONENT
             &2+&
#endif
#ifdef RWALK
             &2+&
#endif
#ifdef LADD_SURR_RHOF
             &3+&
#endif
!#ifdef MD_MAG
             &2+&
!#endif
             &26+&
             &2+&  ! for magdisperisty 
             &1  ! for freeze in friction module
        integer n_blocks,ierror
        integer lengths(n_blocks_max),types(n_blocks_max)
        integer(kind=MPI_ADDRESS_KIND) :: base,addrs(n_blocks_max)&
             &,displs(n_blocks_max)
        type(md_particle_type) :: sample(2)

        n_blocks = 1
        lengths(n_blocks) = 1   ! start of particle in memory
        types(n_blocks) = MPI_LB
        call mpi_get_address(sample(1),addrs(n_blocks),ierror)

        if (present_and_true(x)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! x
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%x(1),addrs(n_blocks),ierror)
        end if

        if (present_and_true(v)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! v
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%v(1),addrs(n_blocks),ierror)
        end if

#ifdef RWALK
        if (present_and_true(v_r)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! v_r
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%v_r(1),addrs(n_blocks),ierror)
        end if

        if (present_and_true(sdx)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1          ! sdx
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%sdx,addrs(n_blocks),ierror)
        end if
#endif

        if (present_and_true(f)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3          ! f
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f(1),addrs(n_blocks),ierror)
        end if

        if (present_and_true(q)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 4       ! q
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%q(0),addrs(n_blocks),ierror)
        end if

        if (present_and_true(w)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! w
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%w(1),addrs(n_blocks),ierror)
        end if

        if (present_and_true(t)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! t
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%t(1),addrs(n_blocks),ierror)
        end if
#ifdef FORCECOMPONENT
        if (present_and_true(f_normal)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_normal 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_normal(1),addrs(n_blocks),ierror)
        end if
        if (present_and_true(f_tangent)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_tangent 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_tangent(1),addrs(n_blocks),ierror)
        end if
        if (present_and_true(f_n)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_n 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_n(1),addrs(n_blocks),ierror)
        end if
        if (present_and_true(f_t)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_t 
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_t(1),addrs(n_blocks),ierror)
        end if
#endif	
        if (present_and_true(mag)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 2 ! mag
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%mag,addrs(n_blocks),ierror)
        end if

!#ifdef MD_MAG
        if (present_and_true(ftmi)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! fmi
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%fmi(1),addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! tmi
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%tmi(1),addrs(n_blocks),ierror)
         end if
!#endif

        if (present_and_true(ft_fluid)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_fluid
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_fluid(1),addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! t_fluid
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%t_fluid(1),addrs(n_blocks),ierror)
        end if

        if (present_and_true(ft_fluid_prev)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! f_fluid_prev
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%f_fluid_prev(1),addrs(n_blocks)&
                &,ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! t_fluid_prev
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%t_fluid_prev(1),addrs(n_blocks)&
                &,ierror)
        end if

        if (present_and_true(vws_fluid)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! v_fluid
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%v_fluid(1),addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! ws_fluid
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%ws_fluid(1),addrs(n_blocks),ierror)
        end if

        if (present_and_true(vws_fluid_avg)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! v_fluid_avg
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%v_fluid_avg(1),addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! ws_fluid_avg
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%ws_fluid_avg(1),addrs(n_blocks)&
                &,ierror)
        end if

        if (present_and_true(vws_fluid_acc)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! v_fluid_acc
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%v_fluid_acc(1),addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! ws_fluid_acc
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%ws_fluid_acc(1),addrs(n_blocks)&
                &,ierror)
        end if

        if (present_and_true(uid)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! uid
           types(n_blocks) = MPI_INTEGER
           call mpi_get_address(sample(1)%uid,addrs(n_blocks),ierror)
        end if

        if (present_and_true(vqw_new)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! vnew
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%vnew(1),addrs(n_blocks),ierror)

           if (use_rotation) then
              n_blocks = n_blocks+1
              lengths(n_blocks) = 4 ! qnew
              types(n_blocks) = MPI_REAL8
              call mpi_get_address(sample(1)%qnew(0),addrs(n_blocks),ierror)

              n_blocks = n_blocks+1
              lengths(n_blocks) = 3 ! wnew
              types(n_blocks) = MPI_REAL8
              call mpi_get_address(sample(1)%wnew(1),addrs(n_blocks),ierror)
           end if
        end if

        if (present_and_true(ws)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 3 ! ws
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%ws(1),addrs(n_blocks),ierror)
        end if

#ifdef PARTICLESTRESS
        if (present_and_true(tau)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 9 ! tau
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%tau(1,1),addrs(n_blocks),ierror)
        end if

#endif
        if (present_and_true(e_pot)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! e_pot
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%e_pot,addrs(n_blocks),ierror)
        end if

        if (present_and_true(master)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! master
           types(n_blocks) = MPI_INTEGER
           call mpi_get_address(sample(1)%master,addrs(n_blocks),ierror)
        end if

        if (present_and_true(R)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! R_orth
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%R_orth,addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! R_para
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%R_para,addrs(n_blocks),ierror)
        end if

#ifdef LADD_SURR_RHOF
        if (present_and_true(rhof)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! rhof
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%rhof,addrs(n_blocks),ierror)
        end if

        if (present_and_true(rhofn_acc)) then
           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! rhof_acc
           types(n_blocks) = MPI_REAL8
           call mpi_get_address(sample(1)%rhof_acc,addrs(n_blocks),ierror)

           n_blocks = n_blocks+1
           lengths(n_blocks) = 1 ! n_acc
           types(n_blocks) = MPI_INTEGER
           call mpi_get_address(sample(1)%n_acc,addrs(n_blocks),ierror)
        end if

#endif

         if (present_and_true(freeze)) then
             n_blocks = n_blocks+1
             lengths(n_blocks) = 3 ! freeze
             types(n_blocks) = MPI_REAL8
             call mpi_get_address(sample(1)%freeze(1),addrs(n_blocks),ierror)
          end if
        n_blocks = n_blocks+1
        lengths(n_blocks) = 1   ! next particle
        types(n_blocks) = MPI_UB
        call mpi_get_address(sample(2),addrs(n_blocks),ierror)

        call mpi_get_address(sample(1),base,ierror) ! base address
        displs(1:n_blocks) = addrs(1:n_blocks) - base

        call mpi_type_create_struct(n_blocks,lengths,displs,types,pmt,ierror)
        call checkmpi(ierror&
             &,'build_particle_mpitype(): mpi_type_create_struct() failed')

        call mpi_type_commit(pmt,ierror)
        call checkmpi(ierror&
             &,'build_particle_mpitype(): mpi_type_commit() failed')
    end subroutine build_particle_mpitype

    !> returns in \c n on process \c dst in \c comm_cart the total
    !> number of md particles on all processes.
    subroutine count_particles(n,dst)
        integer,intent(in) :: dst
        integer,intent(out) :: n
        integer ierror

        call mpi_reduce(nlocal,n,1,MPI_INTEGER,MPI_SUM,dst,comm_cart,ierror)
    end subroutine count_particles

    !> returns in \c n the total number of md particles on all
    !> processes and distributes it to all processes.
    subroutine count_particles_all(n)
        integer,intent(out) :: n
        integer ierror

        call mpi_allreduce(nlocal,n,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
    end subroutine count_particles_all

    subroutine sum_dx_all(systemdx)
        real(kind=rk),intent(out) :: systemdx
        integer ierror

        call mpi_allreduce(domaindx,systemdx,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
    end subroutine sum_dx_all

    !> return the interpolated fluid velocity at global postion  x  in  v  and
    !>  .true.  in  rock  if  x  is enclosed by rock sites.  x  must be at
    !> some place that is covered by the local  N , however, this subroutine
    !> is smart enough to deal with pbc. Velocities of the lbe particle
    !> species r/b/s are weighted as  mask(1/2/3) .
    !> Set  halo_extent>=2  if you use this subroutine for arbitrary positions
    !> inside each local domain!
 subroutine fluid_velocity(N,mask,x,v,rock)
   type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
   real(kind=rk),intent(in) :: mask(3), x(3)
   real(kind=rk),intent(out) :: v(3)
   logical,intent(out) :: rock
   integer l, s
   integer lx(3), lp(3), opp_lp(3) ! lattice point coordinates
   integer :: tz, tx, ty
   real(kind=rk) :: xx(3), v_lp(3), norm, weight
   real*8, dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_r
#ifndef SINGLEFLUID
   real*8, dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_b
#ifndef NOSURFACTANT
   real*8, dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_s
#endif
#endif
   call local_coordinates(x(:),xx(:))
   lx(:) = floor(xx(:))

   v(:) = 0.0_rk
   rock = .true.
   lattice_points: do l=1,n_lp
      lp(:) = lx(:) + lp_sur(:,l)
      ! where  rock_state/=0  there are no velocities
      lp_no_rock: if (N(lp(1),lp(2),lp(3))%rock_state==0.0_rk) then
         rock = .false.
         ! weight each lattice point with the volume of the hypercube
         ! spanned by the particle position and the opposite lattice
         ! point
         opp_lp(:) = lx(:) + opp_lp_sur(:,l)
         weight = abs(product(real(opp_lp(:),kind=rk)-xx(:)))
         v_lp(:) = 0.0_rk
#ifdef SINGLEFLUID
         call md_calculate_sc_forces(N,lp(1),lp(2),lp(3),f_r)
#else
#ifdef NOSURFACTANT
         call md_calculate_sc_forces(N,lp(1),lp(2),lp(3),f_b,f_r)
#else
         call md_calculate_sc_forces(N,lp(1),lp(2),lp(3),f_b,f_r,f_s)
#endif
#endif

              directions: do s=1,nnonrest ! ignore resting lbe particles
#ifdef BOUNCEBACK
                 ! ignore velocities that lead into rock sites, lbe
                 ! particles with these velocities will just bounce in
                 ! the next advection step, so in this timestep they
                 ! are at rest effectively (c(s,:)==0)
                 if (N(lp(1)+cx(s),lp(2)+cy(s),lp(3)+cz(s))%rock_state==0.0_rk) then
#endif
                    v_lp(:) = v_lp(:) + c(s,:)*g(s)*&
                         &( N(lp(1),lp(2),lp(3))%n_r(s)*amass_r*mask(1) &
#ifndef SINGLEFLUID
                         & + N(lp(1),lp(2),lp(3))%n_b(s)*amass_b*mask(2) &
#endif
#ifndef NOSURFACTANT
                         & + N(lp(1),lp(2),lp(3))%n_s(s)*amass_s*mask(3) &
#endif
                         &)
#ifdef BOUNCEBACK
                 end if
#endif
              end do directions

! force correction because Shan Chen forces
              v_lp(:) = v_lp(:)  + (&
                    & - tau_r*f_r(lp(1),lp(2),lp(3),:)*mask(1)/2.0d0&
#ifndef SINGLEFLUID
                    & - tau_b*f_b(lp(1),lp(2),lp(3),:)*mask(2)/2.0d0&
#endif
#ifndef NOSURFACTANT
                    & - tau_s*f_s(lp(1),lp(2),lp(3),:)*mask(3)/2.0d0&
#endif
                    &)


! calculation of mass : 
              norm = sum(g(:)*& ! This sum goes through rest vector, too!
                   &( N(lp(1),lp(2),lp(3))%n_r(:)*amass_r*mask(1)&
#ifndef SINGLEFLUID 
                   & + N(lp(1),lp(2),lp(3))%n_b(:)*amass_b*mask(2)&
#endif   
#ifndef NOSURFACTANT
                   & + N(lp(1),lp(2),lp(3))%n_s(:)*amass_s*mask(3)&
#endif
                   &) )

              v(:) = v(:) + weight*v_lp(:)/max(10.0E-9_rk,norm)
                          
              ! force extra-term
              tx = lp(1) + ccoords(1)*nx
              ty = lp(2) + ccoords(2)*ny
              tz = lp(3) + ccoords(3)*nz
              if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
                   (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
                   (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
!                 if (N(lp(1),lp(2),lp(3))%rock_state==0.0d0) then
                    ! external force term, in the zone where the
                    ! external force works
                    v(1) = v(1) - weight*(g_accn_x/2.0d0)
                    v(2) = v(2) - weight*(g_accn_y/2.0d0)
                    v(3) = v(3) - weight*(g_accn/2.0d0)
!                 end if
              endif

! vel corrections
          end if lp_no_rock
        end do lattice_points
    end subroutine fluid_velocity

    !> return the interpolated fluid velocity and dynamic viscosity at global
    !> postion  x  in  v  and  mu  and  .true.  in  rock  if  x  is enclosed by
    !> rock sites.  x  must be at some place that is covered by the local  N ,
    !> however, this subroutine is smart enough to deal with pbc. Velocities and
    !> dynamic viscosities of the lbe particle species r/b/s are weighted as
    !>  mask(1/2/3) .
    !> Set  halo_extent>=2  if you use this subroutine for arbitrary positions
    !> inside each local domain!
    subroutine fluid_velocity_and_viscosity(N,mask,x,v,mu,rock)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: mask(3),x(3)
        real(kind=rk),intent(out) :: v(3),mu
        logical,intent(out) :: rock
        integer l,s
        integer lx(3),lp(3),opp_lp(3) ! lattice point coordinates
        real(kind=rk) :: xx(3),v_lp(3),norm,weight
        real(kind=rk) :: mu_by_n(3)

        ! this saves some computations later
        mu_by_n(:) = mask(:)*(/amass_r,amass_b,amass_s/)*&
             &(2.0_rk*(/tau_r,tau_b,tau_s/)-1.0_rk)/6.0_rk

        call local_coordinates(x(:),xx(:))
        lx(:) = floor(xx(:))

        v(:) = 0.0_rk
        mu = 0.0_rk
        rock = .true.
        lattice_points: do l=1,n_lp
           lp(:) = lx(:) + lp_sur(:,l)
           ! where rock_state/=0 there are no velocities. Also dynamic
           ! viscosities there are assumed to be zero although this
           ! might be arguable...
           lp_no_rock: if (N(lp(1),lp(2),lp(3))%rock_state==0.0_rk) then
              rock = .false.
              ! weight each lattice point with the volume of the hypercube
              ! spanned by the particle position and the opposite lattice
              ! point
              opp_lp(:) = lx(:) + opp_lp_sur(:,l)
              weight = abs(product(real(opp_lp(:),kind=rk)-xx(:)))
              v_lp(:) = 0.0_rk
              directions: do s=1,nnonrest ! ignore resting lbe particles
#ifdef BOUNCEBACK
                 ! ignore velocities that lead into rock sites, lbe
                 ! particles with these velocities will just bounce in
                 ! the next advection step, so in this timestep they
                 ! are at rest effectively (c(s,:)==0)
                 if (N(lp(1)+cx(s),lp(2)+cy(s),lp(3)+cz(s))%rock_state==0.0_rk)&
                      & then
#endif
                    v_lp(:) = v_lp(:) + c(s,:)*g(s)*&
                         &(N(lp(1),lp(2),lp(3))%n_r(s)*amass_r*mask(1)&
#ifndef SINGLEFLUID
                         &+N(lp(1),lp(2),lp(3))%n_b(s)*amass_b*mask(2)&
#endif
#ifndef NOSURFACTANT
                         &+N(lp(1),lp(2),lp(3))%n_s(s)*amass_s*mask(3)&
#endif
                         &)
#ifdef BOUNCEBACK
                 end if
#endif
              end do directions
              norm = sum(g(:)*& ! sum loops through rest vector, too!
                   &(N(lp(1),lp(2),lp(3))%n_r(:)*amass_r*mask(1)&
#ifndef SINGLEFLUID
                   &+N(lp(1),lp(2),lp(3))%n_b(:)*amass_b*mask(2)&
#endif
#ifndef NOSURFACTANT
                   &+N(lp(1),lp(2),lp(3))%n_s(:)*amass_s*mask(3)&
#endif
                   &))
              v(:) = v(:) + weight*v_lp(:)/max(10.0E-9_rk,norm)

              mu = mu + sum(g(:)*&
                   &(N(lp(1),lp(2),lp(3))%n_r(:)*mu_by_n(1)&
#ifndef SINGLEFLUID
                   &+N(lp(1),lp(2),lp(3))%n_b(:)*mu_by_n(2)&
#endif
#ifndef NOSURFACTANT
                   &+N(lp(1),lp(2),lp(3))%n_s(:)*mu_by_n(3)&
#endif
                   &))*weight
           end if lp_no_rock
        end do lattice_points
    end subroutine fluid_velocity_and_viscosity

    !> returns true if all lattice points surrounding \c pos are rock.
    !>
    !> \c pos must be on the local process!
    logical function inside_rock(N,pos)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: pos(3)
        integer fp(3)

        fp(:) = 1 + floor(pos(:)) - start(:)
        inside_rock = &
             &all(N(fp(1)+lp_sur(1,:),fp(2)+lp_sur(2,:),fp(3)+lp_sur(3,:))&
             &%rock_state/=0.0_rk)
    end function inside_rock

    !> returns the maximum radius in the simulation
    !>
    !> \returns upper boundary for \c R_orth and \c R_para over all particles
    !>
    !> \warning The value is obtained from the input file parameters and
    !> not---in case of polydispersity---from the particles
    !> themselves. This could be a tricky source of problems when
    !> restoring from a checkpoint.
    function maximum_particle_radius()
        real(kind=rk) :: maximum_particle_radius

        select case (initial_radii)
        case ('none')
           maximum_particle_radius = max(R_orth,R_para)
        case ('normal', 'random', 'volume')
           maximum_particle_radius = R_max
        case default
           call error_md('unknown value initial_radii=<'//trim(initial_radii)&
                &//'>---choose either <none>, <normal>, <random>, or <volume>')
        end select
    end function maximum_particle_radius

!> same as  \c lbe_make_filename_output()  in \p lbe_io.F90 except that it uses
!> particle id instead of time.
subroutine md_make_filename_particles(buffer, prefix, suffix, pid)
  implicit none

  character(len=*)   :: buffer, prefix, suffix
  integer            :: pid
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_p) pid, chk_uid

  write(buffer,"('./',A,'/',A,'_',A,'_',A,A)") &
    trim(folder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)

end subroutine md_make_filename_particles

    !> Returns the number density of md particles in \c nd.
    !>
    !> Of course, the density calculated this way is meaningless if the
    !> particles are distributed inhomogeneously in the simulation
    !> space.
    subroutine number_density(nd)
        real(kind=rk),intent(out) :: nd
        integer :: n

        call count_particles_all(n)
        nd = real(n/product(tsize),kind=rk)
    end subroutine number_density

    !> Return in \c o a unit vector that points into the direction of the
    !> symmetry axis of a particle with quaternion \c q in space fixed
    !> coordinates.
    pure function orientation(q)
        real(kind=rk),intent(in) :: q(0:3)
        real(kind=rk) :: orientation(3)

        ! this is A^{-1}\cdot(0 0 1)^t
        ! The sign of second component is changed by me[Xie].
        orientation = (/&
             &2.0_rk*(q(1)*q(3)+q(0)*q(2)),&
             &2.0_rk*(q(2)*q(3)-q(0)*q(1)),&
             &q(0)**2-q(1)**2-q(2)**2+q(3)**2&
             &/)
    end function orientation

    !> returns the half axes of a particle independently on whether \c
    !> polydispersity is set or not
    !>
    !> \param[in] p particle
    !>
    !> \returns particle \c (/R_orth,R_orth,R_para/) as effective for
    !> \c p
    pure function particle_half_axes(p)
        real(kind=rk) :: particle_half_axes(3)
        type(md_particle_type),intent(in) :: p
        real(kind=rk) :: io,ip

        if (polydispersity) then
           particle_half_axes = (/p%R_orth,p%R_orth,p%R_para/)
        else
           particle_half_axes = (/R_orth,R_orth,R_para/)
        end if
    end function particle_half_axes

    !> returns the diagonal of the moment of inertia tensor for a
    !> given particle in its body-fixed coordinate system
    !> independently on whether \c polydispersity is set or not
    !>
    !> \param[in] p particle
    !>
    !> \returns \c (/inertia_orth,inertia_orth,inertia_para/) as
    !> effective for \c p
    pure function particle_inertia(p)
        real(kind=rk),dimension(3) :: particle_inertia
        type(md_particle_type),intent(in) :: p
        real(kind=rk) :: orth,para

        if (polydispersity) then
           orth = p%R_orth**2 + p%R_para**2
           para = 2.0_rk * p%R_orth**2
           particle_inertia = (/orth,orth,para/)*particle_mass(p)/5.0_rk
        else
           particle_inertia = (/inertia_orth,inertia_orth,inertia_para/)
        end if
    end function particle_inertia

    !> returns the inverse of the half axes of a particle
    !> independently on whether \c polydispersity is set or not
    !>
    !> \param[in] p particle
    !>
    !> \returns particle \c 1.0_rk/(/R_orth,R_orth,R_para/) as
    !> effective for \c p
    pure function particle_inverse_half_axes(p)
        real(kind=rk) :: particle_inverse_half_axes(3)
        type(md_particle_type),intent(in) :: p
        real(kind=rk) :: io,ip

        if (polydispersity) then
           io = 1.0_rk/p%R_orth
           ip = 1.0_rk/p%R_para
        else
           io = 1.0_rk/R_orth
           ip = 1.0_rk/R_para
        end if
        particle_inverse_half_axes = (/io,io,ip/)
    end function particle_inverse_half_axes

    !> returns the mass for a given particle independently on whether
    !> \c polydispersity is set or not
    !>
    !> \param[in] p particle
    !>
    !> \returns particle mass as effective for \c p
    pure function particle_mass(p)
        real(kind=rk) :: particle_mass
        type(md_particle_type),intent(in) :: p

        if (polydispersity) then
           particle_mass = rho*4.0_rk*pi*p%R_para*p%R_orth**2/3.0_rk
        else
           particle_mass = mass
        end if
    end function particle_mass

    !> returns the two half axis effective for a given particle
    !> independently on whether \c polydispersity is set or not
    !>
    !> \param[in] p particle
    !>
    !> \returns \c (/R_orth,R_para/) as effective for \c p
    pure function particle_radii(p)
        real(kind=rk),dimension(2) :: particle_radii
        type(md_particle_type),intent(in) :: p

        if (polydispersity) then
           particle_radii = (/p%R_orth,p%R_para/)
        else
           particle_radii = (/R_orth,R_para/)
        end if
    end function particle_radii

    !> returns the quaternion representation of the orientation
    !> specified by the vector \c o
    !>
    !> There might be bugs in the calculation process. Don't use it for
    !> important things without further checks...
    pure function quaternion(o)
      real(kind=rk) :: quaternion(0:3)
      real(kind=rk),intent(in) :: o(3)
      real(kind=rk) :: o3
      real(kind=rk) phi,theta,p2,t2

      ! Don't ask me why, but this seems to make quaternion <-> orientation reversible - SF
      !if (o(2) > 0.0_rk) then
      !  o3 = -o(3)
      !else
      !  o3 =  o(3)
      !endif
      ! End of magic

      !if (o3==0.0_rk) then
      !  theta = 0.5_rk*pi
      !else
      !  theta = atan(sqrt(o(1)**2+o(2)**2)/o3)
      !End of magic

      !if (o(3)==0.0_rk) then
      !  theta = 0.5_rk*pi
   !   else 
   !   theta = atan(sqrt(o(1)**2+o(2)**2)/o3) + pi
      !end if

      !if (o(2)==0.0_rk) then
      !  if (o(1)>0.0_rk) then
      !    phi = 0.5_rk*pi
      !  else
      !    phi = 1.5_rk*pi
      !  end if
      !else
      !  phi = atan(-o(1)/o(2))
      !end if
     !!!New code. The reason is that the range of function atan is [-pi/2,pi/2]. But atan2 
     ! could cover [-pi,pi]
     theta = atan2(sqrt(o(1)**2+o(2)**2),o(3))
     phi = atan2(o(1),-o(2))
      ! build quaternion from Euler-angles (3rd angle is zero because it is
      ! irrelevant because of the rotational symmetry of our particles)
      p2 = 0.5_rk*phi
      t2 = 0.5_rk*theta
      quaternion = (/cos(p2)*cos(t2),cos(p2)*sin(t2),sin(p2)*sin(t2)&
           &,sin(p2)*cos(t2)/)
    end function quaternion

    !> This is just a debug function to check if the quaternion <-> orientation transformation is reversible.
    subroutine check_quaternion()
      real(kind=rk), dimension(8,3) :: testr
      real(kind=rk), dimension(3) :: testo
      real(kind=rk), dimension(4) :: testq
      real(kind=rk) :: s1,s2,s3
      character(len=1), dimension(4) :: display
      integer k,l

      ! Make sure this will be length 1
      s1 = 1.0d0 / sqrt(2.0d0)
      s2 = 1.0d0 / sqrt(4.0d0)
      s3 = 1.0d0 / sqrt(4.0d0)

      testr(1,1) = -s1 ; testr(1,2) = -s2 ; testr(1,3) = -s3
      testr(2,1) = -s1 ; testr(2,2) = -s2 ; testr(2,3) =  s3
      testr(3,1) = -s1 ; testr(3,2) =  s2 ; testr(3,3) = -s3
      testr(4,1) = -s1 ; testr(4,2) =  s2 ; testr(4,3) =  s3
      testr(5,1) =  s1 ; testr(5,2) = -s2 ; testr(5,3) = -s3
      testr(6,1) =  s1 ; testr(6,2) = -s2 ; testr(6,3) =  s3
      testr(7,1) =  s1 ; testr(7,2) =  s2 ; testr(7,3) = -s3
      testr(8,1) =  s1 ; testr(8,2) =  s2 ; testr(8,3) =  s3

      do k =1,8

        testo(1) = testr(k,1)
        testo(2) = testr(k,2)
        testo(3) = testr(k,3)

        do l = 1,3
          if (testo(l) < 0.0d0) then
            display(l) = "-"
          else
            display(l) = "+"
          endif
        enddo

        write(msgstr,"('Testing ',I0,X,3(A,X),X,3(F8.4,X))") k, display(1), display(2), display(3), testo
        call log_msg_md(msgstr)

        testq = quaternion(testo)

        do l = 1,4
          if (testq(l) < 0.0d0) then
            display(l) = "-"
          else
            display(l) = "+"
          endif
        enddo

        write(msgstr,"('        ',I0,X,4(A,X),X,4(F8.4,X))") k, display(1), display(2), display(3), display(4), testq
        call log_msg_md(msgstr)

        testo = orientation(testq)

        do l = 1,3
          if (testo(l) < 0.0d0) then
            display(l) = "-"
          else
            display(l) = "+"
          endif
        enddo

        write(msgstr,"('        ',I0,X,3(A,X),X,3(F8.4,X))") k, display(1), display(2), display(3), testo
        call log_msg_md(msgstr)
      enddo
    end subroutine check_quaternion

    !> root scatters the global rock_state array rs(:,:,:) to the
    !> local lattice chunks N(:,:,:)%rock_state of all processes. Also
    !> the halo with extent halo_extent lu for every process is
    !> scattered.  rs must have the dimensions
    !> 1-halo_extent:(/tnx,tny,tnz/)+halo_extent .
    subroutine scatter_rock_state_halo(rs,N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: &
             &rs(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),allocatable,dimension(:) :: cbuf ! communication buffer
        integer pcoords(3,0:nprocs-1) ! cartesian topology coordinates
        integer os(3)           ! lattice offset for different processors
        integer i,j,k,p,r,ierror,stat,he,h1,csize,status(MPI_STATUS_SIZE)
        integer,parameter :: tag=0

        he = halo_extent
        h1 = halo_extent - 1

        ! amount of data sent/recved
        csize = product(chunksize(:)+2*halo_extent)

        allocate (cbuf(csize),stat=stat)
        call check_allocate(stat,'scatter_rock_state_halo(): cbuf')

        rank0: if (myrankc==0) then
           call find_topology(pcoords)
           ! first copy rock_state to local chunk...
           os(:) = pcoords(:,0)*(/nx,ny,nz/)
           N(1-he:nx+he,1-he:ny+he,1-he:nz+he)%rock_state = rs(&
                &os(1)-h1:os(1)+nx+he,os(2)-h1:os(2)+ny+he,os(3)-h1:os(3)+nz+he&
                &)

           ! ...then send rock_state to the other processes
           rank: do r=1,nprocs-1
              os(:) = pcoords(:,r)*(/nx,ny,nz/)
              p = 0
              do i=os(1)-h1,os(1)+nx+he
                 do j=os(2)-h1,os(2)+ny+he
                    do k=os(3)-h1,os(3)+nz+he
                       p = p + 1
                       cbuf(p) = rs(i,j,k)
                    end do
                 end do
              end do
              call mpi_send(cbuf,csize,MPI_REAL8,r,tag,comm_cart,ierror)
           end do rank
        else rank0
           ! receive rock_state from root process
           call mpi_recv(cbuf,csize,MPI_REAL8,0,tag,comm_cart,status,ierror)

           ! copy it to local chunk
           p = 0
           do i=-h1,nx+he
              do j=-h1,ny+he
                 do k=-h1,nz+he
                    p = p + 1
                    N(i,j,k)%rock_state = cbuf(p)
                 end do
              end do
           end do
        end if rank0

        deallocate(cbuf)
    end subroutine scatter_rock_state_halo

    !> returns a rotation matrix that transforms from the space-fixed
    !> coordinate frame to the body fixed one that is defined by the
    !> quaternion q
    !>
    !> see Allen/Tildesley, eq. (3.36)
    function space_to_body_matrix(q)
        real(kind=rk) :: space_to_body_matrix(3,3)
        real(kind=rk),intent(in) :: q(0:3)

        space_to_body_matrix = reshape((/&
             &q(0)**2+q(1)**2-q(2)**2-q(3)**2,&
             &2.0_rk*(q(1)*q(2)-q(0)*q(3)),&
             &2.0_rk*(q(1)*q(3)+q(0)*q(2)),&

             &2.0_rk*(q(1)*q(2)+q(0)*q(3)),&
             &q(0)**2-q(1)**2+q(2)**2-q(3)**2,&
             &2.0_rk*(q(2)*q(3)-q(0)*q(1)),&

             &2.0_rk*(q(1)*q(3)-q(0)*q(2)),&
             &2.0_rk*(q(2)*q(3)+q(0)*q(1)),&
             &q(0)**2-q(1)**2-q(2)**2+q(3)**2&
             &/),&
             &(/3,3/))
    end function space_to_body_matrix

    !> manually invalidate the neighbor lists, so reneighboring
    !> happens and \c borders() is called instead of \c communicate()
    subroutine trigger_list_update
        list_update_required = .true.
    end subroutine trigger_list_update

    !> prints a warning message to standard output
    subroutine warning(msg)
        character*(*), intent(in) :: msg
        print *,'WARNING: ',msg
    end subroutine warning

#endif
end module lbe_md_helper_module

