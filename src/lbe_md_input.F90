#include "lbe.h"

!> input for md part of the code
module lbe_md_input_module
#ifdef MD
    use lbe_globals_module, only: nd,input_dfile_unit,maxpos,minpos,tsize&
         &,n_spec,myrankc
    use lbe_log_module
    use lbe_md_bc_leesedwards_module, only: md_leesedwards_reproducible_cp&
         &,no_md_leesedwards
    use lbe_md_dynamic_module, only: first_inserted_puid,new_uids
    use lbe_md_fluid_ladd_mc_module, only: pfr, pfb, pfg, global_mass_change, global_mass_target
    use lbe_md_globals_module
    use lbe_md_growing_stage_module, only: growing_stage
    use lbe_md_helper_module, only: log_msg_md,error_md,log_msg_md_hdr
    use lbe_md_init_module, only: init_file,initial_orientations&
         &,initial_rotations,initial_velocities,min_dist,phi,q0,temp0,v0,w0&
         &,x0hi,x0lo
    use lbe_md_lyapunov_module, only: lyapunov
    use lbe_md_module, only: average_vw_fluid,boundary_width,constant_f,constant_fc,cargo_par,cargo_par_id &
         &,constant_t,constant_tc,constant_f_y_min,gaussian_plane_active,gaussian_plane_pos,&
         &gaussian_plane_E,gaussian_plane_c,gaussian_plane_sqwidth,randompartforce,fix_v,fix_w, fix_vx, fix_vx_t,fix_vz, fix_vz_t&
	      &,runtime_particle_mod,rpm_miny,rpm_numpart, rpm_addintervall,rpm_deleteintervall,rpm_minpartdist&
         &,setup_list&
         &,ramping,ramptime,vmax_ramping,forcebalance
    use lbe_md_parallel_module, only: scatter_particles
    use lbe_md_output_module, only: dump_double,dump_forces,dump_format&
         &,dump_ids,dump_quaternions,dump_orientations,dump_pnd&
         &,dump_potentials,dump_positions,dump_radii,dump_rotations&
         &,dump_summary&
#ifdef PARTICLESTRESS
         &,dump_stresses,dump_stresses_b&
#endif
#ifdef LADD_SURR_RHOF
         &,dump_surr_rhof&
#endif
!#ifdef  MD_MAG
        &,dump_mag&
!#endif
         &,dump_time,dump_torques,dump_velocities,msdx_len,msdx_minx,msdx_miny&
         &,msdx_minz,msdx_maxx,msdx_maxy,msdx_maxz,msdx_name,msdx_start&
         &,msdx_subtract_drift,msdy_len,msdy_minx,msdy_miny,msdy_minz,msdy_maxx&
         &,msdy_maxy,msdy_maxz,msdy_name,msdy_start,msdy_subtract_drift&
         &,msdz_len,msdz_minx,msdz_miny,msdz_minz,msdz_maxx,msdz_maxy,msdz_maxz&
         &,msdz_name,msdz_start,msdz_subtract_drift,n_dump,n_dump_substep&
         &,n_dump_substep_start,n_rvz_profile_dump,n_rvz_profile_sample,n_msdx&
         &,n_msdx_max,n_msdx_sample,n_msdy,n_msdy_max,n_msdy_sample,n_msdz&
         &,n_msdz_max,n_msdz_sample
    use lbe_md_potential_module, only: potential
    use lbe_md_magnetic_module, only: cons_mag_dip,cons_mag_fie, gra_mag_fie&
        &,T_switch, T_fre, T_con, md_magnetic
    use lbe_md_rock_module, only: rock, rock_friction
    use lbe_parallel_module, only: check_allocate,comm_cart,nprocs
    use lbe_parms_module, only: boundary,chk_uid,inp_file,nt,restore_string&
         &,sci_profile,arg_input_dfile_set,arg_input_dfile
    use lbe_io_helper_module, only: lbe_make_filename_cp&
         &,lbe_make_filename_restore
    implicit none
    include 'mpif.h'
    private

    public input,md_restore_checkpoint

    namelist /md_input/ alat,average_ft_fluid,average_vw_fluid,boundary_width&
         &,constant_f, constant_fc,cargo_par,cargo_par_id,constant_t,constant_tc,constant_f_y_min&
	 &,gaussian_plane_active,gaussian_plane_pos,gaussian_plane_E,gaussian_plane_c,gaussian_plane_sqwidth,randompartforce&
	&,runtime_particle_mod,rpm_miny,rpm_numpart,rpm_addintervall,rpm_deleteintervall,rpm_minpartdist&
         &,count_periodic,cutoff_v,delta_m,delta_t,delta_x,dfz_maxx,dfz_minx&
         &,diffuse_x,diffuse_y,diffuse_z,dump_double,dump_forces,dump_format&
         &,dump_ids,dump_orientations,dump_pnd,dump_potentials,dump_positions&
         &,dump_quaternions,dump_radii,dump_rotations&
         &, dump_summary&
#ifdef PARTICLESTRESS
         &,dump_stresses,dump_stresses_b&
#endif
#ifdef LADD_SURR_RHOF
         &,dump_surr_rhof&
#endif
!#ifdef MD_MAG
         &,dump_mag&
!#endif
         &,dump_time,dump_torques,dump_velocities,fix_v,fix_w,fix_vx,fix_vx_t,fix_vz, fix_vz_t,growing_stage, md_magnetic&
         &,ineigh,inertia_orth,inertia_para,init_file,initial_orientations&
         &,initial_placing,initial_rotations,initial_velocities,interaction&
         &,lyapunov,mass,md_leesedwards_reproducible_cp,mean_free_path,min_dist&
         &,molec_mass,molec_mass_lbm,msdx_name,msdx_start,msdx_subtract_drift&
         &,n_msdx,n_msdx_sample,msdx_len,msdx_minx,msdx_miny,msdx_minz&
         &,msdx_maxx,msdx_maxy,msdx_maxz,msdy_name,msdy_start&
         &,msdy_subtract_drift,n_msdy,n_msdy_sample,msdy_len,msdy_minx&
         &,msdy_miny,msdy_minz,msdy_maxx,msdy_maxy,msdy_maxz,msdz_name&
         &,msdz_start,msdz_subtract_drift,n_msdz,n_msdz_sample,msdz_len&
         &,msdz_minx,msdz_miny,msdz_minz,msdz_maxx,msdz_maxy,msdz_maxz&
         &,n_rvz_profile_dump,n_dfz,n_dump,n_dump_substep,n_dump_substep_start&
         &,n_rvz_profile_sample,n_stat,new_uids,no_md_leesedwards,np_sphere&
         &,np_sphere_cap,phi,potential,q0,rc,reenter_rdm_max,reenter_rdm_min&
         &,time_output,prtcl_output,reflective_rocks,rho,rock,rock_friction,rs,semiper_max&
         &,semiper_min,temp0,steps_per_lbe_step,temperat,v0,w0,x0hi,x0lo&
         &,ramping,vmax_ramping,ramptime,forcebalance

contains

    !> read input file and broadcast its content
    subroutine input
        integer i,ierror

        call log_msg_md_hdr("Reading MD input")

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_input,err=110)
           close (unit=md_input_file_unit,err=120)
           ! call log_msg_md('Read /md_input/ from file "'//trim(inp_file)//'.md"')
           !write (6,nml=md_input)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_input, IOSTAT = ierror)
          if (ierror .ne. 0) then
              call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('n_dump             = ',I0,', dump_format = <',A,'>')") n_dump, trim(dump_format)
        call log_msg(msgstr)
        write(msgstr,"('dump_double        = ',L1)") dump_double
        call log_msg(msgstr)
        write(msgstr,"('dump_ids           = ',L1,', dump_positions    = ',L1)") dump_ids, dump_positions
        call log_msg(msgstr)
        write(msgstr,"('dump_velocities    = ',L1,', dump_forces       = ',L1)") dump_velocities, dump_forces
        call log_msg(msgstr)
        write(msgstr,"('dump_quaternions   = ',L1,', dump_orientations = ',L1)") dump_quaternions, dump_orientations
        call log_msg(msgstr)
        write(msgstr,"('dump_rotations     = ',L1,', dump_torques      = ',L1)") dump_rotations, dump_torques
        call log_msg(msgstr)
!#ifdef MD_MAG
        write(msgstr,"('dump_mag     = ',L1)") dump_mag
        call log_msg(msgstr)  
!#endif
        write (msgstr&
             &,"('dump_potentials    = ',L1,', dump_pnd          = ',L1)") &
             &dump_potentials,dump_pnd
        call log_msg(msgstr)
        write(msgstr,"('dump_radii         = ',L1)") dump_radii
        call log_msg(msgstr)
#ifdef PARTICLESTRESS
        write(msgstr,"('dump_stresses      = ',L1,', dump_stresses_b   = ',L1)"&
             &) dump_stresses,dump_stresses_b
        call log_msg(msgstr)
#endif
        write(msgstr,"('dump_summary       = ',L1)") dump_summary
        call log_msg(msgstr)
#ifdef LADD_SURR_RHOF
        write(msgstr,"('dump_surr_rhof     = ',L1)") dump_surr_rhof
        call log_msg(msgstr)
#endif

        write(msgstr,"('dump_time          = ',L1)") dump_time
        call log_msg(msgstr)
        write(msgstr,"('n_dump_substep     = ',I0)") n_dump_substep
        call log_msg(msgstr)
        if (n_dump_substep/=0) then
           write(msgstr,"('n_dump_substep_start=',I0)") n_dump_substep_start
           call log_msg(msgstr)
        end if
        write(msgstr,"('n_stat             = ',I0)") n_stat
        call log_msg(msgstr)
        write(msgstr,"('count_periodic     = ',L1)") count_periodic
        call log_msg(msgstr)
        call log_ws()
        write(msgstr,"('n_rvz_profile_dump = ',I0,', n_rvz_profile_sample = ',I0)") n_rvz_profile_dump, n_rvz_profile_sample
        call log_msg(msgstr)
        write(msgstr,"('n_dfz              = ',I0)") n_dfz
        call log_msg(msgstr)
        write(msgstr,"('dfz_minx           = ',I0,', dfz_maxx = ',I0)") dfz_minx, dfz_maxx
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('interaction        = <',A,'>')") trim(interaction)
        call log_msg(msgstr)
        write(msgstr,"('potential          = <',A,'>')") trim(potential)
        call log_msg(msgstr)
        write(msgstr,"('rock               = <',A,'>')") trim(rock)
        call log_msg(msgstr)
        write(msgstr,"('rock_friction               = <',A,'>')") trim(rock_friction)
        call log_msg(msgstr)
        write(msgstr,"('growing_stage      = ',L1)") growing_stage
        call log_msg(msgstr)
        write(msgstr,"('md_magnetic      = ',L1)") md_magnetic
        call log_msg(msgstr)
        write(msgstr,"('lyapunov           = ',L1)") lyapunov
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('no_md_leesedwards              = ',L1)") &
             &no_md_leesedwards
        CALL log_msg(trim(msgstr),.false.)
        write(msgstr,"('md_leesedwards_reproducible_cp = ',L1)") &
             &md_leesedwards_reproducible_cp
        CALL log_msg(trim(msgstr),.false.)
        CALL log_ws(.false.)

        write(msgstr,"('initial_placing    = <',A,'>')") trim(initial_placing)
        call log_msg(msgstr)
        if (trim(initial_placing) .eq. 'file:asc/xvow' .or. &
            trim(initial_placing) .eq. 'file:asc/xvqw' .or. &
            trim(initial_placing) .eq. 'file_parallel:asc/xvow' .or. &
            trim(initial_placing) .eq. 'file_parallel:asc/xvowrm' .or. &
            trim(initial_placing) .eq. 'file_parallel:asc/xvowm' .or. &
            trim(initial_placing) .eq. 'file_parallel:asc/xvowr' ) then
          write(msgstr,"('  init_file        = <',A,'>')") trim(init_file)
          call log_msg(msgstr)
        endif
        if (trim(initial_placing) .eq. 'facex' .or. &
            trim(initial_placing) .eq. 'facey' .or. &
            trim(initial_placing) .eq. 'facez' .or. &
            trim(initial_placing) .eq. 'fcc'   .or. &
            trim(initial_placing) .eq. 'sc' ) then
          write(msgstr,"('  alat             = ',F16.10)") alat
          call log_msg(msgstr)
          write(msgstr,"('  x0lo             = (',2(F16.10,','),F16.10,')')") &
               &x0lo
          call log_msg(msgstr)
          write(msgstr,"('  x0hi             = (',2(F16.10,','),F16.10,')')") &
               &x0hi
          call log_msg(msgstr)
        endif
        if (trim(initial_placing) .eq. 'sphere') then
          write(msgstr,"('  np_sphere        = ',I0)") np_sphere
          call log_msg(msgstr)
          write(msgstr,"('  np_sphere_cap    = ',I0)") np_sphere_cap
          call log_msg(msgstr)
        endif
        if (trim(initial_placing)=='random') then
          write (msgstr,"('  phi              = ',F16.10)") phi
          call log_msg(msgstr)
          if (phi<0.0_rk) then
             write (msgstr,"('  alat             = ',F16.10)") alat
             call log_msg(msgstr)
          end if
          write (msgstr,"('  min_dist         = ',F16.10)") min_dist
          call log_msg(msgstr)
          write(msgstr,"('  x0lo             = (',3(F16.10,:,','),')')") x0lo
          call log_msg(msgstr)
          write(msgstr,"('  x0hi             = (',3(F16.10,:,','),')')") x0hi
          call log_msg(msgstr)
        endif
        call log_ws()

        write(msgstr,"('initial_velocities   = <',A,'>')") trim(initial_velocities)
        call log_msg(msgstr)
        if (trim(initial_velocities) .eq. 'v0') then
          write(msgstr,"('  v0                 = (',F16.10,',',F16.10,',',F16.10,')')") v0(1), v0(2), v0(3)
          call log_msg(msgstr)
        endif
        if (trim(initial_velocities) .eq. 'temp0') then
          write(msgstr,"('  temp0              = ',F16.10)") temp0
          call log_msg(msgstr)
        endif
        write(msgstr,"('initial_orientations = <',A,'>')") trim(initial_orientations)
        call log_msg(msgstr)
        if (trim(initial_orientations) .eq. 'q0') then
          write(msgstr,"('  q0                 = (',F16.10,',',F16.10,',',F16.10,',',F16.10,')')") q0(0), q0(1), q0(2), q0(3)
          call log_msg(msgstr)
        endif
        write(msgstr,"('initial_rotations    = <',A,'>')") trim(initial_rotations)
        call log_msg(msgstr)
        if (trim(initial_rotations) .eq. 'w0') then
          write(msgstr,"('  w0                 = (',F16.10,',',F16.10,',',F16.10,')')") w0(1), w0(2), w0(3)
          call log_msg(msgstr)
        endif
        call log_ws()

        if (boundary=='periodic_inflow') then
           write(msgstr,"('new_uids           = <',A,'>')") trim(new_uids)
           call log_msg(msgstr)
           call log_ws()
        end if

        write(msgstr,"('delta_x            = ',F16.10)") delta_x
        call log_msg(msgstr)
        write(msgstr,"('delta_t            = ',F16.10)") delta_t
        call log_msg(msgstr)
        write(msgstr,"('delta_m            = ',F16.10)") delta_m
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('steps_per_lbe_step = ',I0)") steps_per_lbe_step
        call log_msg(msgstr)
        write(msgstr,"('average_ft_fluid   = ',L1)") average_ft_fluid
        call log_msg(msgstr)
        write(msgstr,"('average_vw_fluid   = ',L1)") average_vw_fluid
        call log_msg(msgstr)

        write(msgstr,"('rho                = ',F16.10)") rho
        call log_msg(msgstr)
        write(msgstr,"('mass               = ',F16.10)") mass
        call log_msg(msgstr)
        write(msgstr,"('inertia_orth       = ',F16.10)") inertia_orth
        call log_msg(msgstr)
        write(msgstr,"('inertia_para       = ',F16.10)") inertia_para
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('ineigh             = ',I0)") ineigh
        call log_msg(msgstr)
        write(msgstr,"('rc                 = ',F16.10)") rc
        call log_msg(msgstr)
        write(msgstr,"('rs                 = ',F16.10)") rs
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('semiper_max        = (',L1,',',L1,',',L1,')')") semiper_max(1), semiper_max(2), semiper_max(3)
        call log_msg(msgstr)
        write(msgstr,"('semiper_min        = (',L1,',',L1,',',L1,')')") semiper_min(1), semiper_min(2), semiper_min(3)
        call log_msg(msgstr)
        write(msgstr,"('diffuse_[xyz]      = (',L1,',',L1,',',L1,')')") diffuse_x, diffuse_y, diffuse_z
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('reenter_rdm_min    = (',L1,',',L1,',',L1,')')") reenter_rdm_min(1), reenter_rdm_min(2), reenter_rdm_min(3)
        call log_msg(msgstr)
        write(msgstr,"('reenter_rdm_max    = (',L1,',',L1,',',L1,')')") reenter_rdm_max(1), reenter_rdm_max(2), reenter_rdm_max(3)
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('reflective_rocks   = ',L1)") reflective_rocks
        call log_msg(msgstr)
        write(msgstr,"('molec_mass         = ',F16.10)") molec_mass
        call log_msg(msgstr)
        write(msgstr,"('molec_mass_lbm     = ',F16.10)") molec_mass_lbm
        call log_msg(msgstr)
        write(msgstr,"('temperat           = ',F16.10)") temperat
        call log_msg(msgstr)
        write(msgstr,"('mean_free_path     = ',F16.10)") mean_free_path
        call log_msg(msgstr)
        call log_ws()

        write(msgstr,"('cutoff_v           = ',F16.10)") cutoff_v
        call log_msg(msgstr)

        write(msgstr,"('fix_v              = ',L1,', fix_w = ',L1, ', fix_vx = ',L1)") fix_v, fix_w, fix_vx
        call log_msg(msgstr)
        write(msgstr,"('fix_vx_t        = ',I0)") fix_vx_t
        call log_msg(msgstr)
        write(msgstr,"('fix_vz   = ',L1)") fix_vz
        call log_msg(msgstr)
        write(msgstr,"('fix_vz_t        = ',I0)") fix_vz_t
        call log_msg(msgstr)

        write(msgstr,"('constant_f         = (',F16.10,',',F16.10,',',F16.10,')')") constant_f(1), constant_f(2), constant_f(3)
        call log_msg(msgstr)
        !> for cargo particle
        write(msgstr,"('cargo_par      = ',L1)") cargo_par
          call log_msg(msgstr)
        write(msgstr,"('cargo_par_id        = ',I0)") cargo_par_id
        call log_msg(msgstr)
        write(msgstr,"('constant_fc         =(',F16.10,',',F16.10,',',F16.10,')')") constant_fc(1), constant_fc(2),constant_fc(3)
        call log_msg(msgstr)
        write(msgstr,"('constant_t         = (',F16.10,',',F16.10,',',F16.10,')')") constant_t(1), constant_t(2), constant_t(3)
        call log_msg(msgstr)
        write(msgstr,"('constant_tc         = (',F16.10,',',F16.10,',',F16.10,')')") constant_tc(1), constant_tc(2), constant_tc(3)
        call log_msg(msgstr)

       
        write(msgstr,"('constant_f_y_min         = ',F16.10)") constant_f_y_min
        call log_msg(msgstr)
        write(msgstr,"('gaussian_plane_active         = ',L1)") gaussian_plane_active
        call log_msg(msgstr)
        write(msgstr,"('gaussian_plane_pos         = ',F16.10)") gaussian_plane_pos
        call log_msg(msgstr)
        write(msgstr,"('gaussian_plane_E         = ',F16.10)") gaussian_plane_E
        call log_msg(msgstr)
        write(msgstr,"('gaussian_plane_c         = ',F16.10)") gaussian_plane_c
        call log_msg(msgstr)
        write(msgstr,"('gaussian_plane_sqwidth         = ',F16.10)") gaussian_plane_sqwidth
        call log_msg(msgstr)
        write(msgstr,"('boundary_width     = ',F16.10)") boundary_width
        call log_msg(msgstr)
        write(msgstr,"('randompartforce         = ',L1)") randompartforce
        call log_msg(msgstr)

        write(msgstr,"('runtime_particle_mod         = ',L1)") runtime_particle_mod
        call log_msg(msgstr)
        write(msgstr,"('rpm_miny         = ',F16.10)") rpm_miny
        call log_msg(msgstr)
        write(msgstr,"('rpm_numpart        = ',I0)") rpm_numpart
        call log_msg(msgstr)
        write(msgstr,"('rpm_addintervall        = ',I0)") rpm_addintervall
        call log_msg(msgstr)
        write(msgstr,"('rpm_deleteintervall        = ',I0)") rpm_deleteintervall
        call log_msg(msgstr)
        write(msgstr,"('rpm_minpartdist         = ',F16.10)") rpm_minpartdist
        call log_msg(msgstr)

        write(msgstr,"('n_msdx             = ',I0)") n_msdx
        call log_msg_md(msgstr)
        if (n_msdx/=0) then
           write (msgstr,"('n_msdx_sample      = ',I0)") n_msdx_sample
           call log_msg_md(msgstr)
           write (msgstr,"('msdx_len           = ',I0)") msdx_len
           call log_msg_md(msgstr)
           if (all(len_trim(msdx_name)==0)) call log_msg_md(&
                &'All mean square x-displacement regions are disabled.')
           do i=1,n_msdx_max
              if (len_trim(msdx_name(i))==0) cycle

              write (msgstr,&
                   &"('msdx region ""',A,'"" referred to as <md-x-',A,'>')"&
                   &) trim(msdx_name(i)),trim(msdx_name(i))
              call log_msg_md(msgstr)
              write (msgstr,"('  msdx_start          = ',I0)") msdx_start(i)
              call log_msg_md(msgstr)
              write (msgstr,"('  msdx_min[xyz]       =',3(X,F16.10))") &
                   &msdx_minx(i),msdx_miny(i),msdx_minz(i)
              call log_msg_md(msgstr)
              write (msgstr,"('  msdx_max[xyz]       =',3(X,F16.10))") &
                   &msdx_maxx(i),msdx_maxy(i),msdx_maxz(i)
              call log_msg_md(msgstr)
              if (any((/msdx_minx(i),msdx_miny(i),msdx_minz(i)&
                   &,msdx_maxx(i),msdx_maxy(i),msdx_maxz(i)/)<0.0_rk)) then
                 write (msgstr,"(' -> defaulting to')")
                 call log_msg_md(msgstr)

                 if (msdx_minx(i)<0.0_rk) msdx_minx(i) = minpos(1)
                 if (msdx_maxx(i)<0.0_rk) msdx_maxx(i) = maxpos(1)
                 if (msdx_miny(i)<0.0_rk) msdx_miny(i) = minpos(2)
                 if (msdx_maxy(i)<0.0_rk) msdx_maxy(i) = maxpos(2)
                 if (msdx_minz(i)<0.0_rk) msdx_minz(i) = minpos(3)
                 if (msdx_maxz(i)<0.0_rk) msdx_maxz(i) = maxpos(3)

                 write (msgstr,"('  msdx_min[xyz]       =',3(X,F16.10))") &
                      &msdx_minx(i),msdx_miny(i),msdx_minz(i)
                 call log_msg_md(msgstr)
                 write (msgstr,"('  msdx_max[xyz]       =',3(X,F16.10))") &
                      &msdx_maxx(i),msdx_maxy(i),msdx_maxz(i)
                 call log_msg_md(msgstr)
              end if
              write (msgstr,"('  msdx_subtract_drift = ',L1)") &
                   &msdx_subtract_drift(i)
              call log_msg_md(msgstr)
           end do
        end if

        write (msgstr,"('n_msdy             = ',I0)") n_msdy
        call log_msg_md(msgstr)
        if (n_msdy/=0) then
           write (msgstr,"('n_msdy_sample      = ',I0)") n_msdy_sample
           call log_msg_md(msgstr)
           write (msgstr,"('msdy_len           = ',I0)") msdy_len
           call log_msg_md(msgstr)
           if (all(len_trim(msdy_name)==0)) call log_msg_md(&
                &'All mean square y-displacement regions are disabled.')
           do i=1,n_msdy_max
              if (len_trim(msdy_name(i))==0) cycle

              write (msgstr,&
                   &"('msdy region ""',A,'"" referred to as <md-y-',A,'>')"&
                   &) trim(msdy_name(i)),trim(msdy_name(i))
              call log_msg_md(msgstr)
              write (msgstr,"('  msdy_start         = ',I0)") msdy_start(i)
              call log_msg_md(msgstr)
              write (msgstr,"('  msdy_min[xyz]      =',3(X,F16.10))") &
                   &msdy_minx(i),msdy_miny(i),msdy_minz(i)
              call log_msg_md(msgstr)
              write (msgstr,"('  msdy_max[xyz]      =',3(X,F16.10))") &
                   &msdy_maxx(i),msdy_maxy(i),msdy_maxz(i)
              call log_msg_md(msgstr)
              if (any((/msdy_minx(i),msdy_miny(i),msdy_minz(i)&
                   &,msdy_maxx(i),msdy_maxy(i),msdy_maxz(i)/)<0.0_rk)) then
                 write (msgstr,"(' -> defaulting to')")
                 call log_msg_md(msgstr)

                 if (msdy_minx(i)<0.0_rk) msdy_minx(i) = minpos(1)
                 if (msdy_maxx(i)<0.0_rk) msdy_maxx(i) = maxpos(1)
                 if (msdy_miny(i)<0.0_rk) msdy_miny(i) = minpos(2)
                 if (msdy_maxy(i)<0.0_rk) msdy_maxy(i) = maxpos(2)
                 if (msdy_minz(i)<0.0_rk) msdy_minz(i) = minpos(3)
                 if (msdy_maxz(i)<0.0_rk) msdy_maxz(i) = maxpos(3)

                 write (msgstr,"('  msdy_min[xyz]      =',3(X,F16.10))") &
                      &msdy_minx(i),msdy_miny(i),msdy_minz(i)
                 call log_msg_md(msgstr)
                 write (msgstr,"('  msdy_max[xyz]      =',3(X,F16.10))") &
                      &msdy_maxx(i),msdy_maxy(i),msdy_maxz(i)
                 call log_msg_md(msgstr)
              end if
              write (msgstr,"('  msdy_subtract_drift = ',L1)") &
                   &msdy_subtract_drift(i)
              call log_msg_md(msgstr)
           end do
        end if

        write (msgstr,"('n_msdz             = ',I0)") n_msdz
        call log_msg_md(msgstr)
        if (n_msdz/=0) then
           write (msgstr,"('n_msdz_sample      = ',I0)") n_msdz_sample
           call log_msg_md(msgstr)
           write (msgstr,"('msdz_len           = ',I0)") msdz_len
           call log_msg_md(msgstr)
           if (all(len_trim(msdz_name)==0)) call log_msg_md(&
                &'All mean square z-displacement regions are disabled.')
           do i=1,n_msdz_max
              if (len_trim(msdz_name(i))==0) cycle

              write (msgstr,&
                   &"('msdz region ""',A,'"" referred to as <md-z-',A,'>')"&
                   &) trim(msdz_name(i)),trim(msdz_name(i))
              call log_msg_md(msgstr)
              write (msgstr,"('  msdz_start         = ',I0)") msdz_start(i)
              call log_msg_md(msgstr)
              write (msgstr,"('  msdz_min[xyz]      =',3(X,F16.10))") &
                   &msdz_minx(i),msdz_miny(i),msdz_minz(i)
              call log_msg_md(msgstr)
              write (msgstr,"('  msdz_max[xyz]      =',3(X,F16.10))") &
                   &msdz_maxx(i),msdz_maxy(i),msdz_maxz(i)
              call log_msg_md(msgstr)
              if (any((/msdz_minx(i),msdz_miny(i),msdz_minz(i)&
                   &,msdz_maxx(i),msdz_maxy(i),msdz_maxz(i)/)<0.0_rk)) then
                 write (msgstr,"(' -> defaulting to')")
                 call log_msg_md(msgstr)

                 if (msdz_minx(i)<0.0_rk) msdz_minx(i) = minpos(1)
                 if (msdz_maxx(i)<0.0_rk) msdz_maxx(i) = maxpos(1)
                 if (msdz_miny(i)<0.0_rk) msdz_miny(i) = minpos(2)
                 if (msdz_maxy(i)<0.0_rk) msdz_maxy(i) = maxpos(2)
                 if (msdz_minz(i)<0.0_rk) msdz_minz(i) = minpos(3)
                 if (msdz_maxz(i)<0.0_rk) msdz_maxz(i) = maxpos(3)

                 write (msgstr,"('  msdz_min[xyz]      =',3(X,F16.10))") &
                      &msdz_minx(i),msdz_miny(i),msdz_minz(i)
                 call log_msg_md(msgstr)
                 write (msgstr,"('  msdz_max[xyz]      =',3(X,F16.10))") &
                      &msdz_maxx(i),msdz_maxy(i),msdz_maxz(i)
                 call log_msg_md(msgstr)
              end if
              write (msgstr,"('  msdz_subtract_drift = ',L1)") &
                   &msdz_subtract_drift(i)
              call log_msg_md(msgstr)
           end do
        end if

        call log_ws()

        call MPI_Bcast(n_dump,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(dump_format,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(dump_double,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_ids,1,MPI_LOGICAL,0,comm_cart,ierror)

        call MPI_Bcast(dump_time,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_positions,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_potentials,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_pnd,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_velocities,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_forces,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_quaternions,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_orientations,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_rotations,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_radii,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_summary,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_torques,1,MPI_LOGICAL,0,comm_cart,ierror)
#ifdef PARTICLESTRESS
        call MPI_Bcast(dump_stresses,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_stresses_b,1,MPI_LOGICAL,0,comm_cart,ierror)
#endif
#ifdef LADD_SURR_RHOF
        call MPI_Bcast(dump_surr_rhof,1,MPI_LOGICAL,0,comm_cart,ierror)
#endif
!#ifdef MD_MAG
       call MPI_Bcast(dump_mag,1,MPI_LOGICAL,0,comm_cart,ierror)
!#endif

        call MPI_Bcast(n_dump_substep,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_dump_substep_start,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_stat,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(count_periodic,1,MPI_LOGICAL,0,comm_cart,ierror)

        call MPI_Bcast(n_rvz_profile_dump,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_rvz_profile_sample,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_dfz,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(dfz_minx,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(dfz_maxx,1,MPI_INTEGER,0,comm_cart,ierror)

        call MPI_Bcast(interaction,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(potential,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(rock,32,MPI_CHARACTER,0,comm_cart,ierror)
		  call MPI_Bcast(rock_friction,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(growing_stage,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(md_magnetic,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(lyapunov,1,MPI_LOGICAL,0,comm_cart,ierror)

        call MPI_Bcast(no_md_leesedwards,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(md_leesedwards_reproducible_cp,1,MPI_LOGICAL,0,comm_cart&
             &,ierror)

        call MPI_Bcast(initial_placing,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(init_file,1024,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(alat,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(phi,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(min_dist,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(np_sphere,1,MPI_INTEGER,0,Comm_cart,ierror)
        call MPI_Bcast(np_sphere_cap,1,MPI_INTEGER,0,Comm_cart,ierror)
        call MPI_Bcast(x0lo,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(x0hi,3,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(initial_velocities,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(v0,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(temp0,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(initial_orientations,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(q0,4,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(initial_rotations,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(w0,3,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(new_uids,32,MPI_CHARACTER,0,comm_cart,ierror)

        call MPI_Bcast(delta_x,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(delta_t,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(delta_m,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(steps_per_lbe_step,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(average_ft_fluid,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(average_vw_fluid,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(rho,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(mass,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(inertia_orth,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(inertia_para,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(ineigh,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(time_output,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(prtcl_output,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(rc,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(rs,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(semiper_max,3,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(semiper_min,3,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(diffuse_x,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(diffuse_y,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(diffuse_z,1,MPI_LOGICAL,0,comm_cart,ierror)

        call MPI_Bcast(reenter_rdm_min,3,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(reenter_rdm_max,3,MPI_LOGICAL,0,comm_cart,ierror)

        call MPI_Bcast(cutoff_v,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(reflective_rocks,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(molec_mass,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(molec_mass_lbm,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(temperat,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(mean_free_path,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(fix_v,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(fix_w,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(fix_vx,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(fix_vx_t,1,MPI_INTEGER,0,comm_cart,ierror)
       call MPI_Bcast(fix_vz,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(fix_vz_t,1,MPI_INTEGER,0,comm_cart,ierror)

        call MPI_Bcast(ramping,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(vmax_ramping,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(ramptime,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(forcebalance,1,MPI_LOGICAL,0,comm_cart,ierror)

        call MPI_Bcast(constant_f,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(cargo_par,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(cargo_par_id,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(constant_fc,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(constant_t,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(constant_tc,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(constant_f_y_min,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(gaussian_plane_active,1,MPI_LOGICAL,0,comm_cart,ierror)


      
        call MPI_Bcast(gaussian_plane_active,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(gaussian_plane_E,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(gaussian_plane_c,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(gaussian_plane_sqwidth,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(gaussian_plane_pos,1,MPI_REAL8,0,comm_cart,ierror)
		  call MPI_Bcast(randompartforce,1,MPI_LOGICAL,0,comm_cart,ierror)

        call MPI_Bcast(runtime_particle_mod,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(rpm_miny,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(rpm_numpart,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(rpm_addintervall,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(rpm_deleteintervall,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(rpm_minpartdist,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(boundary_width,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(n_msdx,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_msdx_sample,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdx_len,1,MPI_INTEGER,0,comm_cart,ierror)

        do i=1,n_msdx_max
           call MPI_Bcast(msdx_name(i),16,MPI_CHARACTER,0,comm_cart,ierror)
        end do
        call MPI_Bcast(msdx_start,n_msdx_max,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdx_minx,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_maxx,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_miny,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_maxy,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_minz,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_maxz,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_subtract_drift,n_msdx_max,MPI_LOGICAL,0,comm_cart&
             &,ierror)

        call MPI_Bcast(n_msdy,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_msdy_sample,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdy_len,1,MPI_INTEGER,0,comm_cart,ierror)

        do i=1,n_msdy_max
           call MPI_Bcast(msdy_name(i),16,MPI_CHARACTER,0,comm_cart,ierror)
        end do
        call MPI_Bcast(msdy_start,n_msdy_max,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdy_minx,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_maxx,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_miny,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_maxy,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_minz,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_maxz,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_subtract_drift,n_msdy_max,MPI_LOGICAL,0,comm_cart&
             &,ierror)

        call MPI_Bcast(n_msdz,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_msdz_sample,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdz_len,1,MPI_INTEGER,0,comm_cart,ierror)

        do i=1,n_msdz_max
           call MPI_Bcast(msdz_name(i),16,MPI_CHARACTER,0,comm_cart,ierror)
        end do
        call MPI_Bcast(msdz_start,n_msdz_max,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdz_minx,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_maxx,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_miny,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_maxy,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_minz,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_maxz,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_subtract_drift,n_msdz_max,MPI_LOGICAL,0,comm_cart&
             &,ierror)

        if (initial_placing=='file_parallel:asc/xvowr') polydispersity = .true.
        if (initial_placing=='file_parallel:asc/xvowrm') then 
        polydispersity = .true.
        magdispersity  = .true.
        end if
        if (initial_placing=='file_parallel:asc/xvowm') then
          magdispersity  = .true.
        end if
        if (sci_profile.and.rc<1.5) call error_md('sci_profile requires rc>1.5')

        return
100     continue
        call error_md('Error opening md input file "'//trim(inp_file)//'.md"')
        return
110     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
        return
120     continue
        call error_md('Error closing md input file "'//trim(inp_file)//'.md"')
        return
    end subroutine input

    !> restores from an md checkpoint file
    subroutine md_restore_checkpoint
        type(md_particle_type),allocatable,dimension(:) :: tp
        character(len=1024) chk_file_name
        integer i,chk_file_handle,ierror,n_global,stat
        real(kind=rk) :: pf(3),tmp(3)

        ! This would be a function pointer in C. In fact,  xdrfdouble  is
        ! an external function, but fortran does not seem to care about the
        ! type anyway...
        integer,external :: xdrfdouble

#ifdef USEXDRF
        rank_0: if (myrankc==0) then
           call lbe_make_filename_restore(chk_file_name,'md-checkpoint','.xdr')

           call xdrfopen(chk_file_handle,trim(chk_file_name),'r',ierror)
           if (ierror==0) then
              ! Checking for XDRSLOPPY probably makes no sense here,
              ! the checkpoint MUST be read.
              call error_md(&
                   &'md_restore_checkpoint(): xdrfopen() failed opening "'&
                   &//trim(chk_file_name)//'"')
           end if

           ! ...at last, read particle data:
           call xdrfint(chk_file_handle,n_global,ierror)

           allocate (tp(n_global),stat=stat)
           call check_allocate(stat,'md_restore_checkpoint(): tp')


           do i=1,n_global
              call xdrfvector(chk_file_handle,tp(i)%x(1:3),3,xdrfdouble,ierror)
              call xdrfvector(chk_file_handle,tp(i)%v(1:3),3,xdrfdouble,ierror)
              if (use_rotation) then
                 call xdrfvector(chk_file_handle,tp(i)%q(0:3),4,xdrfdouble&
                      &,ierror)
                 call xdrfvector(chk_file_handle,tp(i)%w(1:3),3,xdrfdouble&
                      &,ierror)
              end if
              call xdrfint(chk_file_handle,tp(i)%uid,ierror)
              call xdrfvector(chk_file_handle,tp(i)%vnew(1:3),3,xdrfdouble&
                   &,ierror)
              if (use_rotation) then
                 call xdrfvector(chk_file_handle,tp(i)%qnew(0:3),4,xdrfdouble&
                      &,ierror)
                 call xdrfvector(chk_file_handle,tp(i)%wnew(1:3),3,xdrfdouble&
                      &,ierror)
              end if

              call xdrfint(chk_file_handle,tp(i)%master,ierror)

              if (polydispersity) then
                 call xdrfvector(chk_file_handle,tp(i)%R_orth,1,xdrfdouble&
                      &,ierror)
                 call xdrfvector(chk_file_handle,tp(i)%R_para,1,xdrfdouble&
                      &,ierror)
              end if

               if (magdispersity) then
                   call xdrfvector(chk_file_handle,tp(i)%mag,1,xdrfdouble&
                        &,ierror)
                end if

              call xdrfvector(chk_file_handle,tmp,3,xdrfdouble,ierror)
              tp(i)%v_fluid_acc = tmp*real(steps_per_lbe_step-1,kind=rk)
              if (use_rotation) then
                 call xdrfvector(chk_file_handle,tmp,3,xdrfdouble,ierror)
                 tp(i)%ws_fluid_acc = tmp*real(steps_per_lbe_step-1,kind=rk)
              end if

              call xdrfvector(chk_file_handle,tp(i)%v_fluid(1:3),3,xdrfdouble&
                   &,ierror)
              if (use_rotation) then
                 call xdrfvector(chk_file_handle,tp(i)%ws_fluid(1:3),3&
                      &,xdrfdouble,ierror)
              end if

              call xdrfvector(chk_file_handle,tp(i)%f_fluid(1:3),3,xdrfdouble&
                   &,ierror)
              if (use_rotation) then
                 call xdrfvector(chk_file_handle,tp(i)%t_fluid(1:3),3&
                      &,xdrfdouble,ierror)
              end if

#ifdef LADD_SURR_RHOF
              call xdrfvector(chk_file_handle,tp(i)%rhof,1,xdrfdouble,ierror)
#endif
           end do

           ! global variables used only by ladd code
           ladd_r: if (interaction=='ladd') then
              call xdrfvector(chk_file_handle,pf,3,xdrfdouble,ierror)
              pfr = pf(1)
              pfb = pf(2)
              pfg = pf(3)
              call xdrfvector(chk_file_handle,global_mass_change,n_spec,xdrfdouble,ierror)
              call xdrfvector(chk_file_handle,global_mass_target,n_spec,xdrfdouble,ierror)
           end if ladd_r

           periodic_inflow_r: if (boundary=='periodic_inflow') then
              call xdrfint(chk_file_handle,first_inserted_puid,ierror)
           end if periodic_inflow_r

           call xdrfclose(chk_file_handle,ierror)

        end if rank_0

        ladd_b: if (interaction=='ladd') then
           call MPI_Bcast(pfr,1,MPI_REAL8,0,comm_cart,ierror)
           call MPI_Bcast(pfb,1,MPI_REAL8,0,comm_cart,ierror)
           call MPI_Bcast(pfg,1,MPI_REAL8,0,comm_cart,ierror)
           call MPI_Bcast(global_mass_target,n_spec,MPI_REAL8,0,comm_cart,ierror)
           call MPI_Bcast(global_mass_change,n_spec,MPI_REAL8,0,comm_cart,ierror)
           call log_msg_md("Updating global mass change through checkpoint restoration ...")
           do i=1,n_spec
             write(msgstr,"('  global_mass_target(',I0,') = ',E16.8)") i, global_mass_target(i)
             call log_msg_md(msgstr)
             write(msgstr,"('  global_mass_change(',I0,') = ',E16.8, F16.10)") i, global_mass_change(i), 1.0 + global_mass_change(i)/global_mass_target(i)
             call log_msg_md(msgstr)
           end do
        end if ladd_b

        periodic_inflow_b: if (boundary=='periodic_inflow') then
           call MPI_Bcast(first_inserted_puid,1,MPI_INTEGER,0,comm_cart,ierror)
        end if periodic_inflow_b

        call setup_list
        call scatter_particles(n_global,tp,particle_cp_mpitype)

        if (myrankc==0) then
           deallocate (tp)
           write(msgstr&
                &,"('Placed a total of ',I0,' particles from checkpoint.')") &
                &n_global
           call log_msg_md(msgstr)
        end if
#else
        call error_md('failed to restore from checkpoint file'&
             &//' - recompile with USEXDRF set!')
#endif
    end subroutine md_restore_checkpoint

#endif
end module lbe_md_input_module
