# Makefile created by mkmf.pl $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include /nfs/home/6/sanjeevis/lb3d/src/Makefile.template


.DEFAULT:
	-echo $@ does not exist.
all: version lbe
fft_3d.o: ./fft_3d.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	./fft_3d.c
hoshen_kopelman.o: ./hoshen_kopelman.F90 lbe.h lbe_globals.o lbe_log.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./hoshen_kopelman.F90
lbe.o: ./lbe.F90 lbe.h lbe_bc.o lbe_collision.o lbe_force.o lbe_galilean_stabilizer.o lbe_globals.o lbe_helper.o lbe_init.o lbe_io_interface.o lbe_io.o lbe_io_checkpoint.o lbe_io_stress.o lbe_log.o lbe_mass_scaler.o lbe_mean_square_displacement.o lbe_parms.o lbe_parallel.o lbe_timer.o lbe_types.o lbe_elec_interface.o lbe_md_interface.o lbe_collision_simple.o lsuperstruct.o lbe_tracer_interface.o lbe_io_hdf5.o lbe_axisym.o lbe_invasion.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe.F90
lbe_MRT.o: ./lbe_MRT.F90 lbe.h lbe_globals.o lbe_parms.o luxury_rng.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_MRT.F90
lbe_analysis.o: ./lbe_analysis.F90 lbe.h lbe_globals.o lbe_helper.o lbe_parms.o lbe_types.o lbe_md_analysis.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_analysis.F90
lbe_axisym.o: ./lbe_axisym.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_parms.o lbe_parallel.o lbe_timer.o lbe_types.o lsuperstruct_data.o lsuperstruct_interface.o lbe_force_interface.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_axisym.F90
lbe_bc.o: ./lbe_bc.F90 lbe.h lbe_helper.o lbe_globals.o lbe_log.o lbe_invasion.o lbe_leesedwards.o lbe_md_bc_leesedwards.o lbe_parallel.o lbe_parms.o lbe_timer.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_bc.F90
lbe_bdist.o: ./lbe_bdist.F90 lbe.h lbe_globals.o lbe_parallel.o lbe_parms.o luxury_rng.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_bdist.F90
lbe_collision.o: ./lbe_collision.F90 lbe.h lbe_bdist.o lbe_globals.o lbe_helper.o lbe_log.o lbe_leesedwards.o lbe_MRT.o lbe_parms.o lbe_parallel.o lbe_timer.o lbe_types.o lsuperstruct_data.o lsuperstruct_interface.o lbe_force_interface.o lbe_axisym.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_collision.F90
lbe_collision_simple.o: ./lbe_collision_simple.F90 lbe.h lbe_bdist.o lbe_collision.o lbe_globals.o lbe_parallel.o lbe_parms.o lbe_timer.o lbe_types.o lbe_helper.o lbe_md.o lbe_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_collision_simple.F90
lbe_elec_fluxes.o: ./lbe_elec_fluxes.F90 lbe.h lbe_elec_forces.o lbe_elec_globals.o lbe_elec_helper.o lbe_elec_timer.o lbe_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_fluxes.F90
lbe_elec_forces.o: ./lbe_elec_forces.F90 lbe.h lbe_elec_helper.o lbe_elec_globals.o lbe_force_interface.o lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_forces.F90
lbe_elec_globals.o: ./lbe_elec_globals.F90 lbe.h lbe_globals.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_globals.F90
lbe_elec_helper.o: ./lbe_elec_helper.F90 lbe.h lbe_elec_globals.o lbe_helper.o lbe_log.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_helper.F90
lbe_elec_init.o: ./lbe_elec_init.F90 lbe.h lbe_elec_fluxes.o lbe_elec_helper.o lbe_elec_globals.o lbe_elec_output.o lbe_elec_parallel.o lbe_elec_poisson_solver.o lbe_elec_timer.o lbe_elec_timestep.o lbe_globals.o lbe_helper.o lbe_md_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_init.F90
lbe_elec_input.o: ./lbe_elec_input.F90 lbe.h lbe_elec_globals.o lbe_elec_helper.o lbe_globals.o lbe_log.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_input.F90
lbe_elec_interface.o: ./lbe_elec_interface.F90 lbe.h lbe_elec_timestep.o lbe_elec_forces.o lbe_elec_input.o lbe_elec_init.o lbe_elec_output.o lbe_elec_parallel.o lbe_elec_timer.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_interface.F90
lbe_elec_output.o: ./lbe_elec_output.F90 lbe.h lbe_elec_globals.o lbe_elec_helper.o lbe_elec_timer.o lbe_parms.o lbe_io_helper.o lbe_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_output.F90
lbe_elec_parallel.o: ./lbe_elec_parallel.F90 lbe.h lbe_elec_globals.o lbe_elec_helper.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_parallel.F90
lbe_elec_poisson_solver.o: ./lbe_elec_poisson_solver.F90 lbe.h lbe_elec_globals.o lbe_elec_helper.o lbe_elec_parallel.o lbe_elec_timer.o lbe_parallel.o pppm_global.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_poisson_solver.F90
lbe_elec_timer.o: ./lbe_elec_timer.F90 lbe.h lbe_elec_globals.o lbe_elec_helper.o lbe_parallel.o lbe_timer.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_timer.F90
lbe_elec_timestep.o: ./lbe_elec_timestep.F90 lbe.h lbe_elec_fluxes.o lbe_elec_globals.o lbe_elec_poisson_solver.o lbe_elec_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_elec_timestep.F90
lbe_force.o: ./lbe_force.F90 lbe.h lbe_force_constant.o lbe_force_kolmogorov.o lbe_globals.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_force.F90
lbe_force_constant.o: ./lbe_force_constant.F90 lbe.h lbe_force_interface.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_fluid.o lbe_md_globals.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_force_constant.F90
lbe_force_interface.o: ./lbe_force_interface.F90 lbe.h lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_force_interface.F90
lbe_force_kolmogorov.o: ./lbe_force_kolmogorov.F90 lbe.h lbe_force_interface.o lbe_globals.o lbe_log.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_force_kolmogorov.F90
lbe_galilean_stabilizer.o: ./lbe_galilean_stabilizer.F90 lbe.h lbe_analysis.o lbe_bdist.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_interface.o lbe_parallel.o lbe_parms.o lbe_tracer_interface.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_galilean_stabilizer.F90
lbe_globals.o: ./lbe_globals.F90 lbe.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_globals.F90
lbe_helper.o: ./lbe_helper.F90 lbe.h lbe_globals.o lbe_log.o lbe_md_globals.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_helper.F90
lbe_init.o: ./lbe_init.F90 lbe.h lbe_globals.o lbe_helper.o lbe_io_checkpoint.o lbe_io.o lbe_io_xdrf.o lbe_log.o lbe_md_globals.o lbe_parallel.o lbe_parms.o lbe_timer.o lbe_bdist.o lbe_types.o lbe_init_rock.o luxury_rng.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_init.F90
lbe_init_rock.o: ./lbe_init_rock.F90 lbe.h lbe_globals.o lbe_io.o lbe_io_hdf5.o lbe_io_xdrf.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o lbe_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_init_rock.F90
lbe_invasion.o: ./lbe_invasion.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_bdist.o lbe_collision.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_invasion.F90
lbe_io.o: ./lbe_io.F90 lbe.h hoshen_kopelman.o lbe_analysis.o lbe_collision.o lbe_globals.o lbe_helper.o lbe_io_helper.o lbe_log.o lbe_parms.o lbe_parallel.o lbe_types.o lbe_io_hdf5.o lbe_io_xdrf.o lbe_md_globals.o map.o lbe_md_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io.F90
lbe_io_arrstats.o: ./lbe_io_arrstats.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_io_hdf5.o lbe_io_helper.o lbe_io.o lbe_io_xdrf.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io_arrstats.F90
lbe_io_checkpoint.o: ./lbe_io_checkpoint.F90 lbe.h lbe_globals.o lbe_io_arrstats.o lbe_mean_square_displacement.o lbe_log.o lbe_parms.o lbe_parallel.o lbe_helper.o lbe_io_helper.o lbe_leesedwards.o lbe_timer.o lbe_types.o lbe_io_xdrf.o lbe_io_hdf5.o lsuperstruct.o lbe_md_globals.o lbe_md_output.o lbe_tracer_output.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io_checkpoint.F90
lbe_io_hdf5.o: ./lbe_io_hdf5.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_io_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io_hdf5.F90
lbe_io_helper.o: ./lbe_io_helper.F90 lbe.h lbe_log.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io_helper.F90
lbe_io_interface.o: ./lbe_io_interface.F90 lbe.h lbe_bc.o lbe_globals.o lbe_helper.o lbe_io_arrstats.o lbe_io_hdf5.o lbe_io.o lbe_io_stress.o lbe_log.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io_interface.F90
lbe_io_stress.o: ./lbe_io_stress.F90 lbe.h lbe_parms.o lbe_globals.o lbe_log.o lbe_parallel.o lbe_helper.o lbe_io_helper.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io_stress.F90
lbe_io_xdrf.o: ./lbe_io_xdrf.F90 lbe.h lbe_globals.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_io_helper.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_io_xdrf.F90
lbe_leesedwards.o: ./lbe_leesedwards.F90 lbe.h lbe_bdist.o lbe_globals.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_leesedwards.F90
lbe_log.o: ./lbe_log.F90 lbe.h lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_log.F90
lbe_mass_scaler.o: ./lbe_mass_scaler.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_md_interface.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_mass_scaler.F90
lbe_md.o: ./lbe_md.F90 lbe.h lbe_parms.o lbe_bdist.o lbe_helper.o lbe_globals.o lbe_log.o lbe_md_dynamic.o lbe_md_bc_leesedwards.o lbe_md_boundary_condition.o lbe_md_fluid.o lbe_md_globals.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_helper.o lbe_md_memory.o lbe_md_output.o lbe_md_parallel.o lbe_md_potential.o lbe_md_magnetic.o lbe_md_rand_walk.o lbe_md_rock.o lbe_timer.o lbe_parallel.o lbe_types.o luxury_rng.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md.F90
lbe_md_analysis.o: ./lbe_md_analysis.F90 lbe.h lbe_globals.o lbe_md_globals.o lbe_md_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_analysis.F90
lbe_md_bc_leesedwards.o: ./lbe_md_bc_leesedwards.F90 lbe.h lbe_bdist.o lbe_globals.o lbe_helper.o lbe_leesedwards.o lbe_md_fluid_ladd.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_types.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_bc_leesedwards.F90
lbe_md_boundary_condition.o: ./lbe_md_boundary_condition.F90 lbe.h lbe_globals.o lbe_md_dynamic.o lbe_md_globals.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_boundary_condition.F90
lbe_md_debug.o: ./lbe_md_debug.F90 lbe.h lbe_globals.o lbe_io_xdrf.o lbe_md_globals.o lbe_parallel.o lbe_parms.o lbe_types.o lbe_md_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_debug.F90
lbe_md_dynamic.o: ./lbe_md_dynamic.F90 lbe.h lbe_globals.o lbe_log.o lbe_md_globals.o lbe_md_helper.o lbe_md_memory.o lbe_parallel.o lbe_parms.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_dynamic.F90
lbe_md_fluid.o: ./lbe_md_fluid.F90 lbe.h lbe_globals.o lbe_helper.o lbe_md_fluid_friction.o lbe_md_fluid_ladd.o lbe_md_fluid_tracer.o lbe_md_globals.o lbe_md_helper.o lbe_types.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid.F90
lbe_md_fluid_friction.o: ./lbe_md_fluid_friction.F90 lbe.h lbe_force_interface.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid_friction.F90
lbe_md_fluid_ladd.o: ./lbe_md_fluid_ladd.F90 lbe.h lbe_helper.o lbe_bdist.o lbe_force_interface.o lbe_globals.o lbe_leesedwards.o lbe_md_boundary_condition.o lbe_md_dynamic.o lbe_md_fluid_ladd_club.o lbe_md_fluid_ladd_dlub.o lbe_md_fluid_ladd_mc.o lbe_md_fluid_ladd_parms.o lbe_log.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_timer.o lbe_types.o lbe_collision.o luxury_rng.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid_ladd.F90
lbe_md_fluid_ladd_club.o: ./lbe_md_fluid_ladd_club.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_md_fluid_ladd_mc.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid_ladd_club.F90
lbe_md_fluid_ladd_dlub.o: ./lbe_md_fluid_ladd_dlub.F90 lbe.h lbe_globals.o lbe_helper.o lbe_md_boundary_condition.o lbe_md_fluid_ladd_mc.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_timer.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid_ladd_dlub.F90
lbe_md_fluid_ladd_mc.o: ./lbe_md_fluid_ladd_mc.F90 lbe.h lbe_globals.o lbe_log.o lbe_md_globals.o lbe_md_helper.o lbe_helper.o lbe_parallel.o lbe_parms.o lbe_types.o lbe_bdist.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid_ladd_mc.F90
lbe_md_fluid_ladd_parms.o: ./lbe_md_fluid_ladd_parms.F90 lbe.h lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid_ladd_parms.F90
lbe_md_fluid_tracer.o: ./lbe_md_fluid_tracer.F90 lbe.h lbe_globals.o lbe_md_globals.o lbe_md_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_fluid_tracer.F90
lbe_md_globals.o: ./lbe_md_globals.F90 lbe.h lbe_globals.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_globals.F90
lbe_md_growing_stage.o: ./lbe_md_growing_stage.F90 lbe.h lbe_globals.o lbe_init_rock.o lbe_log.o lbe_md_globals.o lbe_md_helper.o lbe_md.o lbe_md_output.o lbe_md_potential.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_growing_stage.F90
lbe_md_helper.o: ./lbe_md_helper.F90 lbe.h lbe_helper.o lbe_globals.o lbe_log.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_parallel.o lbe_parms.o lbe_io_helper.o lbe_types.o lbe_collision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_helper.F90
lbe_md_init.o: ./lbe_md_init.F90 lbe.h lbe_globals.o lbe_helper.o lbe_init_rock.o lbe_log.o lbe_md_bc_leesedwards.o lbe_md_boundary_condition.o lbe_md_debug.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_md_memory.o lbe_md.o lbe_md_output.o lbe_md_parallel.o lbe_md_rand_walk.o lbe_md_rock.o lbe_parallel.o lbe_parms.o lbe_timer.o lbe_types.o luxury_rng.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_init.F90
lbe_md_input.o: ./lbe_md_input.F90 lbe.h lbe_globals.o lbe_log.o lbe_md_bc_leesedwards.o lbe_md_dynamic.o lbe_md_fluid_ladd_mc.o lbe_md_globals.o lbe_md_growing_stage.o lbe_md_helper.o lbe_md_init.o lbe_md_lyapunov.o lbe_md.o lbe_md_parallel.o lbe_md_output.o lbe_md_potential.o lbe_md_magnetic.o lbe_md_rock.o lbe_parallel.o lbe_parms.o lbe_io_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_input.F90
lbe_md_interface.o: ./lbe_md_interface.F90 lbe.h lbe_bc.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_bc_leesedwards.o lbe_md_dynamic.o lbe_md_fluid.o lbe_md_globals.o lbe_md_growing_stage.o lbe_md_helper.o lbe_md_init.o lbe_md_input.o lbe_md_lyapunov.o lbe_md_memory.o lbe_md.o lbe_md_output.o lbe_md_potential.o lbe_md_rock.o lbe_timer.o lbe_parallel.o lbe_parms.o lbe_types.o lbe_md_magnetic.o lbe_md_fluid_ladd_mc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_interface.F90
lbe_md_lyapunov.o: ./lbe_md_lyapunov.F90 lbe.h lbe_globals.o lbe_init_rock.o lbe_io_helper.o lbe_log.o lbe_md_globals.o lbe_md_helper.o lbe_md.o lbe_md_parallel.o lbe_timer.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_lyapunov.F90
lbe_md_magnetic.o: ./lbe_md_magnetic.F90 lbe.h lbe_globals.o lbe_helper.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_log.o lbe_md_helper.o lbe_md_output.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_magnetic.F90
lbe_md_memory.o: ./lbe_md_memory.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_memory.F90
lbe_md_output.o: ./lbe_md_output.F90 lbe.h lbe_globals.o lbe_helper.o lbe_md_dynamic.o lbe_md_fluid.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_mc.o lbe_md_globals.o lbe_md_helper.o lbe_md_parallel.o lbe_mean_square_displacement.o lbe_timer.o lbe_parallel.o lbe_parms.o lbe_types.o lbe_io_helper.o lbe_io_xdrf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_output.F90
lbe_md_parallel.o: ./lbe_md_parallel.F90 lbe.h lbe_globals.o lbe_md_bc_leesedwards.o lbe_md_globals.o lbe_md_helper.o lbe_md_memory.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_parallel.F90
lbe_md_potential.o: ./lbe_md_potential.F90 lbe.h lbe_md_globals.o lbe_md_helper.o lbe_md_potential_bprh.o lbe_md_potential_bprhag.o lbe_md_potential_gb.o lbe_md_potential_hertz.o lbe_md_potential_lj.o lbe_md_potential_dlvo.o lbe_md_potential_spring.o lbe_md_potential_lubrication.o lbe_md_fluid_ladd_parms.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential.F90
lbe_md_potential_bprh.o: ./lbe_md_potential_bprh.F90 lbe.h lbe_helper.o lbe_globals.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_md_output.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_bprh.F90
lbe_md_potential_bprhag.o: ./lbe_md_potential_bprhag.F90 lbe.h lbe_helper.o lbe_globals.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_bprhag.F90
lbe_md_potential_dlvo.o: ./lbe_md_potential_dlvo.F90 lbe.h lbe_globals.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_dlvo.F90
lbe_md_potential_gb.o: ./lbe_md_potential_gb.F90 lbe.h lbe_helper.o lbe_globals.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_gb.F90
lbe_md_potential_hertz.o: ./lbe_md_potential_hertz.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_hertz.F90
lbe_md_potential_lj.o: ./lbe_md_potential_lj.F90 lbe.h lbe_globals.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_lj.F90
lbe_md_potential_lubrication.o: ./lbe_md_potential_lubrication.F90 lbe.h lbe_globals.o lbe_md_fluid_ladd.o lbe_md_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_lubrication.F90
lbe_md_potential_spring.o: ./lbe_md_potential_spring.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_potential_spring.F90
lbe_md_rand_walk.o: ./lbe_md_rand_walk.F90 lbe.h lbe_helper.o lbe_globals.o lbe_md_globals.o lbe_md_helper.o lbe_types.o luxury_rng.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_rand_walk.F90
lbe_md_rock.o: ./lbe_md_rock.F90 lbe.h lbe_globals.o lbe_helper.o lbe_md_globals.o lbe_md_helper.o lbe_md_rock_lj.o lbe_md_rock_dlvo.o lbe_md_rock_bprh.o lbe_md_rock_friction.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_rock.F90
lbe_md_rock_bprh.o: ./lbe_md_rock_bprh.F90 lbe.h lbe_bc.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_rock_bprh.F90
lbe_md_rock_dlvo.o: ./lbe_md_rock_dlvo.F90 lbe.h lbe_bc.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_rock_dlvo.F90
lbe_md_rock_friction.o: ./lbe_md_rock_friction.F90 lbe.h lbe_bc.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_rock_friction.F90
lbe_md_rock_lj.o: ./lbe_md_rock_lj.F90 lbe.h lbe_bc.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_md_helper.o lbe_parallel.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_rock_lj.F90
lbe_md_scale_vel.o: ./lbe_md_scale_vel.F90 lbe.h lbe_md_globals.o lbe_md_helper.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_md_scale_vel.F90
lbe_mean_square_displacement.o: ./lbe_mean_square_displacement.F90 lbe.h lbe_analysis.o lbe_globals.o lbe_helper.o lbe_io_helper.o lbe_io_xdrf.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o map.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_mean_square_displacement.F90
lbe_parallel.o: ./lbe_parallel.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_parms.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_parallel.F90
lbe_parms.o: ./lbe_parms.F90 lbe.h lbe_version.h lbe_globals.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_parms.F90
lbe_timer.o: ./lbe_timer.F90 lbe.h lbe_globals.o lbe_log.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_timer.F90
lbe_tracer.o: ./lbe_tracer.F90 lbe.h lbe_globals.o lbe_parallel.o lbe_timer.o lbe_tracer_globals.o lbe_tracer_helper.o lbe_tracer_memory.o lbe_tracer_output.o lbe_tracer_recoloring.o lbe_types.o lbe_md_globals.o lbe_tracer_md_fluid_ladd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer.F90
lbe_tracer_globals.o: ./lbe_tracer_globals.F90 lbe.h lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_globals.F90
lbe_tracer_helper.o: ./lbe_tracer_helper.F90 lbe.h lbe_collision.o lbe_globals.o lbe_helper.o lbe_log.o lbe_md_fluid.o lbe_parallel.o lbe_parms.o lbe_tracer_globals.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_helper.F90
lbe_tracer_init.o: ./lbe_tracer_init.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_timer.o lbe_tracer_globals.o lbe_tracer_helper.o lbe_tracer_memory.o lbe_tracer.o lbe_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_init.F90
lbe_tracer_input.o: ./lbe_tracer_input.F90 lbe.h lbe_globals.o lbe_io_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_tracer_globals.o lbe_tracer_helper.o lbe_tracer_init.o lbe_tracer.o lbe_tracer_parallel.o lbe_tracer_output.o lbe_tracer_recoloring.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_input.F90
lbe_tracer_interface.o: ./lbe_tracer_interface.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_parms.o lbe_timer.o lbe_tracer_globals.o lbe_tracer_helper.o lbe_tracer_init.o lbe_tracer_input.o lbe_tracer_memory.o lbe_tracer.o lbe_tracer_output.o lbe_types.o lbe_md_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_interface.F90
lbe_tracer_md_fluid_ladd.o: ./lbe_tracer_md_fluid_ladd.F90 lbe.h lbe_helper.o lbe_md_boundary_condition.o lbe_md_fluid_ladd.o lbe_md_fluid_ladd_parms.o lbe_md_globals.o lbe_tracer_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_md_fluid_ladd.F90
lbe_tracer_memory.o: ./lbe_tracer_memory.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_tracer_globals.o lbe_tracer_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_memory.F90
lbe_tracer_output.o: ./lbe_tracer_output.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o lbe_io_helper.o lbe_mean_square_displacement.o lbe_parallel.o lbe_parms.o lbe_tracer_globals.o lbe_tracer_helper.o lbe_types.o lbe_io_xdrf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_output.F90
lbe_tracer_parallel.o: ./lbe_tracer_parallel.F90 lbe.h lbe_globals.o lbe_parallel.o lbe_tracer_globals.o lbe_tracer_helper.o lbe_tracer_memory.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_parallel.F90
lbe_tracer_recoloring.o: ./lbe_tracer_recoloring.F90 lbe.h lbe_tracer_globals.o lbe_tracer_helper.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_tracer_recoloring.F90
lbe_types.o: ./lbe_types.F90 lbe.h lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lbe_types.F90
lextobj.o: ./lextobj.F90 lbe.h lbe_globals.o lbe_parallel.o lbe_helper.o lbe_log.o lmesh.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lextobj.F90
lmesh.o: ./lmesh.F90 lbe.h lbe_globals.o lbe_helper.o lbe_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lmesh.F90
lsuperstruct.o: ./lsuperstruct.F90 lbe.h lbe_globals.o lbe_parms.o lbe_types.o lbe_timer.o lsuperstruct_IBM.o lsuperstruct_init.o lsuperstruct_interface.o lsuperstruct_timeloop.o lsuperstruct_parallel.o lsuperstruct_dump.o lextobj.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct.F90
lsuperstruct_IBM.o: ./lsuperstruct_IBM.F90 lbe.h lbe_force_interface.o lbe_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o lsuperstruct_data.o lextobj.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_IBM.F90
lsuperstruct_data.o: ./lsuperstruct_data.F90 lbe.h lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_data.F90
lsuperstruct_dump.o: ./lsuperstruct_dump.F90 lbe.h lbe_globals.o lbe_helper.o lbe_io_helper.o lbe_io.o lbe_io_xdrf.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o lsuperstruct_data.o lsuperstruct_helper.o lmesh.o lextobj.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_dump.F90
lsuperstruct_helper.o: ./lsuperstruct_helper.F90 lbe.h lbe_bdist.o lbe_globals.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o lsuperstruct_data.o lextobj.o lmesh.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_helper.F90
lsuperstruct_init.o: ./lsuperstruct_init.F90 lbe.h lbe_globals.o lbe_helper.o lbe_init_rock.o lbe_io_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o lsuperstruct_data.o lsuperstruct_dump.o lsuperstruct_IBM.o lsuperstruct_interface.o lsuperstruct_parallel.o lsuperstruct_helper.o lsuperstruct_timeloop.o lextobj.o lmesh.o lbe_io_xdrf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_init.F90
lsuperstruct_interface.o: ./lsuperstruct_interface.F90 lbe.h lbe_bdist.o lbe_force_interface.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_globals.o lbe_types.o lsuperstruct_data.o lsuperstruct_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_interface.F90
lsuperstruct_parallel.o: ./lsuperstruct_parallel.F90 lbe.h lbe_globals.o lbe_log.o lbe_parallel.o lbe_parms.o lmesh.o lsuperstruct_data.o lsuperstruct_helper.o lextobj.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_parallel.F90
lsuperstruct_timeloop.o: ./lsuperstruct_timeloop.F90 lbe.h hoshen_kopelman.o lbe_force_interface.o lbe_helper.o lbe_log.o lbe_parallel.o lbe_parms.o lbe_types.o lsuperstruct_data.o lsuperstruct_helper.o lextobj.o lsuperstruct_interface.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./lsuperstruct_timeloop.F90
luxury_rng.o: ./luxury_rng.F90 lbe.h lbe_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./luxury_rng.F90
map.o: ./map.F90 lbe.h lbe_globals.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I.	./map.F90
mesh_3d.o: ./mesh_3d.F90 lbe_log.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./mesh_3d.F90
pppm.o: ./pppm.F90 pppm_global.o lbe_parallel.o lbe_log.o lbe_elec_timer.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./pppm.F90
pppm_coeff.o: ./pppm_coeff.F90 pppm_global.o lbe_log.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./pppm_coeff.F90
pppm_global.o: ./pppm_global.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./pppm_global.F90
pppm_remap.o: ./pppm_remap.F90 pppm_global.o lbe_parallel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./pppm_remap.F90
SRC = ./lbe_elec_interface.F90 ./lbe_tracer_helper.F90 ./lbe_md_potential_gb.F90 ./lbe_md_globals.F90 ./lbe_elec_timer.F90 ./lbe_tracer_md_fluid_ladd.F90 ./lbe_elec_poisson_solver.F90 ./lbe_md_fluid_friction.F90 ./lbe_md_bc_leesedwards.F90 ./hoshen_kopelman.F90 ./lbe_tracer_memory.F90 ./lbe.F90 ./lbe_md_potential_bprhag.F90 ./lbe_axisym.F90 ./lbe_md_parallel.F90 ./lbe_analysis.F90 ./lbe_io_stress.F90 ./lbe_elec_input.F90 ./lbe_elec_init.F90 ./lbe_md_dynamic.F90 ./lsuperstruct_helper.F90 ./lsuperstruct_data.F90 ./lbe_parms.F90 ./lbe_io_interface.F90 ./lbe_tracer_input.F90 ./lbe_force.F90 ./lbe_md_input.F90 ./lbe_md_rock.F90 ./lbe_md_potential_spring.F90 ./lbe_io_xdrf.F90 ./lbe_io_helper.F90 ./lbe_md_potential.F90 ./pppm_coeff.F90 ./lbe_init_rock.F90 ./lbe_mean_square_displacement.F90 ./lbe_io_arrstats.F90 ./lbe_elec_globals.F90 ./lbe_collision.F90 ./lbe_tracer.F90 ./lbe_parallel.F90 ./lbe_timer.F90 ./lbe_collision_simple.F90 ./lbe_md_fluid_ladd_parms.F90 ./lbe_md_rock_bprh.F90 ./lbe_leesedwards.F90 ./pppm_remap.F90 ./lbe_md_debug.F90 ./lbe_md_potential_hertz.F90 ./lbe_MRT.F90 ./lbe_elec_timestep.F90 ./lbe_md_fluid_ladd.F90 ./lbe_md_init.F90 ./lbe_mass_scaler.F90 ./lbe_md_fluid_ladd_mc.F90 ./fft_3d.c ./lextobj.F90 ./lbe_tracer_init.F90 ./lbe_elec_output.F90 ./pppm_global.F90 ./lbe_elec_fluxes.F90 ./lbe_io.F90 ./lsuperstruct_init.F90 ./lbe_bc.F90 ./lbe_md_rock_lj.F90 ./lbe_elec_helper.F90 ./lbe_md_boundary_condition.F90 ./mesh_3d.F90 ./luxury_rng.F90 ./lbe_helper.F90 ./lsuperstruct_interface.F90 ./lbe_md_potential_dlvo.F90 ./lbe_tracer_recoloring.F90 ./lbe_md_fluid_ladd_dlub.F90 ./lsuperstruct_dump.F90 ./lbe_md_analysis.F90 ./lbe_md_fluid_ladd_club.F90 ./lbe_elec_parallel.F90 ./lbe_md_rock_friction.F90 ./lbe_md_scale_vel.F90 ./lbe_md_memory.F90 ./lsuperstruct_IBM.F90 ./lbe_log.F90 ./lbe_io_hdf5.F90 ./map.F90 ./lbe_types.F90 ./lbe_md_magnetic.F90 ./lbe_bdist.F90 ./lbe_md_rock_dlvo.F90 ./lsuperstruct_parallel.F90 ./lbe_tracer_interface.F90 ./lbe_tracer_output.F90 ./lbe_md_potential_lubrication.F90 ./lsuperstruct_timeloop.F90 ./lbe_md_growing_stage.F90 ./lbe_md.F90 ./lbe_elec_forces.F90 ./lbe_invasion.F90 ./lbe_force_interface.F90 ./lbe_tracer_parallel.F90 ./lbe_md_potential_lj.F90 ./lbe_md_potential_bprh.F90 ./lbe_galilean_stabilizer.F90 ./lbe_init.F90 ./lbe_md_output.F90 ./lbe_md_helper.F90 ./lbe_md_fluid.F90 ./lbe_force_kolmogorov.F90 ./lmesh.F90 ./lbe_tracer_globals.F90 ./lbe_globals.F90 ./lbe_md_lyapunov.F90 ./lbe_force_constant.F90 ./lbe_md_fluid_tracer.F90 ./lbe_md_rand_walk.F90 ./lsuperstruct.F90 ./lbe_md_interface.F90 ./pppm.F90 ./lbe_io_checkpoint.F90 lbe_version.h lbe.h
OBJ = lbe_elec_interface.o lbe_tracer_helper.o lbe_md_potential_gb.o lbe_md_globals.o lbe_elec_timer.o lbe_tracer_md_fluid_ladd.o lbe_elec_poisson_solver.o lbe_md_fluid_friction.o lbe_md_bc_leesedwards.o hoshen_kopelman.o lbe_tracer_memory.o lbe.o lbe_md_potential_bprhag.o lbe_axisym.o lbe_md_parallel.o lbe_analysis.o lbe_io_stress.o lbe_elec_input.o lbe_elec_init.o lbe_md_dynamic.o lsuperstruct_helper.o lsuperstruct_data.o lbe_parms.o lbe_io_interface.o lbe_tracer_input.o lbe_force.o lbe_md_input.o lbe_md_rock.o lbe_md_potential_spring.o lbe_io_xdrf.o lbe_io_helper.o lbe_md_potential.o pppm_coeff.o lbe_init_rock.o lbe_mean_square_displacement.o lbe_io_arrstats.o lbe_elec_globals.o lbe_collision.o lbe_tracer.o lbe_parallel.o lbe_timer.o lbe_collision_simple.o lbe_md_fluid_ladd_parms.o lbe_md_rock_bprh.o lbe_leesedwards.o pppm_remap.o lbe_md_debug.o lbe_md_potential_hertz.o lbe_MRT.o lbe_elec_timestep.o lbe_md_fluid_ladd.o lbe_md_init.o lbe_mass_scaler.o lbe_md_fluid_ladd_mc.o fft_3d.o lextobj.o lbe_tracer_init.o lbe_elec_output.o pppm_global.o lbe_elec_fluxes.o lbe_io.o lsuperstruct_init.o lbe_bc.o lbe_md_rock_lj.o lbe_elec_helper.o lbe_md_boundary_condition.o mesh_3d.o luxury_rng.o lbe_helper.o lsuperstruct_interface.o lbe_md_potential_dlvo.o lbe_tracer_recoloring.o lbe_md_fluid_ladd_dlub.o lsuperstruct_dump.o lbe_md_analysis.o lbe_md_fluid_ladd_club.o lbe_elec_parallel.o lbe_md_rock_friction.o lbe_md_scale_vel.o lbe_md_memory.o lsuperstruct_IBM.o lbe_log.o lbe_io_hdf5.o map.o lbe_types.o lbe_md_magnetic.o lbe_bdist.o lbe_md_rock_dlvo.o lsuperstruct_parallel.o lbe_tracer_interface.o lbe_tracer_output.o lbe_md_potential_lubrication.o lsuperstruct_timeloop.o lbe_md_growing_stage.o lbe_md.o lbe_elec_forces.o lbe_invasion.o lbe_force_interface.o lbe_tracer_parallel.o lbe_md_potential_lj.o lbe_md_potential_bprh.o lbe_galilean_stabilizer.o lbe_init.o lbe_md_output.o lbe_md_helper.o lbe_md_fluid.o lbe_force_kolmogorov.o lmesh.o lbe_tracer_globals.o lbe_globals.o lbe_md_lyapunov.o lbe_force_constant.o lbe_md_fluid_tracer.o lbe_md_rand_walk.o lsuperstruct.o lbe_md_interface.o pppm.o lbe_io_checkpoint.o
clean: neat
	-rm -f .lbe.cppdefs $(OBJ) lbe
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lbe: $(OBJ) 
	$(LD) $(OBJ) -o lbe  $(LDFLAGS)
VTMP="lbe_version.tmp"
VH="lbe_version.h"

version:
	@./../get-version.sh -h > $(VTMP)
	@echo "#define LBE_FLAGS \"$(DEFFLAGS) $(MAKEFFLAGS)\"" >> $(VTMP)
	@echo "#define LBE_PLATFORM \"$(PLATFORM)\"" >> $(VTMP)
	@diff $(VTMP) $(VH) > /dev/null 2>&1 ; if [ "$$?" -eq 0 ] ; then rm -rf $(VTMP) ; else echo "Creating new $(VH)"; mv $(VTMP) $(VH); fi
