#include "lbe.h"

!> delegates to magnetic force
module lbe_md_magnetic_module
#ifdef MD
!#ifdef MD_MAG

    use lbe_globals_module, only: pi,tsize, input_dfile_unit,myrankc
    use lbe_helper_module, only: cross_product
    use lbe_md_boundary_condition_module, only: rdst
    use lbe_md_fluid_ladd_module, only: particle_lubrication&
        &,rc_lubrication_particle
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_md_globals_module
    use lbe_log_module
    use lbe_md_helper_module, only: count_particles_all,error_md,log_msg_md&
         &,number_density,log_msg_md_hdr
    use lbe_md_output_module, only: dump_potentials,n_dump
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: nt, inp_file,arg_input_dfile_set,arg_input_dfile
    implicit none
    private

    include 'mpif.h'
    public  input_magnetic, force_and_torque_magnetic,force_and_torque_magnetic_inter
    !> switches magnetic force on/off
    logical,save,public :: md_magnetic=.false.
    logical,save,public :: md_mag_o=.false. !< magnetic magnetic dipole direction same as particle orientation (false) or orthogonal

    real(kind=rk), save,public :: cons_mag_dip = 0.0_rk !< magnitude of dispole moment
    real(kind=rk), save,public :: cons_mag_fie(3) = (/0.0_rk,0.0_rk,0.0_rk/) !< external constant magnetic field
    real(kind=rk), save,public :: gra_mag_fie(3) = (/0.0_rk,0.0_rk,0.0_rk/) !< external linear magnetic field
    real(kind=rk), parameter :: mu0 = 0.0000001_rk
    integer,save,public :: T_switch = 0 !< switch on oscillating magnetic field after this time step
    real(kind=rk), save, public :: T_fre = 0.0_rk  !< Frequency of ocillating mangetic field
    real(kind=rk), save, public :: T_con(3) = (/0.0_rk,0.0_rk,0.0_rk/) !< constant coefficient of ocillating magnetic field         
    
    namelist /md_magnetic_para/ cons_mag_dip,cons_mag_fie,gra_mag_fie, T_switch,T_fre, T_con, md_mag_o
contains

!>   read magnetic section from md input file
  subroutine input_magnetic
        integer ierror
        call log_msg_md_hdr("Reading MD Magnetic Input")
        use_rotation = .true.                                                                                                                                                                
        calculate_orientations = .true. 
        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_magnetic_para,err=100)
           close (unit=md_input_file_unit,err=100)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
         read(UNIT = input_dfile_unit, NML = md_magnetic_para, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          call log_ws()
        end if

        write(msgstr,"('cons_mag_dip            = ',F24.10)") cons_mag_dip
        call log_msg(msgstr) 

        write(msgstr,"('cons_mag_fie            = (',F16.10,',',F16.10,',',F16.10,')')") cons_mag_fie(1), cons_mag_fie(2), cons_mag_fie(3)
        call log_msg(msgstr) 
 
        write(msgstr,"('gra_mag_fie            = (',F16.10,',',F16.10,',',F16.10,')')") gra_mag_fie(1), gra_mag_fie(2), gra_mag_fie(3)
        call log_msg(msgstr)
        write(msgstr,"('T_switch        = ',I0)") T_switch
         call log_msg(msgstr)

        write(msgstr,"('T_fre            = ',F16.10)") T_fre
        call log_msg(msgstr)
        
        write(msgstr,"('T_con            = (',F16.10,',',F16.10,',',F16.10,')')") T_con(1), T_con(2), T_con(3)
        call log_msg(msgstr)
       
		  write(msgstr,"('md_mag_o)      = ',L1)") md_mag_o
          call log_msg(msgstr)

        call log_ws()
       
        call MPI_Bcast(cons_mag_dip,1,MPI_REAL8,0,comm_cart,ierror)      
        call MPI_Bcast(cons_mag_fie,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(gra_mag_fie,3,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(T_switch,1,MPI_INTEGER,0,Comm_cart,ierror)
        call MPI_Bcast(T_fre,1,MPI_REAL8,0,Comm_cart,ierror)
        call MPI_Bcast(T_con,3,MPI_REAL8,0,Comm_cart,ierror)
        call MPI_Bcast(md_mag_o,1,MPI_LOGICAL,0,comm_cart,ierror)
        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_magnetic


    subroutine force_and_torque_magnetic
        integer i,ii
        real(kind=rk) :: fm(3),tm(3),po(3)
        real(kind=rk),dimension(3) :: mag_fie
        real(kind=rk) :: gra_x(3)=(/0.0_rk,0.0_rk,0.0_rk/)
        real(kind=rk) :: gra_y(3)=(/0.0_rk,0.0_rk,0.0_rk/)
        real(kind=rk) :: gra_z(3)=(/0.0_rk,0.0_rk,0.0_rk/)
        
        calculate_orientations = .true.   
        gra_x(1)=gra_mag_fie(1)
        gra_y(2)=gra_mag_fie(2)
        gra_z(3)=gra_mag_fie(3)
        i = atompnt
        do ii = 1,nlocal
 			if (nt >= T_switch) then
                    mag_fie(1)=gra_mag_fie(1)*P(i)%x(1)+cons_mag_fie(1)+T_con(1)*sin(nt*2*pi*T_fre)
                    mag_fie(2)=gra_mag_fie(2)*P(i)%x(2)+cons_mag_fie(2)+T_con(2)*sin(nt*2*pi*T_fre+pi/2)
                    mag_fie(3)=gra_mag_fie(3)*P(i)%x(3)+cons_mag_fie(3)+T_con(3)*sin(nt*2*pi*T_fre+pi)
  				else 
                      mag_fie(1)=gra_mag_fie(1)*P(i)%x(1)+cons_mag_fie(1)
                      mag_fie(2)=gra_mag_fie(2)*P(i)%x(2)+cons_mag_fie(2)
                      mag_fie(3)=gra_mag_fie(3)*P(i)%x(3)+cons_mag_fie(3)
 			endif

 		 if(md_mag_o) then
							po = mag_o(P(i))
			else 
							po = P(i)%o 
		endif
      !print *, "orientation is : %f\n", P(i)%o
      !print *, "magnetic dipole orientation is : %f\n", po
       if (magdispersity) then 

                      fm(1) = P(i)%mag(1)*dot_product(po,gra_x)
                      fm(2) = P(i)%mag(1)*dot_product(po,gra_y)
                      fm(3) = P(i)%mag(1)*dot_product(po,gra_z)
                      ! torque on particle i
                      tm(:) = P(i)%mag(1)*cross_product(po,mag_fie)
               else 
                    ! force on particle i 
                    fm(1) = cons_mag_dip*dot_product(po,gra_x)
                    fm(2) = cons_mag_dip*dot_product(po,gra_y)
                    fm(3) = cons_mag_dip*dot_product(po,gra_z) 
                    ! torque on particle i

                    tm(:) = cons_mag_dip*cross_product(po,mag_fie)
        end if
                    P(i)%fmi = P(i)%fmi +fm
                    P(i)%tmi = P(i)%tmi +tm                    
                    P(i)%f = P(i)%f + fm
                    P(i)%t = P(i)%t + tm
           i = list(i)
        enddo
      end subroutine force_and_torque_magnetic
      subroutine force_and_torque_magnetic_inter
        integer i,j,k,ii
        real(kind=rk) :: arij,bracket(3),e,e_pot_2,eomr,f(3),fracm,fracm_arij&
             &,fracp,fracp_arij,grad_ri_r(3),omr,omxoioj,opxoioj,r,r_X2,r2&
             &,rij(3),rji(3),rij2,rijoi,rijoimrijoj,rijoiprijoj,rijoj,s2,sqrt1,sqrt2&
             &,ti(3),tj(3),urij(3),urji(3),xoioj
        real(kind=rk) :: fbi(3),tbi(3),fbj(3),tbj(3), mimj
        i = atompnt
        do ii = 1,nlocal

           do k = nnlist(ii),nnlist(ii+1)-1

              j = nlist(k)
              rij = rdst(P(j),P(i))
              rji = rdst(P(i),P(j))
              ! (absolute value of rij)**2
              rij2 = dot_product(rij,rij)
              sqrt1= sqrt(rij2)
              urij(:) = rij(:)/sqrt1
              urji(:) = rji(:)/sqrt1

             if (magdispersity) then 
                   mimj = P(i)%mag(1)*P(j)%mag(1)
					else
					mimj =cons_mag_dip**2
					end if

               !force of particle j on particle i
              fbi= 3*mu0*mimj/(rij2**2)&
                 &*(dot_product(P(i)%o,urji)*P(j)%o(:)&
                 &+dot_product(P(j)%o,urji)*P(i)%o(:)&
                 &-(5*dot_product(P(j)%o,urji)*dot_product(P(i)%o,urji)-dot_product(P(j)%o,P(i)%o))*urji(:))
               !torque of particle j on particle i
              tbi= mu0*mimj/(rij2*sqrt1)&
                 &*(3*dot_product(P(j)%o,urji)*cross_product(P(i)%o,urji)-cross_product(P(i)%o,P(j)%o))
       !print *, "magnetic dipole is : %f\n", P(i)%mag

        P(i)%fmi = P(i)%fmi + fbi
        P(i)%tmi = P(i)%tmi + tbi 
       P(i)%f = P(i)%f + fbi 
       P(i)%t = P(i)%t + tbi 

       if (j.le.npmax) then
          !force of particle i on particle j is opposite w.r.t the force of particle j on particle i
          P(j)%fmi = P(j)%fmi -fbi
          P(j)%f = P(j)%f - fbi
          !torque of particle i on particle j
                  tbj= mu0*mimj/(rij2*sqrt1)&
                          &*(3*dot_product(P(i)%o,urij)*cross_product(P(j)%o,urij)-cross_product(P(j)%o,P(i)%o))
          P(j)%tmi = P(j)%tmi+tbj
          P(j)%t = P(j)%t + tbj
   end if  
! It is necessary to check lubrication using rij or rji
!if (lubrication) call particle_lubrication(rij2,rij,i,j)
enddo
i = list(i)
enddo
end subroutine  force_and_torque_magnetic_inter

  !> returns the orientation of magnetic dipole of particle \i
      pure function mag_o(p)
          type(md_particle_type),intent(in) :: p
          real(kind=rk) :: mag_o(3),theta
   theta = p%mag(2)      
   mag_o(1) = theta/abs(theta)*sqrt(p%o(3)**2/((p%o(1)+p%o(2)*tan(theta))**2+p%o(3)*p%o(3)*(1+tan(theta)**2)))
   mag_o(2) = mag_o(1)*tan(theta) 
   mag_o(3) = theta/abs(theta)*sqrt(1- mag_o(1)**2- mag_o(2)**2 )        
            
   end function mag_o


#endif
end module lbe_md_magnetic_module

