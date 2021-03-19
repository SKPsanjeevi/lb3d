#ifdef P3M

! LAMMPS 2001 - Molecular Dynamics Simulator
! Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
! Steve Plimpton, sjplimp@sandia.gov
!
! Copyright (1998-2001) Sandia Corporation.  Under the terms of Contract
! DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
! certain rights in this software.  This software is distributed under 
! the GNU General Public License.
!
! See the README file in the top-level LAMMPS directory.

! Contributing authors: Roy Pollock (LLNL), Paul Crozier (Sandia)

! -------------------------------------------------------------------------
! setup of 2 decompositions (3-d brick and 2-d FFT) and various coeffs
!   for PPPM calculation
! iflag = 0 section is only computed once since those quantities
!   depend only on grid size, it is called at beginning of a run
! NOTE: if box volume changes enough that a different resolution grid
!   is needed then need to put in calls to recompute iflag = 0 section
! iflag = 1 section must be recomputed whenever simulation volume changes,
!   it is called during a run when a volume change occurs

      subroutine pppm_coeff(iflag)
      use lammps_global_module
      use mpi
      use lbe_log_module
      use lbe_parallel_module, only: comm_cart

      implicit none

! argument variables

      integer iflag

! local variables

      integer istatus(mpi_status_size)
      integer nlo,nhi,itotal,ntotal,ierror,nplanes,jflag
      integer kflag,nx,ny,nz,idelx,idely,idelz,nxx,nyy
      integer npey_fft,npez_fft,me_y,me_z,itmp,i,j,k,kper,l
      integer lper,m,mper,nzz,nfft_brick,icount
      real*8 delxinv,delyinv,delzinv,pshift,cuthalf
      real*8 unitk,sqk,vterm

! adjustment of z dimension for 2d slab PPPM
! 3d PPPM just uses zprd since volfactor = 1.0

      zprd_slab = zprd*slab_volfactor

      if (iflag == 0) then

! PPPM grid size and cutoff parameter from estimation routine

        call get_p3m_parms

! nlo_in,nhi_in = lower/upper limits of the 3-d sub-brick of
!   global PPPM grid that I own in each dimension WITHOUT ghost cells
! global indices range from 0 to N-1
! for slab PPPM, assign z grid as if it were not extended

        nxlo_in = me(1)*nx_pppm/pgrid(1)
        nxhi_in = (me(1)+1)*nx_pppm/pgrid(1) - 1
        nylo_in = me(2)*ny_pppm/pgrid(2)
        nyhi_in = (me(2)+1)*ny_pppm/pgrid(2) - 1
        nzlo_in = me(3)*nint(nz_pppm/slab_volfactor)/pgrid(3)
        nzhi_in = (me(3)+1)*nint(nz_pppm/slab_volfactor)/pgrid(3) - 1

! nlower,nupper = stencil size for mapping particles to PPPM grid
        nlower = -(orderflag-1)/2
        nupper = orderflag/2

! nlo_out,nhi_out = lower/upper limits of the 3-d sub-brick of
!   global PPPM grid that I own in each dimension WITH ghost cells
! ghost 3-d brick boundaries = inner + stencil + particles moving skin/2.0
! nlo/nhi = global coords of grid pt to "lower left" of smallest/largest
!           position a particle in my box can be at
! add/subtract 4096 to avoid problem of int(-0.75) = 0 when it needs to be -1

        delxinv = nx_pppm/xprd
        delyinv = ny_pppm/yprd
        delzinv = nz_pppm/zprd_slab

        if (mod(orderflag,2).ne.0) then
          pshift = 4096.5
        else
          pshift = 4096.0
        endif

        cuthalf = skin/2.0

        nlo = int((border(1,1)-cuthalf-box(1,1))*delxinv+pshift) - 4096
        nhi = int((border(2,1)+cuthalf-box(1,1))*delxinv+pshift) - 4096
        nxlo_out = nlo + nlower
        nxhi_out = nhi + nupper

        nlo = int((border(1,2)-cuthalf-box(1,2))*delyinv+pshift) - 4096
        nhi = int((border(2,2)+cuthalf-box(1,2))*delyinv+pshift) - 4096
        nylo_out = nlo + nlower
        nyhi_out = nhi + nupper

        nlo = int((border(1,3)-cuthalf-box(1,3))*delzinv+pshift) - 4096
        nhi = int((border(2,3)+cuthalf-box(1,3))*delzinv+pshift) - 4096
        nzlo_out = nlo + nlower
        nzhi_out = nhi + nupper

! for slab-geometry-corrected calculations, change the grid boundary for
!   processors at +z end to include the empty volume between periodically
!   repeating slabs
! for slab PPPM, want charge data communicated from -z proc to +z proc,
!   but not vice versa, also want field data communicated from +z proc to
!   -z proc, but not vice versa
! this is accomplished by nzhi_in = nzhi_out on +z end (no ghost cells)

        if (slabflag == 1 .and. me(3)+1 == pgrid(3)) then
          nzhi_in =  nz_pppm - 1
          nzhi_out = nz_pppm - 1
        endif

       if(.not.ghost_cells)then
!gjp code to remove ghost cells, required when using poisson solver alone
         nxlo_out = nxlo_in
         nxhi_out = nxhi_in
         nylo_out = nylo_in
         nyhi_out = nyhi_in
         nzlo_out = nzlo_in
         nzhi_out = nzhi_in
       endif

! exchange messages to determine how many planes I will send,recv in each dir
! if no neighbor proc exists, values comes from self
!   since I have ghosts regardless

        nplanes = nxlo_in - nxlo_out
        if (mpart(1,1).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,1),0, &
     &         comm_cart,ierror)
          call mpi_recv(nxhi_ghost,1,mpi_integer,mpart(2,1),0, &
     &         comm_cart,istatus,ierror)
        else
          nxhi_ghost = nplanes
        endif

        nplanes = nxhi_out - nxhi_in
        if (mpart(2,1).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,1),0, &
     &         comm_cart,ierror)
          call mpi_recv(nxlo_ghost,1,mpi_integer,mpart(1,1),0, &
     &         comm_cart,istatus,ierror)
        else
          nxlo_ghost = nplanes
        endif

        nplanes = nylo_in - nylo_out
        if (mpart(1,2).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,2),0, &
     &         comm_cart,ierror)
          call mpi_recv(nyhi_ghost,1,mpi_integer,mpart(2,2),0, &
     &         comm_cart,istatus,ierror)
        else
          nyhi_ghost = nplanes
        endif

        nplanes = nyhi_out - nyhi_in
        if (mpart(2,2).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,2),0, &
     &         comm_cart,ierror)
          call mpi_recv(nylo_ghost,1,mpi_integer,mpart(1,2),0, &
     &         comm_cart,istatus,ierror)
        else
          nylo_ghost = nplanes
        endif

        nplanes = nzlo_in - nzlo_out
        if (mpart(1,3).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,3),0, &
     &         comm_cart,ierror)
          call mpi_recv(nzhi_ghost,1,mpi_integer,mpart(2,3),0, &
     &         comm_cart,istatus,ierror)
        else
          nzhi_ghost = nplanes
        endif

        nplanes = nzhi_out - nzhi_in
        if (mpart(2,3).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,3),0, &
     &         comm_cart,ierror)
          call mpi_recv(nzlo_ghost,1,mpi_integer,mpart(1,3),0, &
     &         comm_cart,istatus,ierror)
        else
          nzlo_ghost = nplanes
        endif

! test that ghost overlap is not bigger than my sub-domain

        jflag = 0
        if (nxlo_ghost.gt.nxhi_in-nxlo_in+1) jflag = 1
        if (nxhi_ghost.gt.nxhi_in-nxlo_in+1) jflag = 1
        if (nylo_ghost.gt.nyhi_in-nylo_in+1) jflag = 1
        if (nyhi_ghost.gt.nyhi_in-nylo_in+1) jflag = 1
        if (nzlo_ghost.gt.nzhi_in-nzlo_in+1) jflag = 1
        if (nzhi_ghost.gt.nzhi_in-nzlo_in+1) jflag = 1

        call mpi_allreduce(jflag,kflag,1,mpi_integer,mpi_sum, &
     &       comm_cart,ierror)

        if (kflag.gt.0) &
     &       call error('PPPM stencil extends beyond neighbor proc'// &
     &       ' - reduce PPPM order')

! decomposition of FFT mesh
! proc owns entire x-dimension, clump of columns in y,z dimensions
! npey_fft,npez_fft = # of procs in y,z dims
! if ncores is small enough, proc can own 1 or more entire xy planes,
!   else proc owns 2-d sub-blocks of yz plane
! me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
! nlo_fft,nhi_fft = lower/upper limit of the section
!   of the global FFT mesh that I own in each dimension
! indices range from 0 to N-1

        if (nz_pppm >= ncores) then
          npey_fft = 1
          npez_fft = ncores
        else
          call proc2grid2d(ncores,ny_pppm,nz_pppm,npey_fft,npez_fft)
        endif

        me_y = mod(node,npey_fft)
        me_z = node/npey_fft

        nxlo_fft = 0
        nxhi_fft = nx_pppm - 1
        nylo_fft = me_y*ny_pppm/npey_fft
        nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1
        nzlo_fft = me_z*nz_pppm/npez_fft
        nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1

! allocate all memory needed for PPPM solution

! allocate grid arrays for PPPM mesh on this proc, including ghosts

        maxgrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) * &
     &       (nzhi_out-nzlo_out+1)

        call mpi_allreduce(maxgrid,ntotal,1,mpi_integer, &
     &       mpi_max,comm_cart,ierror)

        write (msgstr,*) 'PPPM max grid size =',ntotal
        call log_msg(msgstr)

        allocate (density_brick(maxgrid))
!gjp phi_brick is new to LAMMPS code
        allocate (phi_brick(maxgrid))
        allocate (vdx_brick(maxgrid))
        allocate (vdy_brick(maxgrid))
        allocate (vdz_brick(maxgrid))

! allocate FFT arrays for PPPM mesh on this proc, without ghosts
! nfft = FFT points in FFT decomposition on this proc
! nfft_brick = FFT points in 3-d brick-decomposition on this proc

        nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) * &
     &       (nzhi_fft-nzlo_fft+1)
        nfft_brick = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) * &
     &       (nzhi_in-nzlo_in+1)
        maxfft = max(nfft,nfft_brick)

        call mpi_allreduce(maxfft,ntotal,1,mpi_integer, &
     &       mpi_max,comm_cart,ierror)

        write (msgstr,*) 'PPPM max FFT size =',ntotal
        call log_msg(msgstr)

        allocate (density_fft(maxfft))
        allocate (greensfn(maxfft))
        allocate (workvec1(maxfft))
        allocate (workvec2(maxfft))
!gjp workvec_E is required when using the FFT E solver
        allocate (workvec_E(maxfft))

        allocate (vg(6,maxfft))

! allocate per-particle arrays

        allocate (partgrid(3,maxown))
        allocate (ek(3,maxown))

! allocate space for fkvecs_xyz factors

        allocate (fkvecs_x(nxlo_fft:nxhi_fft))
        allocate (fkvecs_y(nylo_fft:nyhi_fft))
        allocate (fkvecs_z(nzlo_fft:nzhi_fft))

! allocate buffer space for use in brick2fft and fillbrick
! idel = max # of ghost planes to send or recv in +/- direction of each dim
! nxx,nyy,nzz = max # of cells to send in each dim
! maxpbuf = max in any dim, factor of 3 for components of vd_xyz in fillbrick

        nx = nxhi_out - nxlo_out + 1
        ny = nyhi_out - nylo_out + 1
        nz = nzhi_out - nzlo_out + 1

        idelx = max(nxlo_ghost,nxhi_ghost)
        idelx = max(idelx,nxhi_out-nxhi_in)
        idelx = max(idelx,nxlo_in-nxlo_out)

        idely = max(nylo_ghost,nyhi_ghost)
        idely = max(idely,nyhi_out-nyhi_in)
        idely = max(idely,nylo_in-nylo_out)

        idelz = max(nzlo_ghost,nzhi_ghost)
        idelz = max(idelz,nzhi_out-nzhi_in)
        idelz = max(idelz,nzlo_in-nzlo_out)

        nxx = idelx * ny * nz
        nyy = idely * nx * nz
        nzz = idelz * nx * ny

        itotal = max(nxx,nyy)
        itotal = max(itotal,nzz)
        maxpbuf = 3*itotal

        call mpi_allreduce(maxpbuf,ntotal,1,mpi_integer, &
     &       mpi_max,comm_cart,ierror)
        write (msgstr,*) 'PPPM max buffer size =',ntotal
        call log_msg(msgstr)

        allocate (pbuf1(maxpbuf),pbuf2(maxpbuf))

! create 2 FFT plans and remap plan
! 1st FFT plan keeps data in FFT decompostion
! 2nd FFT plan returns data in 3-d brick decomposition
! remap plan takes data from 3-d brick to FFT decomposition
! all indices must be converted to range from 1 to N

        call fft_3d_create_plan(comm_cart, &
     &       nx_pppm,ny_pppm,nz_pppm, &
     &       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1, &
     &       nzlo_fft+1,nzhi_fft+1, &
     &       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1, &
     &       nzlo_fft+1,nzhi_fft+1, &
     &       0,0,itmp,plan1_fft)

        call fft_3d_create_plan(comm_cart, &
     &       nx_pppm,ny_pppm,nz_pppm, &
     &       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1, &
     &       nzlo_fft+1,nzhi_fft+1, &
     &       nxlo_in+1,nxhi_in+1,nylo_in+1,nyhi_in+1, &
     &       nzlo_in+1,nzhi_in+1, &
     &       0,0,itmp,plan2_fft)

        call remap_3d_create_plan(comm_cart, &
     &       nxlo_in+1,nxhi_in+1,nylo_in+1,nyhi_in+1, &
     &       nzlo_in+1,nzhi_in+1, &
     &       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1, &
     &       nzlo_fft+1,nzhi_fft+1, &
     &       1,0,0,2,plan_remap)

      endif

! computations that must be re-done whenever simulation volume changes

! pre-compute fkvecs_xyz(*) for my FFT grid pts

      unitk = (2.0*3.141592654/xprd)
      do k = nxlo_fft,nxhi_fft
        kper = k - nx_pppm*int(2*k/nx_pppm)
        fkvecs_x(k) = unitk*kper
      enddo

      unitk = (2.0*3.141592654/yprd)
      do l = nylo_fft,nyhi_fft
        lper = l - ny_pppm*int(2*l/ny_pppm)
        fkvecs_y(l) = unitk*lper
      enddo

      unitk = (2.0*3.141592654/zprd_slab)
      do m = nzlo_fft,nzhi_fft
        mper = m - nz_pppm*int(2*m/nz_pppm)
        fkvecs_z(m) = unitk*mper
      enddo

! pre-compute virial coefficients

      icount = 0
      do k = nzlo_fft,nzhi_fft
        do j = nylo_fft,nyhi_fft
          do i = nxlo_fft,nxhi_fft
            icount = icount + 1
            sqk = fkvecs_x(i)*fkvecs_x(i) + &
     &            fkvecs_y(j)*fkvecs_y(j) + &
     &            fkvecs_z(k)*fkvecs_z(k)
            if (sqk == 0.0) then
              vg(1,icount) = 0.0
              vg(2,icount) = 0.0
              vg(3,icount) = 0.0
              vg(4,icount) = 0.0
              vg(5,icount) = 0.0
              vg(6,icount) = 0.0
            else
              vterm = -2.0*(1.0/sqk + 0.25/(gewald*gewald))
              vg(1,icount) = 1.0 + vterm*fkvecs_x(i)*fkvecs_x(i)
              vg(2,icount) = 1.0 + vterm*fkvecs_y(j)*fkvecs_y(j)
              vg(3,icount) = 1.0 + vterm*fkvecs_z(k)*fkvecs_z(k)
              vg(4,icount) = vterm*fkvecs_x(i)*fkvecs_y(j)
              vg(5,icount) = vterm*fkvecs_x(i)*fkvecs_z(k)
              vg(6,icount) = vterm*fkvecs_y(j)*fkvecs_z(k)
            endif
          enddo
        enddo
      enddo

! pre-compute Green's function coeffs

      if(jflag==0)then
        call make_nonadapted_greensfn(orderflag,gewald, &
     &     nx_pppm,ny_pppm,nz_pppm,xprd,yprd,zprd_slab, &
     &     nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft, &
     &     greensfn)
      else
        call make_greensfn(orderflag,gewald, &
     &     nx_pppm,ny_pppm,nz_pppm,xprd,yprd,zprd_slab, &
     &     nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft, &
     &     greensfn)

      endif
      return
      end

! -------------------------------------------------------------------------
! get p3m parameters
!   INPUT:
!           orderflag------------------p3m assignment scheme order
!           long_prec------------------desired relative error in forces
!           cutcoul--------------------real space cutoff
!           xprd,yprd,zprd_slab--------edge lengths of periodic cell
!           natoms---------------------total number of atoms
!           qsqsum---------------------sum of the squared charges
!   OUTPUT:
!           gewald---------------------Ewald parameter
!           nx_pppm,ny_pppm,nz_pppm----P3m grid sizes
!

      subroutine get_p3m_parms
      use lammps_global_module
!gjp including this for subroutine error
      use lbe_log_module

      implicit none

! local variables (old method)

      integer maxorder
      real*8 a,b,gdel
      parameter (maxorder=16)
      dimension a(maxorder),b(maxorder)
      data a /5.0,5.78,6.93,7.76,8.55,9.10,9.66,10.17,10.60, &
     &     11.1,14,15,16,17,18,19/
      data b /2.0,2.25,3.41,4.62,5.85,7.05,8.20,9.28,10.34, &
     &     11.3,12,13,14,15,16,17/

! local variables (new method)

      logical factorable,usenew
      integer nfacs,factors(3),ncount,m
      real*8 small,large,acons(7,0:6),h,h1,h2,er,er1,er2,gew1,gew2
      real*8 hx,hy,hz,summation,spr,lpr,lprx,lpry,lprz,q2

      parameter (small=0.00001,large=10000)

! comment one of next 2 lines out
! 1st is for power-of-2 FFTs
! 2nd is for non-power-of-2 FFTs

!      data nfacs /1/, factors /2,0,0/
      data nfacs /3/, factors /2,3,5/

! if problems arise with the "new" method, set usenew = .false.
      usenew = .false. ! gjp changed from LAMMPS true to use "old method" instead

      if (usenew) then

        q2 = qsqsum/dielectric

! see JCP 109, pg. 7698 for derivation of coefficients
! higher order coefficients may be computed if needed

        acons(1,0) = 2./3.
        acons(2,0) = 1./50.
        acons(2,1) = 5./294.
        acons(3,0) = 1./588.
        acons(3,1) = 7./1440.
        acons(3,2) = 21./3872.
        acons(4,0) = 1./4320.
        acons(4,1) = 3./1936.
        acons(4,2) = 7601./2271360.
        acons(4,3) = 143./28800.
        acons(5,0) = 1./23232.
        acons(5,1) = 7601./13628160.
        acons(5,2) = 143./69120.
        acons(5,3) = 517231./106536960.
        acons(5,4) = 106640677./11737571328.
        acons(6,0) = 691./68140800.
        acons(6,1) = 13./57600.
        acons(6,2) = 47021./35512320.
        acons(6,3) = 9694607./2095994880.
        acons(6,4) = 733191589./59609088000.
        acons(6,5) = 326190917./11700633600.
        acons(7,0) = 1./345600.
        acons(7,1) = 3617./35512320.
        acons(7,2) = 745739./838397952.
        acons(7,3) = 56399353./12773376000.
        acons(7,4) = 25091609./1560084480.
        acons(7,5) = 1755948832039./36229939200000.
        acons(7,6) = 4887769399./37838389248.

! error check, if higher order desired, add needed acons values above,
! or set usenew = .false.

        if (orderflag > 7) &
     &      call error('PPPM order cannot be greater than 7')

! compute necessary gewald
! based on desired error and real space cutoff
! fluid-occupied volume used to estimate real-space error

        gewald = sqrt(-log(long_prec* &
     &    sqrt(natoms*cutcoul*xprd*yprd*zprd)/(2*q2)))/cutcoul

! compute optimal nx_pppm,ny_pppm,nz_pppm
! based on assignment order and desired error

        h = 1.0
        h1 = 2.0
        ncount = 0
        er = large
        er1 = 0
        do while (abs(er) > small)
           summation = 0.0
           do m = 0,orderflag-1
              summation = summation + acons(orderflag,m)*(h*gewald)**(2*m)
           end do
           lpr = q2*(h*gewald)**orderflag * &
     &          sqrt(gewald*xprd*sqrt(2.0D0*3.141592654)*summation/natoms)/ &
     &          (xprd*xprd)
           er = log(lpr) - log(long_prec)
           er2 = er1
           er1 = er
           h2 = h1
           h1 = h
           if ((er1 - er2) == 0.0) then
             h = h1 + er1
           else
             h = h1 + er1*(h2 - h1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call error('can''t get x-mesh spacing')
        end do
        nx_pppm = xprd/h + 1

        ncount = 0
        er = large
        er1 = 0
        do while (abs(er) > small)
           summation = 0.0
           do m = 0,orderflag-1
              summation = summation + acons(orderflag,m)*(h*gewald)**(2*m)
           end do
           lpr = q2*(h*gewald)**orderflag * &
     &          sqrt(gewald*yprd*sqrt(2.0D0*3.141592654)*summation/natoms)/ &
     &          (yprd*yprd)
           er = log(lpr) - log(long_prec)
           er2 = er1
           er1 = er
           h2 = h1
           h1 = h
           if ((er1 - er2) == 0.0) then
             h = h1 + er1
           else
             h = h1 + er1*(h2 - h1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call error('can''t get y-mesh spacing')
        end do
        ny_pppm = yprd/h + 1

        ncount = 0
        er = large
        er1 = 0
        do while (abs(er) > small)
           summation = 0.0
           do m = 0,orderflag-1
              summation = summation + acons(orderflag,m)*(h*gewald)**(2*m)
           end do
           lpr = q2*(h*gewald)**orderflag * &
     &          sqrt(gewald*zprd_slab*sqrt(2.0*3.141592654)*summation/natoms)/ &
     &          (zprd_slab*zprd_slab)
           er = log(lpr) - log(long_prec)
           er2 = er1
           er1 = er
           h2 = h1
           h1 = h
           if ((er1 - er2) == 0.0) then
             h = h1 + er1
           else
             h = h1 + er1*(h2 - h1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call error('can''t get z-mesh spacing')
        end do
        nz_pppm = zprd_slab/h + 1

! use the old method to find the needed parameters

      else

! error check - should be maxorder, but has bug for > 6

!        if (orderflag > maxorder)
!     $       call error('PPPM order cannot be larger than maxorder')
        if (orderflag > 6) &
     &       call error('PPPM order cannot be larger than 6')

! compute necessary gewald based on desired error and real space cutoff
! compute optimal nx_pppm,ny_pppm,nz_pppm based on assignment order
! and desired error

        gewald = sqrt(2.0)*(1.35-.15*log(long_prec))/cutcoul
        gdel = (long_prec/(2.5*exp(-a(orderflag))))**(1.0/b(orderflag))
        nx_pppm = gewald*xprd/gdel
        ny_pppm = gewald*yprd/gdel
        nz_pppm = gewald*zprd_slab/gdel
        lpr = long_prec

      endif

! convert PPPM grid size into "good" numbers - factorable by 2 or 2,3,5

      do while (.not.factorable(nx_pppm,nfacs,factors))
        nx_pppm = nx_pppm + 1
      enddo
      do while (.not.factorable(ny_pppm,nfacs,factors))
        ny_pppm = ny_pppm + 1
      enddo
      do while (.not.factorable(nz_pppm,nfacs,factors))
        nz_pppm = nz_pppm + 1
      enddo

! overwrite PPPM grid size with input values

      if (meshflag.eq.1) then
        nx_pppm = nx_pppm_input
        ny_pppm = ny_pppm_input
        nz_pppm = nz_pppm_input
      endif

! print-out actual PPPM grid to be used

      write (msgstr,*) 'Actual PPPM grid =', nx_pppm,ny_pppm,nz_pppm
      call log_msg(msgstr)

      if (usenew) then

! now that the actual grid has been determined, find best ewald param

        gew1 = gewald + 0.01
        ncount = 0
        er = large
        er1 = 0
        hx = xprd/nx_pppm
        hy = yprd/ny_pppm
        hz = zprd_slab/nz_pppm
        do while (abs(er) > small)
           summation = 0.0
           do m = 0,orderflag-1
              summation = summation + acons(orderflag,m)*(hx*gewald)**(2*m)
           enddo
           lprx = q2*(hx*gewald)**orderflag * &
     &          sqrt(gewald*xprd*sqrt(2.0D0*3.141592654)*summation/natoms)/ &
     &          (xprd*xprd)

           summation = 0.0
           do m = 0,orderflag-1
              summation = summation + acons(orderflag,m)*(hy*gewald)**(2*m)
           enddo
           lpry = q2*(hy*gewald)**orderflag * &
     &          sqrt(gewald*yprd*sqrt(2.0D0*3.141592654)*summation/natoms)/ &
     &          (yprd*yprd)

           summation = 0.0
           do m = 0,orderflag-1
              summation = summation + acons(orderflag,m)*(hz*gewald)**(2*m)
           enddo
           lprz = q2*(hz*gewald)**orderflag * &
     &          sqrt(gewald*zprd_slab*sqrt(2.0D0*3.141592654)* &
     &          summation/natoms)/(zprd_slab*zprd_slab)
           lpr = sqrt(lprx*lprx + &
     &                lpry*lpry + lprz*lprz)/sqrt(3.)
           spr = 2.*q2*exp(-gewald*gewald*cutcoul*cutcoul)/ &
     &              sqrt(natoms*cutcoul*xprd*yprd*zprd_slab)
           er = log(lpr) - log(spr)
           er2 = er1
           er1 = er
           gew2 = gew1
           gew1 = gewald
           if ((er1 - er2) == 0.0) then
             gewald = gew1 + er1
           else
             gewald = gew1 + er1*(gew2 - gew1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call error('Trouble finding Ewald G')
        end do

! h*gewald must be small in order for this "new" error estimation routine
!  to be completely valid. See JCP 109, 7698. However, this routine
!  may produce a reasonable guess even if h*gewald > 1

        if ((hx*gewald > 1) .or. (hy*gewald > 1) .or. &
     &       (hz*gewald > 1)) then
          write (msgstr,*) 'WARNING: Invalid assumption ', &
     &         'that h*gewald is small was used in choice of PPPM G'
          call log_msg(msgstr)
        endif

      endif

      write (msgstr,*) 'PPPM G =',gewald
      call log_msg(msgstr)
      write (msgstr,*) 'Expected RMS precision =', lpr
      call log_msg(msgstr)

      return
      end


! -------------------------------------------------------------------------
! gjp new code based on make_greensfn below
! set up NON-ADAPTED (Hockney-Eastwood) Coulomb Green's Fn
!        (replaces 1/k**2) for k-vectors on this PE
!        NB: no 4pi in numerator
!gjp arg list kept the same as make_greensfn below
!     INPUT:
!         na_ord=order of assignment scheme
!         gewald=(G-Ewald) width parameter for Gaussians
!         nx_pppm,ny_pppm,nz_pppm=grid array dimensions
!         cellz,celly,cellz=grid sizes
!         nx1,nx2 = lower, upper x limit of k-space arrays for this PE
!         ny1,ny2 = lower, upper y limit of k-space arrays for this PE
!         nz1,nz2 = lower, upper z limit of k-space arrays for this PE
!     OUTPUT:
!         greensfn=non-modified Coulomb Green's Fn for k-vectors on this PE

      subroutine make_nonadapted_greensfn(na_ord,gewald,nx_pppm,ny_pppm, &
     &     nz_pppm,cellx,celly,cellz,nx1,nx2,ny1,ny2,nz1,nz2,greensfn)
      implicit none

      integer na_ord, nx_pppm,ny_pppm,nz_pppm ! gobal mesh
      integer nx1,nx2,ny1,ny2,nz1,nz2
      real*8 gewald, cellx,celly,cellz  ! cell dimensions 
      real*8 greensfn ! 3D greesfn for output
      dimension greensfn(*)

!     lcaol variables
      integer indx,m,mper,l,lper,k,kper
      real*8 unitkx,unitky,unitkz,sqk

      unitkx = (2.*3.141592654/cellx)
      unitky = (2.*3.141592654/celly)
      unitkz = (2.*3.141592654/cellz)
      indx = 1
      do m = nz1,nz2
        mper = m-nz_pppm*int(2*m/nz_pppm)
        do l = ny1,ny2
          lper = l-ny_pppm*int(2*l/ny_pppm)
          do k = nx1,nx2
            kper = k-nx_pppm*int(2*k/nx_pppm)
            sqk = ((unitkx*kper)**2+(unitky*lper)**2+(unitkz*mper)**2)
            if (sqk.ne.0.0) then
!gjp NB no 4pi in numerator, as this is electrostic greens fn, not gravitational
              greensfn(indx) = 1.0/sqk
            else
              greensfn(indx) = 0.0
            endif
            indx = indx+1
          enddo
        enddo
      enddo

      return

      end subroutine make_nonadapted_greensfn

      subroutine make_greensfn(na_ord,gewald,nx_pppm,ny_pppm, &
     &     nz_pppm,cellx,celly,cellz,nx1,nx2,ny1,ny2,nz1,nz2,greensfn)
      implicit none

! argument variables

      integer na_ord,nx_pppm,ny_pppm,nz_pppm
      integer nx1,nx2,ny1,ny2,nz1,nz2
      real*8 gewald,cellx,celly,cellz,greensfn
      dimension greensfn(*)

! local variables

      integer nbx,nby,nbz,indx,m,mper,l,lper,k,kper,nx,ny
      integer nz
      real*8 epshoc,shpfn,arg,unitkx,unitky,unitkz,form
      real*8 snz2,sny2,snx2,sqk,sqk_perp,sum1,qx,sx,wx
      real*8 argx,qy,sy,wy,argy,qz,sz,wz,argz,dot1,dot2
      real*8 sum2sq
      data epshoc /1.0D-7/

! implicit function

      shpfn(arg) = exp(-.25*(arg/gewald)**2)

      unitkx = (2.*3.141592654/cellx)
      unitky = (2.*3.141592654/celly)
      unitkz = (2.*3.141592654/cellz)

      nbx = (gewald*cellx/(3.141592654*nx_pppm))*((-log(epshoc))**.25)
      nby = (gewald*celly/(3.141592654*ny_pppm))*((-log(epshoc))**.25)
      nbz = (gewald*cellz/(3.141592654*nz_pppm))*((-log(epshoc))**.25)

      form = 1.0
      indx = 1
      do m = nz1,nz2
        mper = m-nz_pppm*int(2*m/nz_pppm)
        snz2 = (sin(.5*unitkz*mper*cellz/nz_pppm))**2
        do l = ny1,ny2
          lper = l-ny_pppm*int(2*l/ny_pppm)
          sny2 = (sin(.5*unitky*lper*celly/ny_pppm))**2
          do k = nx1,nx2

            kper = k-nx_pppm*int(2*k/nx_pppm)
            snx2 = (sin(.5*unitkx*kper*cellx/nx_pppm))**2
            sqk = ((unitkx*kper)**2+(unitky*lper)**2+(unitkz*mper)**2)
            sqk_perp = ((unitkx*kper)**2+(unitky*lper)**2)

            if (sqk.ne.0.0) then
              greensfn(indx) = form*12.5663706/sqk
              sum1 = 0.0
              do nx = -nbx,nbx
                qx = unitkx*(kper+nx_pppm*nx)
                sx = shpfn(qx)
                wx = 1.0
                argx = .5*qx*cellx/nx_pppm
                if (argx.ne.0.0) wx = (sin(argx)/argx)**na_ord
                do ny = -nby,nby
                  qy = unitky*(lper+ny_pppm*ny)
                  sy = shpfn(qy)
                  wy = 1.0
                  argy = .5*qy*celly/ny_pppm
                  if (argy.ne.0.0) wy = (sin(argy)/argy)**na_ord
                  do nz = -nbz,nbz
                    qz = unitkz*(mper+nz_pppm*nz)
                    sz = shpfn(qz)
                    wz = 1.0
                    argz = .5*qz*cellz/nz_pppm
                    if (argz.ne.0.0) wz = (sin(argz)/argz)**na_ord

                    dot1 = unitkx*kper*qx + &
     &                   unitky*lper*qy+unitkz*mper*qz
                    dot2 = qx*qx+qy*qy+qz*qz
                    sum1 = sum1+(dot1/dot2)*sx*sy*sz*(wx*wy*wz)**2
                  enddo
                enddo
              enddo
              call gf_denom(na_ord,snx2,sny2,snz2,sum2sq)
              greensfn(indx) = greensfn(indx)*sum1/sum2sq
            else
              greensfn(indx) = 0.0
            endif

            indx = indx+1
          enddo
        enddo
      enddo

      return
      end

! -------------------------------------------------------------------------
! returns denominator for Hockney-Eastwood Green's function
!               inf                 n-1
!      S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
!              j=-inf               l=0
!
!       = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z=sin(x)
!     INPUT:
!           n ====== order of assignment scheme
!           x,y,z=== sin(kx*deltax/2), etc.
!           b ====== denominator expansion coefficients /Gamma(2n)
!     OUTPUT:
!           s ====== denominator of Hockney-Eastwood expression

      subroutine gf_denom(n,x,y,z,s)
      implicit none

! argument variables

      integer n
      real*8 x,y,z,s

! local variables

      integer maxorder,nlast,k,ifact,l,m
      real*8 b,gaminv,sx,sy,sz
      parameter (maxorder=16)
      dimension b(0:maxorder-1)
      save b
      data nlast /0/
      save nlast

      if (n.ne.nlast) then
        ifact = 1
        do k = 1,2*n-1
          ifact = ifact*k
        enddo
        gaminv = 1./float(ifact)
        do l = 0,n-1
          b(l) = 0.0
        enddo
        b(0) = 1.0

        do m = 1,n-1
          do l = m,1,-1
            b(l) = 4.0*(b(l)*(l-m)*(l-m-.5)-b(l-1)*(l-1-m)**2)
          enddo
          b(0) = 4.0*(b(0)*(l-m)*(l-m-.5))
        enddo
        do l = 0,n-1
          b(l) = b(l)*gaminv
        enddo
        nlast = n
      endif

      sx = 0.0
      sy = 0.0
      sz = 0.0
      do l = n-1,0,-1
        sx = b(l) + sx*x
        sy = b(l) + sy*y
        sz = b(l) + sz*z
      enddo
      s = sx*sy*sz
      s = s*s

      return
      end

! -------------------------------------------------------------------------
! check if all the factors of n are in factors(1:nfacs)
! return TRUE if they are
! return FALSE if they are not

      logical function factorable(n,nfacs,factors)
      implicit none

! argument variables

      integer n,nfacs,factors(nfacs)

! local variables

      integer nremain,i,j,iflag
      integer nlist,list(20)

! build list of factors of n
! nlist = # of factors
! after inner while loop, i is a factor of n
! nremain = remaining portion of n that has not been factored

      nlist = 0
      nremain = n

      do while (nremain.gt.1)
        i = 2
        do while (mod(nremain,i).ne.0)
          i = i + 1
        enddo
        nlist = nlist + 1
        list(nlist) = i
        nremain = nremain/i
      enddo

! check if every factor of n in list is also is in input factors
! as soon as one list item is not in factors, return FALSE

      do i = 1,nlist
        iflag = 0
        do j = 1,nfacs
          if (list(i).eq.factors(j)) iflag = 1
        enddo
        if (iflag.eq.0) then
          factorable = .FALSE.
          return
        endif
      enddo

      factorable = .TRUE.

      return
      end

! -------------------------------------------------------------------------
! deallocate all PPPM memory, called at end of a run

      subroutine pppm_deallocate
      use lammps_global_module
      implicit none

      if (allocated(density_brick)) then
        deallocate (density_brick)
!gjp phi_brick is new arravy to LAMMPS code
        deallocate (phi_brick)
        deallocate (vdx_brick,vdy_brick,vdz_brick)
        deallocate (density_fft,greensfn)
        deallocate (workvec1,workvec2,vg)
!gjp workvec is new and is necessary when using the FFT E solver
        deallocate (workvec_E)
        deallocate (partgrid,ek)
        deallocate (fkvecs_x,fkvecs_y,fkvecs_z)
        deallocate (pbuf1,pbuf2)
      endif

      return
      end

#endif

