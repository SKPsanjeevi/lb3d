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
! PPPM driver

      subroutine pppm(iflag)
      use lammps_global_module
      use mpi

      use lbe_parallel_module, only: comm_cart

      implicit none

! argument variables
      integer iflag

!gjp addition logical flag 
      logical phi_and_E

! local variables
      integer i,ierror
      real*8 time1,time2,tmp,pi,virtmp(6)

! find grid points for all my particles

      call particle_map

! map my particle charge onto my local 3-d density grid

      time1 = mpi_wtime()
      call make_rho(density_brick)
      time2 = mpi_wtime()
      time_rho = time_rho + time2-time1

! all procs communicate density from their 3-d brick (with ghosts)
!  to others FFT columns - result is filled/summed grid in FFT decomposition

      call brick2fft(density_brick,density_fft)

! compute potential gradient on my FFT grid and
!   portion of e_long on this proc's FFT grid
! return gradients in 3-d brick decomposition

!gjp phi_and_E = .false. then compute phi only (new gjp code)
!gjp phi_and_E = .true., then compute Ex, Ey and Ez as well (original code)
      phi_and_E=.true.
      time1 = mpi_wtime()
      call poisson(phi_brick,vdx_brick,vdy_brick,vdz_brick,iflag, phi_and_E)
      time2 = mpi_wtime()
      time_poiss = time_poiss + time2-time1

! all procs communicate electric fields from their FFT columns 
!  others 3-d bricks (with ghosts) - result is filled grid in 3-d decomp

      call fillbrick(vdx_brick,vdy_brick,vdz_brick)

! calculate the e-fields on my particles

      time1 = mpi_wtime()
      call electric_field(vdx_brick,vdy_brick,vdz_brick)
      time2 = mpi_wtime()
      time_field = time_field + time2-time1

! convert e-field to force

      do i = 1,nlocal
        f(1,i) = f(1,i) + q(i)*coulpre*ek(1,i)
        f(2,i) = f(2,i) + q(i)*coulpre*ek(2,i)
        f(3,i) = f(3,i) + q(i)*coulpre*ek(3,i)
      enddo

! sum long-range component of energy

      tmp = e_long
      call mpi_allreduce(tmp,e_long,1,mpi_double_precision,mpi_sum, &
     &     comm_cart,ierror)

      pi = 4.0D0*atan(1.0D0)
      e_long = e_long - gewald*qsqsum/1.772453851D0 - &
     &     0.5D0*pi*qsum*qsum/(gewald*gewald*xprd*yprd*zprd_slab)
      e_long = coulpre*e_long

! sum long-range component of virial

      if (iflag >= 1) then
        call mpi_allreduce(vir_long,virtmp,6, &
     &       mpi_double_precision,mpi_sum,comm_cart,ierror)
        do i = 1,6
          vir_long(i) = coulpre*virtmp(i)
        enddo
      endif

!gjp commented out as not required when employing poisson solver only 
!gjp and slabcorr is not part of LAMMPS' P3M routines 
! correct for slab geometry if needed
!      if (slabflag == 1) call slabcorr
 
      return
      end


! -------------------------------------------------------------------------
! find center grid pt for each of my particles
! check that full stencil for the particle will fit in my 3-d brick
! store central grid pt in partgrid(3,nlocal) array

      subroutine particle_map
      use lammps_global_module
      use mpi

      use lbe_log_module
      use lbe_parallel_module, only: comm_cart

      implicit none

! local variables

      integer icount,i,nx,ny,nz,jcount,ierror
      real*8 delxinv,delyinv,delzinv,pshift,xtmp,ytmp,ztmp
      
      delxinv = nx_pppm/xprd
      delyinv = ny_pppm/yprd
      delzinv = nz_pppm/zprd_slab

      if (mod(orderflag,2) /= 0) then
        pshift = 4096.5
      else
        pshift = 4096.0
      endif

      icount = 0
      do i = 1,nlocal

! (nx,ny,nz) = global coords of grid pt to "lower left" of charge
! current particle coord can be outside global and local box
! add/subtract 4096 to avoid problem of int(-0.75) = 0 when it needs to be -1

        nx = int((x(1,i)-box(1,1))*delxinv+pshift) - 4096
        ny = int((x(2,i)-box(1,2))*delyinv+pshift) - 4096
        nz = int((x(3,i)-box(1,3))*delzinv+pshift) - 4096

! store grid pt

        partgrid(1,i) = nx
        partgrid(2,i) = ny
        partgrid(3,i) = nz

! check that entire stencil around nx,ny,nz will fit in my 3-d brick

        if (nx+nlower < nxlo_out .or. nx+nupper > nxhi_out .or. &
     &       ny+nlower < nylo_out .or. ny+nupper > nyhi_out .or. &
     &       nz+nlower < nzlo_out .or. nz+nupper > nzhi_out) then
          write (msgstr,*) 'Out-of-range Atom:',node,x(1,i),x(2,i),x(3,i), &
     &         nx,nxlo_out,nxhi_out,ny,nylo_out,nyhi_out, &
     &         nz,nzlo_out,nzhi_out
          icount = icount + 1
          call log_msg(msgstr,.true.)
        endif

      enddo

      call mpi_allreduce(icount,jcount,1,mpi_integer,mpi_sum, &
     &     comm_cart,ierror)
      if (jcount > 0) call error('Cannot compute PPPM')

      return
      end


! -------------------------------------------------------------------------
! create discretized "density" on section of global grid due to my particles
! density(x,y,z) = charge "density" at grid points of my 3-d brick
! (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
!  in global grid

      subroutine make_rho(density)
      use lammps_global_module
      implicit none

! argument variables

      real*8 density(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)

! local variables

      integer maxorder,k,j,i,nx,ny,nz,n,mz,m,my,l,mx
      real*8 rho1d,delxinv,delyinv,delzinv,delvolinv
      real*8 qshift,dx,dy,dz,a0,aa,aaa
      parameter (maxorder=16)
      dimension rho1d(-maxorder/2:maxorder/2,3)
      
      delxinv = nx_pppm/xprd
      delyinv = ny_pppm/yprd
      delzinv = nz_pppm/zprd_slab
      delvolinv = delxinv*delyinv*delzinv

      if (mod(orderflag,2).ne.0) then
        qshift = 0.0
      else
        qshift = 0.5
      endif

! clear density array

      do k = nzlo_out,nzhi_out
        do j = nylo_out,nyhi_out
          do i = nxlo_out,nxhi_out
            density(i,j,k) = 0.0
          enddo
        enddo
      enddo

! loop over my charges, add their contribution to nearby grid points
! (nx,ny,nz) = global coords of grid pt to "lower left" of charge
! (dx,dy,dz) = distance to "lower left" grid pt
! (mx,my,mz) = global coords of moving stencil pt

      do i = 1,nlocal

        nx = partgrid(1,i)
        ny = partgrid(2,i)
        nz = partgrid(3,i)
        dx = nx+qshift - (x(1,i)-box(1,1))*delxinv
        dy = ny+qshift - (x(2,i)-box(1,2))*delyinv
        dz = nz+qshift - (x(3,i)-box(1,3))*delzinv

        call calc_rho1d(dx,dy,dz,rho1d(-maxorder/2,1),orderflag)

        a0 = delvolinv * q(i)
        do n = nlower,nupper
          mz = n+nz
          aa = a0*rho1d(n,3)
          do m = nlower,nupper
            my = m+ny
            aaa = aa*rho1d(m,2)
            do l = nlower,nupper
              mx = l+nx
              density(mx,my,mz) = density(mx,my,mz) + aaa*rho1d(l,1)
            enddo
          enddo
        enddo

      enddo

      return
      end


! -------------------------------------------------------------------------
! Poisson solver - FFT based

!gjp first arg is an addition to LAMMPS code 
      subroutine poisson(phi,vdx,vdy,vdz,iflag,phi_and_E)
      use lammps_global_module
      use mpi
      use lbe_parallel_module, only: comm_cart
      use lbe_elec_timer_module
      implicit none

! argument variables

!gjp phi is an addition to LAMMPS code
      real*8 phi(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 vdx(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 vdy(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 vdz(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      integer iflag

!gjp addition logical flag to compute both phi and E
      logical phi_and_E

! local variables

      integer indx,icount,k,j,i,ix,iy,iz
      real*8 scaleinv,eindx

      complex*16 zim
      real*8 zero,minus_one

      zero = 0.0
      minus_one = -1.0
      zim = cmplx(zero,minus_one)

! transform charge density (r -> k) 

      do indx = 1,nfft
        workvec1(indx) = cmplx(density_fft(indx),zero)
      enddo

      call fft_3d(workvec1,workvec1,1,plan1_fft)

! scale by 1/total-grid-pts to get rho(k)
! multiply by Green's function to get V(k)
! e_long = portion of k-space potential energy on this proc's FFT grid
! compute long-range contribution to virial if requested

      scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm)
      e_long = 0.0

      if (iflag == 0) then

        do indx = 1,nfft
          workvec2(indx) = scaleinv*workvec1(indx)
          workvec1(indx) = greensfn(indx)*workvec2(indx)
          e_long = e_long + dble(workvec1(indx)*conjg(workvec2(indx)))
        enddo
        e_long = 0.5*xprd*yprd*zprd_slab * e_long

      else

        do i = 1,6
          vir_long(i) = 0.0
        enddo
        do indx = 1,nfft
          workvec2(indx) = scaleinv*workvec1(indx)
          workvec1(indx) = greensfn(indx)*workvec2(indx)
          eindx = dble(workvec1(indx)*conjg(workvec2(indx)))
          e_long = e_long + eindx
          do i = 1,6
            vir_long(i) = vir_long(i) + eindx*vg(i,indx)
          enddo
        enddo
        e_long = 0.5*xprd*yprd*zprd_slab * e_long
        do i = 1,6
          vir_long(i) = 0.5*xprd*yprd*zprd_slab * vir_long(i)
        enddo

      endif

!gjp compute potential
!gjp if using the FFT E solver, we need to store workvec1
      if(phi_and_E)workvec_E=workvec1
      call fft_3d(workvec1,workvec1,-1,plan2_fft)
      icount = 0
      do iz = nzlo_in,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            phi(ix,iy,iz) = dble(workvec1(icount))
          enddo
        enddo
      enddo

!gjp compute E field using FFTs
      if(phi_and_E)then
        call start_timer(ti_elec_p3m_calc_E)

! compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
! FFT leaves data in 3-d brick decomposition
! copy it into inner portion of vd_xyz arrays

! x direction gradient

        icount = 0
        do k = nzlo_fft,nzhi_fft
          do j = nylo_fft,nyhi_fft
            do i = nxlo_fft,nxhi_fft
              icount = icount + 1
              workvec2(icount) = zim*fkvecs_x(i)*workvec_E(icount)
            enddo
          enddo
        enddo
  
        call fft_3d(workvec2,workvec2,-1,plan2_fft)
  
        icount = 0
        do iz = nzlo_in,nzhi_in
          do iy = nylo_in,nyhi_in
            do ix = nxlo_in,nxhi_in
              icount = icount + 1
              vdx(ix,iy,iz) = dble(workvec2(icount))
            enddo
          enddo
        enddo
  
! y direction gradient

        icount = 0
        do k = nzlo_fft,nzhi_fft
          do j = nylo_fft,nyhi_fft
            do i = nxlo_fft,nxhi_fft
              icount = icount + 1
              workvec2(icount) = zim*fkvecs_y(j)*workvec_E(icount)
            enddo
          enddo
        enddo
  
        call fft_3d(workvec2,workvec2,-1,plan2_fft)
  
        icount = 0
        do iz = nzlo_in,nzhi_in
          do iy = nylo_in,nyhi_in
            do ix = nxlo_in,nxhi_in
              icount = icount + 1
              vdy(ix,iy,iz) = dble(workvec2(icount))
            enddo
          enddo
        enddo
  
! z direction gradient

        icount = 0
        do k = nzlo_fft,nzhi_fft
          do j = nylo_fft,nyhi_fft
            do i = nxlo_fft,nxhi_fft
              icount = icount + 1
              workvec2(icount) = zim*fkvecs_z(k)*workvec_E(icount)
            enddo
          enddo
        enddo
  
        call fft_3d(workvec2,workvec2,-1,plan2_fft)
  
        icount = 0
        do iz = nzlo_in,nzhi_in
          do iy = nylo_in,nyhi_in
            do ix = nxlo_in,nxhi_in
              icount = icount + 1
              vdz(ix,iy,iz) = dble(workvec2(icount))
            enddo
          enddo
        enddo
        call stop_timer(ti_elec_p3m_calc_E)
!gjp  
      endif ! computing E using FFTs

      return
      end


! -------------------------------------------------------------------------
! interpolate grid values to get electric field acting on my particles

      subroutine electric_field(vdx,vdy,vdz)
      use lammps_global_module
      implicit none

! argument variables

      real*8 vdx(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 vdy(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 vdz(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)

! local variables

      integer maxorder,i,nx,ny,nz,n,mz,m,my,l,mx
      real*8 rho1d,delxinv,delyinv,delzinv,qshift,dx,dy,dz
      real*8 z0,y0,x0
      parameter (maxorder=16)
      dimension rho1d(-maxorder/2:maxorder/2,3)

      delxinv = nx_pppm/xprd
      delyinv = ny_pppm/yprd
      delzinv = nz_pppm/zprd_slab

      if (mod(orderflag,2).ne.0) then
        qshift = 0.0
      else
        qshift = 0.5
      endif
      
! clear e-fields

      do i = 1,nlocal
        ek(1,i) = 0.0
        ek(2,i) = 0.0
        ek(3,i) = 0.0
      enddo

! loop over my charges, interpolate electric field from nearby grid points
! (nx,ny,nz) = global coords of grid pt to "lower left" of charge
! (dx,dy,dz) = distance to "lower left" grid pt
! (mx,my,mz) = global coords of moving stencil pt

      do i = 1,nlocal

        nx = partgrid(1,i)
        ny = partgrid(2,i)
        nz = partgrid(3,i)
        dx = nx+qshift - (x(1,i)-box(1,1))*delxinv
        dy = ny+qshift - (x(2,i)-box(1,2))*delyinv
        dz = nz+qshift - (x(3,i)-box(1,3))*delzinv
        call calc_rho1d(dx,dy,dz,rho1d(-maxorder/2,1),orderflag)

        do n = nlower,nupper
          mz = n+nz
          z0 = rho1d(n,3)
          do m = nlower,nupper
            my = m+ny
            y0 = z0*rho1d(m,2)
            do l = nlower,nupper
              mx = l+nx
              x0 = y0*rho1d(l,1)
              ek(1,i) = ek(1,i) - x0*vdx(mx,my,mz)
              ek(2,i) = ek(2,i) - x0*vdy(mx,my,mz)
              ek(3,i) = ek(3,i) - x0*vdz(mx,my,mz)
            enddo
          enddo
        enddo

      enddo

      return
      end


! -------------------------------------------------------------------------
! cacl_rho1d - calculate 1d assignment functions
!   INPUT:
!         dx,dy,dz = distance of particle  from "lower, left" grid point
!         n ======== order of assignment scheme
!   OUTPUT:
!         w = assignment functions for particle (3*n values)

      subroutine calc_rho1d(dx,dy,dz,w,n)
      implicit none

! argument variables

      integer n
      real*8 dx,dy,dz,w

! local variables

      integer maxorder,n_old,j,k,l,kk
      real*8 acoef,a,s

      parameter (maxorder=16)
      dimension a(0:maxorder-1,-maxorder:maxorder), &
     &     acoef(0:maxorder-1,(1-maxorder)/2:maxorder/2), &
     &     w(-maxorder/2:maxorder/2,3)
      save acoef
      data n_old /0/
      save n_old

! generate the coeffients for the weight function of order n
!              (n-1)
!  Wn(x) =     Sum    wn(k,x) , Sum is over every other integer
!           k=-(n-1)
!  For k=-(n-1),-(n-1)+2, ....., (n-1)-2,n-1
!  Note k is odd integers if n is even
!      and  even integers if n is old
!              ---
!             | n-1
!             | Sum a(l,j)*(x-k/2)**l   if abs(x-k/2) < 1/2
!  wn(k,x) = <  l=0
!             |
!             |  0                       otherwise
!              ---
!  The a coeffients are packed into the array acoef to eliminate zeros
!  acoef(l,((k+mod(n+1,2))/2) = a(l,k)

      if (n.ne.n_old) then
        do k = -n,n
          do l = 0,n-1
            a(l,k) = 0.0d0
          enddo
        enddo
        
        a(0,0) = 1.0d0
        do j = 1,n-1
          do k = -j,j,2
            s = 0.0d0
            do l = 0,j-1
              a(l+1,k) = (a(l,k+1)-a(l,k-1))/dble(l+1)
              s = s+0.5d0**(l+1)*(a(l,k-1)+(-1.0d0)**l*a(l,k+1))/ &
     &             dble(l+1)
            enddo
            a(0,k) = s
          enddo
        enddo
        
        kk = (1-n)/2
        do k = -(n-1),n-1,2
          do l = 0,n-1
            acoef(l,kk) = a(l,k)
          enddo
          kk = kk + 1
        enddo
        n_old = n
      endif

! if acoef is already calculated then just form the {w}s

      do k = (1-n)/2,n/2
        w(k,1)  =   0.0
        w(k,2)  =   0.0
        w(k,3)  =   0.0
        do l = n-1,0,-1
          w(k,1) = acoef(l,k) + w(k,1)*dx
          w(k,2) = acoef(l,k) + w(k,2)*dy
          w(k,3) = acoef(l,k) + w(k,3)*dz
        enddo
      enddo
      
      return
      end

#endif

