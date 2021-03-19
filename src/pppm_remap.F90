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

! -------------------------------------------------------------------------
! map density in 3-d brick (with ghosts) to density in FFT decomposition
! density_in = 3-d brick decomp
! density_out = FFT decomp

      subroutine brick2fft(density_in,density_out)
      use lammps_global_module
      use mpi
      use lbe_parallel_module, only: comm_cart
      implicit none

! argument variables

      real*8 density_in(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 density_out(*)

! local variables

      integer icount,iz,iy,ix,i,irequest,ierror
      integer istatus(mpi_status_size)

! pack my ghosts for +x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in+1,nxhi_out
            icount = icount + 1
            pbuf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

! pass data to self or +x processor

      if (mpart(2,1).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,1),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,1),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_in,nxlo_in+nxlo_ghost-1
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + pbuf2(icount)
          enddo
        enddo
      enddo

! pack my ghosts for -x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_out,nxlo_in-1
            icount = icount + 1
            pbuf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

! pass data to self or -x processor

      if (mpart(1,1).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,1),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,1),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in-nxhi_ghost+1,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + pbuf2(icount)
          enddo
        enddo
      enddo

! pack my ghosts for +y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in+1,nyhi_out
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            pbuf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

! pass data to self or +y processor

      if (mpart(2,2).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,2),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,2),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_in,nylo_in+nylo_ghost-1
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + pbuf2(icount)
          enddo
        enddo
      enddo

! pack my ghosts for -y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nylo_in-1
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            pbuf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

! pass data to self or -y processor

      if (mpart(1,2).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,2),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,2),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in-nyhi_ghost+1,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + pbuf2(icount)
          enddo
        enddo
      enddo

! pack my ghosts for +z processor

      icount = 0
      do iz = nzhi_in+1,nzhi_out
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            pbuf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

! pass data to self or +z processor

      if (mpart(2,3).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,3),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,3),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_in,nzlo_in+nzlo_ghost-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + pbuf2(icount)
          enddo
        enddo
      enddo

! pack my ghosts for -z processor

      icount = 0
      do iz = nzlo_out,nzlo_in-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            pbuf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

! pass data to self or -z processor

      if (mpart(1,3).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,3),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,3),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack and sum recv data into my real cells

      icount = 0
      do iz = nzhi_in-nzhi_ghost+1,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + pbuf2(icount)
          enddo
        enddo
      enddo

! remap from 3-d brick decomposition to FFT decomposition
! copy done first to grab only inner portion of density from 3-d brick
! remap could be done as pre-stage of FFT,
!  but this one works optimally on only real*8 values, not complex values

      icount = 0
      do iz = nzlo_in,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_out(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

      call remap_3d(density_out,density_out,workvec1,plan_remap)

      return
      end


! ----------------------------------------------------------------------
! fill-in ghost values of potential gradients on each 3-d brick

      subroutine fillbrick(vdx,vdy,vdz)
      use lammps_global_module
      use mpi
      use lbe_parallel_module, only: comm_cart

      implicit none

! argument variables

      real*8 vdx(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 vdy(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)
      real*8 vdz(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)

! local variables

      integer icount,iz,iy,ix,i,irequest,ierror
      integer istatus(mpi_status_size)

! pack my real cells for +z processor

      icount = 0
      do iz = nzhi_in-nzhi_ghost+1,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            pbuf1(icount+2) = vdy(ix,iy,iz)
            pbuf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

! pass data to self or +z processor

      if (mpart(2,3).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,3),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,3),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzlo_in-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            vdy(ix,iy,iz) = pbuf2(icount+2)
            vdz(ix,iy,iz) = pbuf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

! pack my real cells for -z processor

      icount = 0
      do iz = nzlo_in,nzlo_in+nzlo_ghost-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            pbuf1(icount+2) = vdy(ix,iy,iz)
            pbuf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

! pass data to self or +z processor

      if (mpart(1,3).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,3),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,3),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzhi_in+1,nzhi_out
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            vdy(ix,iy,iz) = pbuf2(icount+2)
            vdz(ix,iy,iz) = pbuf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

! pack my real cells for +y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in-nyhi_ghost+1,nyhi_in
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            pbuf1(icount+2) = vdy(ix,iy,iz)
            pbuf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

! pass data to self or +y processor

      if (mpart(2,2).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,2),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,2),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nylo_in-1
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            vdy(ix,iy,iz) = pbuf2(icount+2)
            vdz(ix,iy,iz) = pbuf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

! pack my real cells for -y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_in,nylo_in+nylo_ghost-1
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            pbuf1(icount+2) = vdy(ix,iy,iz)
            pbuf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

! pass data to self or -y processor

      if (mpart(1,2).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,2),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,2),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in+1,nyhi_out
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            vdy(ix,iy,iz) = pbuf2(icount+2)
            vdz(ix,iy,iz) = pbuf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

! pack my real cells for +x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in-nxhi_ghost+1,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            pbuf1(icount+2) = vdy(ix,iy,iz)
            pbuf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

! pass data to self or +x processor

      if (mpart(2,1).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,1),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,1),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_out,nxlo_in-1
            vdx(ix,iy,iz) = pbuf2(icount+1)
            vdy(ix,iy,iz) = pbuf2(icount+2)
            vdz(ix,iy,iz) = pbuf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

! pack my real cells for -x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_in,nxlo_in+nxlo_ghost-1
            pbuf1(icount+1) = vdx(ix,iy,iz)
            pbuf1(icount+2) = vdy(ix,iy,iz)
            pbuf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

! pass data to self or +x processor

      if (mpart(1,1).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,1),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,1),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in+1,nxhi_out
            vdx(ix,iy,iz) = pbuf2(icount+1)
            vdy(ix,iy,iz) = pbuf2(icount+2)
            vdz(ix,iy,iz) = pbuf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

      return
      end

! ----------------------------------------------------------------------
! gjp this routine is based upon fillbrick
! gjp used to redistribute phi only
! fill-in ghost values of phi on each 3-d brick

      subroutine fillbrick_phi(vdx)
      use lammps_global_module
      use mpi
      use lbe_parallel_module, only: comm_cart

      implicit none

! argument variables

      real*8 vdx(nxlo_out:nxhi_out,nylo_out:nyhi_out, &
     &     nzlo_out:nzhi_out)

! local variables

      integer icount,iz,iy,ix,i,irequest,ierror
      integer istatus(mpi_status_size)

! pack my real cells for +z processor

      icount = 0
      do iz = nzhi_in-nzhi_ghost+1,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            icount = icount + 1
          enddo
        enddo
      enddo

! pass data to self or +z processor

      if (mpart(2,3).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,3),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,3),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzlo_in-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            icount = icount + 1
          enddo
        enddo
      enddo

! pack my real cells for -z processor

      icount = 0
      do iz = nzlo_in,nzlo_in+nzlo_ghost-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            icount = icount + 1
          enddo
        enddo
      enddo

! pass data to self or +z processor

      if (mpart(1,3).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,3),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,3),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzhi_in+1,nzhi_out
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            icount = icount + 1
          enddo
        enddo
      enddo

! pack my real cells for +y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in-nyhi_ghost+1,nyhi_in
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            icount = icount + 1
          enddo
        enddo
      enddo

! pass data to self or +y processor

      if (mpart(2,2).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,2),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,2),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nylo_in-1
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            icount = icount + 1
          enddo
        enddo
      enddo

! pack my real cells for -y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_in,nylo_in+nylo_ghost-1
          do ix = nxlo_in,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            icount = icount + 1
          enddo
        enddo
      enddo

! pass data to self or -y processor

      if (mpart(1,2).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,2),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,2),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in+1,nyhi_out
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = pbuf2(icount+1)
            icount = icount + 1
          enddo
        enddo
      enddo

! pack my real cells for +x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in-nxhi_ghost+1,nxhi_in
            pbuf1(icount+1) = vdx(ix,iy,iz)
            icount = icount + 1
          enddo
        enddo
      enddo

! pass data to self or +x processor

      if (mpart(2,1).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(1,1),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(2,1),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_out,nxlo_in-1
            vdx(ix,iy,iz) = pbuf2(icount+1)
            icount = icount + 1
          enddo
        enddo
      enddo

! pack my real cells for -x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_in,nxlo_in+nxlo_ghost-1
            pbuf1(icount+1) = vdx(ix,iy,iz)
            icount = icount + 1
          enddo
        enddo
      enddo

! pass data to self or +x processor

      if (mpart(1,1).eq.node) then
        do i = 1,icount
          pbuf2(i) = pbuf1(i)
        enddo
      else
        call mpi_irecv(pbuf2,maxpbuf,mpi_double_precision, &
     &       mpart(2,1),0,comm_cart,irequest,ierror)
        call mpi_send(pbuf1,icount,mpi_double_precision, &
     &       mpart(1,1),0,comm_cart,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

! unpack recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in+1,nxhi_out
            vdx(ix,iy,iz) = pbuf2(icount+1)
            icount = icount + 1
          enddo
        enddo
      enddo

      return
      end

#endif
