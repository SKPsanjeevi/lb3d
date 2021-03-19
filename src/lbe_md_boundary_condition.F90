#include "lbe.h"

module lbe_md_boundary_condition_module
#ifdef MD
    use lbe_globals_module, only: maxpos,minpos,tsize
    use lbe_md_dynamic_module, only: first_inserted_puid
    use lbe_md_globals_module
!!$    use lbe_parallel_module, only: ccoords,last_periodic_z,maxpos_pi,start,tnx&
!!$         &,tny,tnz,tsize_pi
    use lbe_parallel_module, only: start,tnx,tny,tnz
    use lbe_parms_module, only: nx,ny,nz
    implicit none
    private

    !> Periodic boundaries can split up a particle volume in at most 8
    !> contiguous chunks
    integer,parameter :: n_max_local_chunks = 8

    !> type describing one contiguous lattice chunk related to a particle
    type local_chunk_type
       real(kind=rk) :: xc(3)   !<
       integer :: minx(3)       !< minimum lattice position
       integer :: maxx(3)       !< maximum lattice position
!!$       !> maybe replace this with a particle uid some time???
!!$       logical :: periodic_particle
    end type local_chunk_type

!!$    public calculate_chunksize_pi_z,clip_chunks_minz_global&
!!$         &,copy_clip_chunks_to_real_sites&
!!$         &,copy_split_chunks_into_halo_intersections,is_periodic&
!!$         &,local_chunk_type,local_chunks&
!!$         &,md_clip_periodic,n_boxes_pi_downward_z,n_boxes_pi_upward_z&
!!$         &,n_max_local_chunks,rdst,rdstsq
    public closest_image,copy_clip_chunks_to_real_sites&
         &,copy_split_chunks_into_halo_intersections&
         &,local_chunk_type,local_chunks,md_clip_periodic,min_distance&
         &,n_max_local_chunks,rdst,rdstsq

contains

    !> clips chunks to a new minimum z boundary specified in global
    !> coordinates
    !>
    !> It is possible that empty chunks are returned.
    !>
    !> \param[in] nio number of chunks
    !>
    !> \param[in,out] cio chunks
    !>
    !> \param[in] mzg minimu z coordinate in global coordinates
    subroutine clip_chunks_minz_global(nio,cio,mzg)
        integer,intent(in) :: nio
        type(local_chunk_type),intent(inout) ::cio(n_max_local_chunks)
        integer,intent(in) :: mzg
        integer lmzg,i

        lmzg = mzg+1-start(3)

        do i=1,nio
           cio(i)%minx(3) = max(cio(i)%minx(3),lmzg)
        end do
    end subroutine clip_chunks_minz_global

    !> finds the closest periodic image of one position with respect
    !> to another one
    !>
    !> \param[in] pos position for which to find the image
    !>
    !> \param[in] ref reference position
    !>
    !> \returns the image position of pos which is closest to ref
    pure function closest_image(pos,ref)
        real(kind=rk),dimension(3) :: closest_image
        real(kind=rk),intent(in) :: pos(3),ref(3)
!!$#ifdef PERIODIC_INFLOW
!!$
!!$        call error_md('closest_image() not implemented yet for PERIODIC_INFLOW')
!!$#else
        real(kind=rk) :: dst(3),img(3)
        integer k

        img = pos
        dst = pos-ref
        do k=1,3
           ! "normal" minimum image criterion, result is always valid
           if (abs(dst(k))>0.5_rk*tsize(k)) then
              if (dst(k)<0.0_rk) then
                 img(k) = img(k) + tsize(k)
              else
                 img(k) = img(k) - tsize(k)
              endif
           endif
        end do
        closest_image = img
!!$#endif
    end function closest_image

    !> Returns a clipped copy of a set of chunks so the new one
    !> contains only real (non-halo) sites
    !>
    !> \param[in] ni input number of chunks
    !>
    !> \param[in] ci input chunks
    !>
    !> \param[out] no output number of chunks
    !>
    !> \param[out] co output chunks
    !>
    !> At the moment \c no==ni, so some of the chunks might be empty
    !> but this might change in the future.
    subroutine copy_clip_chunks_to_real_sites(ni,ci,no,co)
        integer,intent(in) :: ni
        type(local_chunk_type),intent(in) ::ci(n_max_local_chunks)
        integer,intent(out) :: no
        type(local_chunk_type),intent(out) ::co(n_max_local_chunks)
        integer :: ls(3)
        integer :: i
!!$#ifdef PERIODIC_INFLOW
!!$        integer :: or(3),ts(3),lor(3),lts(3)
!!$
!!$        if (ni==0) then
!!$           no = 0
!!$           return
!!$        end if
!!$
!!$        ! this assumes that all chunks in ci belong to the same
!!$        ! particle (or at least have the same value of
!!$        ! periodic_particle) and that there is at least one chunk
!!$        if (ci(1)%periodic_particle) then
!!$           or = 1
!!$           ts = (/tnx,tny,last_periodic_z/)
!!$        else
!!$           or = (/1,1,last_periodic_z+1/)
!!$           ts = (/tnx,tny,tnz/)
!!$        end if
!!$        lor = or+1-start
!!$        lts = ts+1-start
!!$#endif

        ls = (/nx,ny,nz/)

        no = 1
        do i=1,ni
           co(no)%minx = max(ci(i)%minx,1)
           co(no)%maxx = min(ci(i)%maxx,ls)
!!$#ifdef PERIODIC_INFLOW
!!$           if (ci(i)%periodic_particle) then
!!$              co(no)%maxx = min(co(no)%maxx,lts)
!!$           else
!!$              co(no)%minx = max(co(no)%minx,lor)
!!$           end if
!!$#endif
           if (all(co(no)%minx<=co(no)%maxx)) then
              co(no)%xc = ci(i)%xc
!!$#ifdef PERIODIC_INFLOW
!!$              co(no)%periodic_particle = ci(i)%periodic_particle
!!$#endif
              no = no+1
           end if
        end do
        no = no-1
    end subroutine copy_clip_chunks_to_real_sites

    !> splits lattice chunks into a set of chunks describing all
    !> intersections of chunks and halo planes
    !>
    !> Halo planes here are defined as planes for which calling \c
    !> update_halo_node() in \c lbe_md_fluid_ladd_module is
    !> desired. This is a bit tricky in case of \c PERIODIC_INFLOW.
    !>
    !> \param[in] ni number of input chunks
    !>
    !> \param[in] ci input chunks
    !>
    !> \param[out] no number of generated output chunks
    !>
    !> \param[out] co generated output chunks
    subroutine copy_split_chunks_into_halo_intersections(ni,ci,no,co)
        integer,intent(in) :: ni
        type(local_chunk_type),intent(in) ::ci(n_max_local_chunks)
        integer,intent(out) :: no
        type(local_chunk_type),intent(out) ::co(6*n_max_local_chunks)
!!$#ifdef PERIODIC_INFLOW
!!$
!!$        if (ni==0) then
!!$           no = 0
!!$           return
!!$        end if
!!$
!!$        ! this assumes that all chunks in ci belong to the same
!!$        ! particle (or at least have the same value of
!!$        ! periodic_particle) and that there is at least one chunk
!!$        if (ci(1)%periodic_particle) then
!!$           call copy_split_chunks_into_halo_intersections_pi_per(ni,ci,no,co)
!!$        else
!!$           call copy_split_chunks_into_halo_intersections_pi_non_per&
!!$                &(ni,ci,no,co)
!!$        end if
!!$#else
        integer :: ls(3),max_mask(3),maxx(3),min_mask(3),minx(3)&
             &,other(3,2),tmp(3),i,k

        ls = (/nx,ny,nz/)

        other(1,:) = (/2,3/)
        other(2,:) = (/1,3/)
        other(3,:) = (/1,2/)

        ! Look at one face of the halo at a time (as far as it has
        ! intersections with one of the images of the chunk). After
        ! one direction is finished, narrow the corresponding chunk
        ! boundaries to the local chunk without halo, so the result
        ! contains no duplicate nodes.
        min_mask = 0
        max_mask = ls+1

        no = 1
        do k=1,3
           do i=1,ni
              minx = max(ci(i)%minx,min_mask)
              maxx = min(ci(i)%maxx,max_mask)
              if ((minx(k)==0).and.(maxx(k)>=0)) then
                 tmp(k) = 0
                 tmp(other(k,1)) = minx(other(k,1))
                 tmp(other(k,2)) = minx(other(k,2))
                 co(no)%minx = tmp
                 tmp(other(k,1)) = maxx(other(k,1))
                 tmp(other(k,2)) = maxx(other(k,2))
                 co(no)%maxx = tmp
                 co(no)%xc = ci(i)%xc
                 no = no+1
              end if
              if ((minx(k)<=ls(k)+1).and.(maxx(k)==ls(k)+1)) then
                 tmp(k) = ls(k)+1
                 tmp(other(k,1)) = minx(other(k,1))
                 tmp(other(k,2)) = minx(other(k,2))
                 co(no)%minx = tmp
                 tmp(other(k,1)) = maxx(other(k,1))
                 tmp(other(k,2)) = maxx(other(k,2))
                 co(no)%maxx = tmp
                 co(no)%xc = ci(i)%xc
                 no = no+1
              end if
           end do
           min_mask(k) = 1
           max_mask(k) = ls(k)
        end do
        no = no-1
!!$#endif
    end subroutine copy_split_chunks_into_halo_intersections

!!$    !> Does the same as \c copy_split_chunks_into_halo_intersections()
!!$    !> but for the case of PERIODIC_INFLOW and periodic particles
!!$    !>
!!$    !> \param[in] ni number of input chunks
!!$    !>
!!$    !> \param[in] ci input chunks
!!$    !>
!!$    !> \param[out] no number of generated output chunks
!!$    !>
!!$    !> \param[out] co generated output chunks
!!$    !>
!!$    !> The main difference to \c
!!$    !> copy_split_chunks_into_halo_intersections() is that the \c z halo
!!$    !> planes are assumed at \c z==0 and \c z==last_periodic_z+1 .
!!$    subroutine copy_split_chunks_into_halo_intersections_pi_per(ni,ci,no,co)
!!$        integer,intent(in) :: ni
!!$        type(local_chunk_type),intent(in) ::ci(n_max_local_chunks)
!!$        integer,intent(out) :: no
!!$        type(local_chunk_type),intent(out) ::co(6*n_max_local_chunks)
!!$        integer :: ls(3),max_mask(3),maxx(3),min_mask(3),minx(3)&
!!$             &,other(3,2),tmp(3),i,k,ltnz_periodic
!!$
!!$        ltnz_periodic = last_periodic_z-(start(3)-1)
!!$
!!$        ! there is no upper z halo plane for periodic particles---this
!!$        ! plane is updated by the non-periodic slaves, still the x and
!!$        ! y edges around z=last_periodic_z are regarded as halo
!!$        ls = (/nx,ny,min(ltnz_periodic-1,nz)/)
!!$
!!$        other(1,:) = (/2,3/)
!!$        other(2,:) = (/1,3/)
!!$        other(3,:) = (/1,2/)
!!$
!!$        ! Look at one face of the halo at a time (as far as it has
!!$        ! intersections with one of the images of the chunk). After
!!$        ! one direction is finished, narrow the corresponding chunk
!!$        ! boundaries to the local chunk without halo, so the result
!!$        ! contains no duplicate nodes.
!!$        min_mask = 0
!!$        max_mask = ls+1
!!$
!!$        no = 1
!!$        do k=1,3
!!$           do i=1,ni
!!$              minx = max(ci(i)%minx,min_mask)
!!$              maxx = min(ci(i)%maxx,max_mask)
!!$              if ((minx(k)==0).and.(maxx(k)>=0)) then
!!$                 tmp(k) = 0
!!$                 tmp(other(k,1)) = minx(other(k,1))
!!$                 tmp(other(k,2)) = minx(other(k,2))
!!$                 co(no)%minx = tmp
!!$                 tmp(other(k,1)) = maxx(other(k,1))
!!$                 tmp(other(k,2)) = maxx(other(k,2))
!!$                 co(no)%maxx = tmp
!!$                 co(no)%xc = ci(i)%xc
!!$                 co(no)%periodic_particle = ci(i)%periodic_particle
!!$                 no = no+1
!!$              end if
!!$              if ((minx(k)<=ls(k)+1).and.(maxx(k)==ls(k)+1)&
!!$                   &.and.(k/=3.or.ls(k)+1/=ltnz_periodic)) then
!!$                 tmp(k) = ls(k)+1
!!$                 tmp(other(k,1)) = minx(other(k,1))
!!$                 tmp(other(k,2)) = minx(other(k,2))
!!$                 co(no)%minx = tmp
!!$                 tmp(other(k,1)) = maxx(other(k,1))
!!$                 tmp(other(k,2)) = maxx(other(k,2))
!!$                 co(no)%maxx = tmp
!!$                 co(no)%xc = ci(i)%xc
!!$                 co(no)%periodic_particle = ci(i)%periodic_particle
!!$                 no = no+1
!!$              end if
!!$           end do
!!$           min_mask(k) = 1
!!$           max_mask(k) = ls(k)
!!$        end do
!!$        no = no-1
!!$    end subroutine copy_split_chunks_into_halo_intersections_pi_per
!!$
!!$    subroutine copy_split_chunks_into_halo_intersections_pi_non_per&
!!$         &(ni,ci,no,co)
!!$        integer,intent(in) :: ni
!!$        type(local_chunk_type),intent(in) ::ci(n_max_local_chunks)
!!$        integer,intent(out) :: no
!!$        type(local_chunk_type),intent(out) ::co(6*n_max_local_chunks)
!!$        integer :: lor(3),ls(3),max_mask(3),maxx(3),min_mask(3),minx(3)&
!!$             &,other(3,2),tmp(3),i,k,ltnz,ltnz_periodic
!!$
!!$        ltnz_periodic = last_periodic_z+1-start(3)
!!$        ltnz = tnz+1-start(3)
!!$
!!$        ! ltnz_periodic does not need to be updated as halo as there
!!$        ! are no periodic z-boundaries, only the edge nodes of the
!!$        ! plane last_periodic_z+1 will be updated as halo here (not
!!$        ! the z-halo itself!)
!!$        lor = (/1,1,max(ltnz_periodic+2,1)/)
!!$
!!$        ! ltnz+1 does not need to be updated as halo as there are no
!!$        ! periodic z-boundaries, only the edge nodes of the plane ltnz
!!$        ! will be updated as halo here (not the z-halo itself!)
!!$        ls = (/nx,ny,min(ltnz-1,nz)/)
!!$
!!$        other(1,:) = (/2,3/)
!!$        other(2,:) = (/1,3/)
!!$        other(3,:) = (/1,2/)
!!$
!!$        ! Look at one face of the halo at a time (as far as it has
!!$        ! intersections with one of the images of the chunk). After
!!$        ! one direction is finished, narrow the corresponding chunk
!!$        ! boundaries to the local chunk without halo, so the result
!!$        ! contains no duplicate nodes.
!!$        min_mask = lor-1
!!$        max_mask = ls+1
!!$
!!$        no = 1
!!$        do k=1,3
!!$           do i=1,ni
!!$              minx = max(ci(i)%minx,min_mask)
!!$              maxx = min(ci(i)%maxx,max_mask)
!!$              if (minx(k)<=lor(k)-1.and.maxx(k)>=lor(k)-1&
!!$                   &.and.(k/=3.or.lor(k)/=ltnz_periodic+1)) then
!!$                 tmp(k) = lor(k)-1
!!$                 tmp(other(k,1)) = minx(other(k,1))
!!$                 tmp(other(k,2)) = minx(other(k,2))
!!$                 co(no)%minx = tmp
!!$                 tmp(other(k,1)) = maxx(other(k,1))
!!$                 tmp(other(k,2)) = maxx(other(k,2))
!!$                 co(no)%maxx = tmp
!!$                 co(no)%xc = ci(i)%xc
!!$                 co(no)%periodic_particle = ci(i)%periodic_particle
!!$                 no = no+1
!!$              end if
!!$              if (minx(k)<=ls(k)+1.and.maxx(k)>=ls(k)+1&
!!$                   &.and.(k/=3.or.nz<ltnz)) then
!!$                 tmp(k) = ls(k)+1
!!$                 tmp(other(k,1)) = minx(other(k,1))
!!$                 tmp(other(k,2)) = minx(other(k,2))
!!$                 co(no)%minx = tmp
!!$                 tmp(other(k,1)) = maxx(other(k,1))
!!$                 tmp(other(k,2)) = maxx(other(k,2))
!!$                 co(no)%maxx = tmp
!!$                 co(no)%xc = ci(i)%xc
!!$                 co(no)%periodic_particle = ci(i)%periodic_particle
!!$                 no = no+1
!!$              end if
!!$           end do
!!$           min_mask(k) = 1
!!$           max_mask(k) = ls(k)
!!$        end do
!!$        no = no-1
!!$    end subroutine copy_split_chunks_into_halo_intersections_pi_non_per

!!$    logical function is_periodic(p)
!!$        type(md_particle_type),intent(in) :: p
!!$
!!$        is_periodic = p%uid<first_inserted_puid
!!$    end function is_periodic

    !> Finds all lattice chunks that (thanks to pbc) contain parts of
    !> a cubic bounding box around a particle
    !>
    !> \param[in[ p particle
    !>
    !> \param[in] r radius around \c p%x to be contained by the chunks
    !>
    !> \param[in] he halo extent
    !>
    !> \param[out] n_chunks number of chunks
    !>
    !> \param[out] chunks chunks
    !>
    !> \c he does not need to match the actual halo extent. For
    !> example, if it is set to \c 0, only real sites will be
    !> returned, even if also halo sites would lie within the bounding
    !> box.
    subroutine local_chunks(p,r,he,n_chunks,chunks)
        type(md_particle_type),intent(in) :: p
        real(kind=rk),intent(in) :: r
        integer,intent(in) :: he
        integer,intent(out) :: n_chunks
        type(local_chunk_type),intent(out) ::chunks(n_max_local_chunks)
        integer :: minb(3),maxb(3)
        integer :: minc(2,3),maxc(2,3)
        real(kind=rk) :: pc(2,3)
        integer :: ls(3),origin(3),ts(3)
        integer :: i,j,k

        ! start and end indices of box around particle in global coordinates
        minb = floor(p%x-r)
        maxb = ceiling(p%x+r)

        ls = (/nx,ny,nz/)

!!$#ifdef PERIODIC_INFLOW
!!$        if (is_periodic(p)) then
!!$           origin = 1
!!$           ts = (/tnx,tny,last_periodic_z/)
!!$
!!$           ! The box might be split into two chunks in each direction
!!$           ! because of pbc. To each chunk belongs a different image of
!!$           ! the particle center. Key to the cartoons: '[xXx]' box around
!!$           ! particle center 'X', '|' global boundaries ( 1 or ts )
!!$           where (minb<origin)       ! |Xx]     [x|
!!$              minc(1,:) = origin-he
!!$              maxc(1,:) = maxb
!!$              pc(1,:) = p%x
!!$              minc(2,:) = minb+ts
!!$              maxc(2,:) = ts+he
!!$              pc(2,:) = p%x+ts
!!$           elsewhere (maxb>ts)  ! |x]     [xX|
!!$              minc(1,:) = minb
!!$              maxc(1,:) = ts+he
!!$              pc(1,:) = p%x
!!$              minc(2,:) = origin-he
!!$              maxc(2,:) = maxb-ts
!!$              pc(2,:) = p%x-ts
!!$           elsewhere            ! |  [xXx]   |
!!$              minc(1,:) = minb
!!$              maxc(1,:) = maxb
!!$              pc(1,:) = p%x
!!$              minc(2,:) = 0   ! disable the loop through the 2nd chunk
!!$              maxc(2,:) = -1
!!$           end where
!!$        else
!!$           origin = (/1,1,last_periodic_z+1/)
!!$           ts = (/tnx,tny,tnz/)
!!$
!!$           ! The box might be split into two chunks in each direction
!!$           ! because of pbc. To each chunk belongs a different image of
!!$           ! the particle center. Key to the cartoons: '[xXx]' box around
!!$           ! particle center 'X', '|' global boundaries ( 1 or ts )
!!$
!!$           ! plain pbc in x and y direction
!!$           where (minb(1:2)<origin(1:2)) ! |Xx]     [x|
!!$              minc(1,1:2) = origin(1:2)-he
!!$              maxc(1,1:2) = maxb(1:2)
!!$              pc(1,1:2) = p%x(1:2)
!!$              minc(2,1:2) = minb(1:2)+ts(1:2)
!!$              maxc(2,1:2) = ts(1:2)+he
!!$              pc(2,1:2) = p%x(1:2)+ts(1:2)
!!$           elsewhere (maxb(1:2)>ts(1:2)) ! |x]     [xX|
!!$              minc(1,1:2) = minb(1:2)
!!$              maxc(1,1:2) = ts(1:2)+he
!!$              pc(1,1:2) = p%x(1:2)
!!$              minc(2,1:2) = origin(1:2)-he
!!$              maxc(2,1:2) = maxb(1:2)-ts(1:2)
!!$              pc(2,1:2) = p%x(1:2)-ts(1:2)
!!$           elsewhere            ! |  [xXx]   |
!!$              minc(1,1:2) = minb(1:2)
!!$              maxc(1,1:2) = maxb(1:2)
!!$              pc(1,1:2) = p%x(1:2)
!!$              minc(2,1:2) = 0   ! disable loop through the 2nd chunk
!!$              maxc(2,1:2) = -1
!!$           end where
!!$
!!$           ! treat only main chunk in z-direction
!!$           if (minb(3)<origin(3)) then ! |Xx]       |
!!$              minc(1,3) = origin(3)-he
!!$              maxc(1,3) = maxb(3)
!!$              pc(1,3) = p%x(3)
!!$           else if (maxb(3)>ts(3)) then ! |       [xX|
!!$              minc(1,3) = minb(3)
!!$              maxc(1,3) = ts(3)+he
!!$              pc(1,3) = p%x(3)
!!$           else                 ! |  [xXx]   |
!!$              minc(1,3) = minb(3)
!!$              maxc(1,3) = maxb(3)
!!$              pc(1,3) = p%x(3)
!!$           end if
!!$           minc(2,3) = 0 ! always disable loop through the 2nd z-chunk
!!$           maxc(2,3) = -1
!!$        end if
!!$#else
        origin = 1
        ts = (/tnx,tny,tnz/)

        ! The box might be split into two chunks in each direction
        ! because of pbc. To each chunk belongs a different image of
        ! the particle center. Key to the cartoons: '[xXx]' box around
        ! particle center 'X', '|' global boundaries ( 1 or ts )
        where (minb<origin)       ! |Xx]     [x|
           minc(1,:) = origin-he
           maxc(1,:) = maxb
           pc(1,:) = p%x
           minc(2,:) = minb+ts
           maxc(2,:) = ts+he
           pc(2,:) = p%x+ts
        elsewhere (maxb>ts)  ! |x]     [xX|
           minc(1,:) = minb
           maxc(1,:) = ts+he
           pc(1,:) = p%x
           minc(2,:) = origin-he
           maxc(2,:) = maxb-ts
           pc(2,:) = p%x-ts
        elsewhere            ! |  [xXx]   |
           minc(1,:) = minb
           maxc(1,:) = maxb
           pc(1,:) = p%x
           minc(2,:) = 0      ! disable the loop through the 2nd chunk
           maxc(2,:) = -1
        end where
!!$#endif

        do i=1,2
!!$#ifdef PERIODIC_INFLOW
!!$           ! global boundaries might differ from local chunk
!!$           ! boundaries - clip to global boundaries
!!$           maxc(i,:) = min(maxc(i,:),ts+he)
!!$           minc(i,:) = max(minc(i,:),origin-he)
!!$
!!$#endif
           ! transform to local coordinates
           minc(i,:) = minc(i,:)+1-start
           maxc(i,:) = maxc(i,:)+1-start
           pc(i,:) = pc(i,:)+real(1-start,kind=rk)

           ! clip to local lattice chunk
           minc(i,:) = max(minc(i,:),1-he)
           maxc(i,:) = min(maxc(i,:),ls+he)
        end do

        ! prepare list of non-nil chunks on the local domain
        n_chunks = 0
        x: do i=1,2
           if (minc(i,1)<=maxc(i,1)) then
              y: do j=1,2
                 if (minc(j,2)<=maxc(j,2)) then
                    z: do k=1,2
                       if (minc(k,3)<=maxc(k,3)) then
                          n_chunks = n_chunks+1
                          chunks(n_chunks)%xc = (/pc(i,1),pc(j,2),pc(k,3)/)
                          chunks(n_chunks)%minx = &
                               &(/minc(i,1),minc(j,2),minc(k,3)/)
                          chunks(n_chunks)%maxx = &
                               &(/maxc(i,1),maxc(j,2),maxc(k,3)/)
!!$#ifdef PERIODIC_INFLOW
!!$                          chunks(n_chunks)%periodic_particle = is_periodic(p)
!!$#endif
                       end if
                    end do z
                 end if
              end do y
           end if
        end do x
    end subroutine local_chunks

    !> Minimum distance vector between two global coordinate positions
    !> taking into account the actual boundary conditions
    !>
    !> \param[in] x1 position 1
    !>
    !> \param[in] x2 position 2
    !>
    !> \returns distance vector from position 2 to position 1
    function min_distance(x1,x2)
        real(kind=rk),dimension(3) :: min_distance
        real(kind=rk),intent(in) :: x1(3),x2(3)
!!$#ifdef PERIODIC_INFLOW
!!$
!!$        call error_md('min_distance() not implemented yet for PERIODIC_INFLOW')
!!$#else
        real(kind=rk) :: r12(3)
        integer k

        r12 = x1-x2
        do k=1,3
           ! "normal" minimum image criterion, result is always valid
           if (abs(r12(k))>0.5_rk*tsize(k)) then
              if (r12(k)<0.0_rk) then
                 r12(k) = r12(k) + tsize(k)
              else
                 r12(k) = r12(k) - tsize(k)
              endif
           endif
        end do
        min_distance = r12
!!$#endif
    end function min_distance

    subroutine md_clip_periodic(k,r)
        integer,intent(in) :: k
        real(kind=rk),intent(inout) :: r

!!$#ifdef PERIODIC_INFLOW
!!$        if (r>=maxpos_pi(k)) then
!!$           r = r-tsize_pi(k)
!!$        else if (r<minpos(k)) then
!!$           r = r+tsize_pi(k)
!!$        end if
!!$#else
        if (r>=maxpos(k)) then
           r = r-tsize(k)
        else if (r<minpos(k)) then
           r = r+tsize(k)
        end if
!!$#endif
    end subroutine md_clip_periodic

!!$    integer function calculate_chunksize_pi_z(cz)
!!$        integer,intent(in) :: cz
!!$        integer c,cdims_pi_z
!!$
!!$        cdims_pi_z = ceiling(real(last_periodic_z,kind=rk)/nz)
!!$
!!$        c = mod(cz,cdims_pi_z)
!!$        do
!!$           if (c>=0) exit
!!$           c = c+cdims_pi_z
!!$        end do
!!$
!!$        if (c==cdims_pi_z-1) then
!!$           calculate_chunksize_pi_z = last_periodic_z-c*nz
!!$        else
!!$           calculate_chunksize_pi_z = nz
!!$        end if
!!$    end function calculate_chunksize_pi_z
!!$
!!$    integer function n_boxes_pi_downward_z(n)
!!$        integer,intent(in) :: n
!!$        integer cdims_pi_z,i,ret
!!$
!!$        cdims_pi_z = ceiling(real(last_periodic_z,kind=rk)/nz)
!!$
!!$        ret = 0
!!$        do i=1,n
!!$           ret = ret+calculate_chunksize_pi_z(ccoords(3)+1-i)
!!$        end do
!!$        n_boxes_pi_downward_z = ret
!!$    end function n_boxes_pi_downward_z
!!$
!!$    integer function n_boxes_pi_upward_z(n)
!!$        integer,intent(in) :: n
!!$        integer cdims_pi_z,i,ret
!!$
!!$        cdims_pi_z = ceiling(real(last_periodic_z,kind=rk)/nz)
!!$
!!$        ret = 0
!!$        do i=1,n
!!$           ret = ret+calculate_chunksize_pi_z(ccoords(3)+i-1)
!!$        end do
!!$        n_boxes_pi_upward_z = ret
!!$    end function n_boxes_pi_upward_z

    !> Distance vector between two particles taking into account the
    !> actual boundary conditions
    !>
    !> An optimized minimum image criterion is used - the result is
    !> useless if its square is greater than \c cutsq1.
    !>
    !> \param[in] p1 particle 1
    !>
    !> \param[in] p2 particle 2
    !>
    !> \returns distance vector from particle 2 to particle 1
    function rdst(p1,p2)
        real(kind=rk),dimension(3) :: rdst
        type(md_particle_type),intent(in) :: p1,p2
!!$#ifdef PERIODIC_INFLOW
!!$
!!$        rdst = rdst_pi(p1,p2)
!!$#else
        real(kind=rk) :: r12(3)
        integer k

        r12 = p1%x-p2%x
        do k=1,3
           if (abs(r12(k))>tsize_mrc(k)) then
              if (r12(k)<0.0_rk) then
                 r12(k) = r12(k) + tsize(k)
              else
                 r12(k) = r12(k) - tsize(k)
              endif
           endif
        end do
        rdst = r12
!!$#endif
    end function rdst

!!$    !> Does the same as rdst but takes into account that in
!!$    !> z-direction, the system is split into a periodic and a
!!$    !> non-periodic sub-volume
!!$    !>
!!$    !> An optimized minimum image criterion is used - the result is
!!$    !> useless if its square is greater than \c cutsq1.
!!$    !>
!!$    !> \param[in] p1 particle 1
!!$    !>
!!$    !> \param[in] p2 particle 2
!!$    !>
!!$    !> \returns distance vector from particle 2 to particle 1
!!$    function rdst_pi(p1,p2)
!!$        real(kind=rk),dimension(3) :: rdst_pi
!!$        type(md_particle_type),intent(in) :: p1,p2
!!$        real(kind=rk) :: r12(3)
!!$        integer k
!!$
!!$        ! no interaction between particles from different sub-volumes
!!$        if (is_periodic(p1).neqv.is_periodic(p2)) then
!!$           rdst_pi = 2.0_rk*rc
!!$           return
!!$        end if
!!$
!!$        r12 = p1%x-p2%x
!!$
!!$        ! no special treatment for x- and y-direction
!!$        do k=1,2
!!$           if (abs(r12(k))>tsize_mrc(k)) then
!!$              if (r12(k)<0.0_rk) then
!!$                 r12(k) = r12(k) + tsize(k)
!!$              else
!!$                 r12(k) = r12(k) - tsize(k)
!!$              endif
!!$           endif
!!$        end do
!!$
!!$        ! z-direction: both particles are either in the periodic or in
!!$        ! the non-periodic sub-volume. Periodic particles experience
!!$        ! periodic wrapping at tsize_periodic(3), for non-periodic
!!$        ! particles there are no periodic z-boundaries.
!!$        if (is_periodic(p1)) then
!!$           if (abs(r12(3))>tsize_mrc_pi(3)) then
!!$              if (r12(3)<0.0_rk) then
!!$                 r12(3) = r12(3) + tsize_pi(3)
!!$              else
!!$                 r12(3) = r12(3) - tsize_pi(3)
!!$              endif
!!$           endif
!!$        end if
!!$
!!$        rdst_pi = r12
!!$    end function rdst_pi

    !> squared distance between two particles taking into account pbc.
    !>
    !> \warning "optimized" minimum image criterion - return value is
    !> useless if it is greater than cutsq1
    !>
    !> \param[in] p1 particle 1
    !>
    !> \param[in] p2 particle 2
    !>
    !> \returns squared distance from particle 2 to particle 1
    function rdstsq(p1,p2)
        real(kind=rk) :: rdstsq
        type(md_particle_type),intent(in) :: p1,p2
        real(kind=rk) :: del(3)

!!$#ifdef PERIODIC_INFLOW
!!$        if (is_periodic(p1).neqv.is_periodic(p2)) then
!!$           rdstsq = 2.0_rk*rc*rc
!!$           return
!!$        end if
!!$
!!$#endif
        del = rdst(p1,p2)
        rdstsq = dot_product(del,del)
    end function rdstsq

#endif
end module lbe_md_boundary_condition_module
