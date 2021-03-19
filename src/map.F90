#include "lbe.h"

!> routines and datatypes implementing an associative array
!>
!> At the moment, this module contains only code necessary to
!> implement something equivalent to (in C++) an \c std::map<int,int>
!> or an \c std::map<int,double[]> with, in the latter case, each
!> value consisting of 6 double values. Thus, all types, subroutines
!> and functions were either given the prefix \c Mii_ or \c Mi6r_.
!>
!> The data is stored in an allocatable array of user defined
!> type. The interface is defined in a way to allow seamless usage of
!> the Mii_ type for uid2i, which is used by the MD routines to map
!> unique MD particle IDs to indices in the local array \c P. The
!> focus should be on fast access to content values based on a given
!> key. The caller has to make sure that indices, keys or the number
!> of inserted elements do not exceed the available index or key
!> values or the capacity of the map. The behavior is undefined if
!> keys are not unique.
module map_module
    use lbe_globals_module, only: rk

    implicit none
    private

    !> one data element of a map
    type Mii_data_type
       integer key              !< unique key
       integer value            !< associated value
    end type Mii_data_type

    !> represents a single map object with integer keys and values
    type Mii_type
       private
       type(Mii_data_type),allocatable :: data(:) !< actual data
       integer count                              !< actual number of elements
       integer position                           !< current data pointer
       integer,allocatable :: remove_list(:)      !< list of keys to remove
       integer n_remove                           !< length of remove_list
    end type Mii_type

    public Mii_at,Mii_clear,Mii_count,Mii_commit,Mii_cur_key,Mii_cur_value&
         &,Mii_destroy,Mii_exist,Mii_exist_value,Mii_forward,Mii_init&
         &,Mii_insert,Mii_invert,Mii_map,Mii_preinsert,Mii_preremove&
         &,Mii_provide_capacity,Mii_purge,Mii_remove,Mii_rewind,Mii_safe_assign&
         &,Mii_step,Mii_type

    !> one data element of a map
    type Mi6r_data_type
       integer key               !< unique key
       real(kind=rk) :: value(6) !< associated value
    end type Mi6r_data_type

    !> represents a single map object with integer keys and values of
    !> real(kind=rk),dimension(6)
    type Mi6r_type
       private
       type(Mi6r_data_type),allocatable :: data(:) !< actual data
       integer count                              !< actual number of elements
       integer position                           !< current data pointer
       integer,allocatable :: remove_list(:)      !< list of keys to remove
       integer n_remove                           !< length of remove_list
    end type Mi6r_type

    public Mi6r_at,Mi6r_clear,Mi6r_count,Mi6r_commit,Mi6r_cur_key&
         &,Mi6r_cur_value,Mi6r_destroy,Mi6r_exist,Mi6r_exist_value,Mi6r_forward&
         &,Mi6r_init,Mi6r_insert,Mi6r_invert,Mi6r_map,Mi6r_preinsert&
         &,Mi6r_preremove,Mi6r_provide_capacity,Mi6r_purge,Mi6r_remove&
         &,Mi6r_rewind,Mi6r_safe_assign,Mi6r_step,Mi6r_type

contains

    !> access on a value based not on the respective key but of the
    !> internal position within \c data.
    !>
    !> \param[in] m map object
    !>
    !> \param[in] p index (starting with 1)
    !>
    !> \returns resulting value
    integer function Mii_at(m,p)
        type(Mii_type),intent(in) :: m
        integer,intent(in) :: p

        Mii_at = m%data(p)%value
    end function Mii_at

    !> empty a map
    !>
    !> \note This does not free the memory allocated for \c m.
    !> \param[in,out] m map object
    subroutine Mii_clear(m)
        type(Mii_type),intent(inout) :: m

        m%count = 0
        m%n_remove = 0
        call Mii_rewind(m)
    end subroutine Mii_clear

    !> return the number of elements in a map
    !>
    !> \param[in] m map object
    !>
    !> \returns number of elements
    integer function Mii_count(m)
        type(Mii_type),intent(in) :: m

        Mii_count = m%count
    end function Mii_count

    !> prepare a map for usage after elements have been added with
    !> \c Mii_preinsert()
    !>
    !> \warning This needs to be called before \c m is accessed by any other
    !> means than a further call to \c Mii_preinsert().
    !>
    !> \param[in] m map object
    subroutine Mii_commit(m)
        type(Mii_type),intent(inout) :: m

        m%count = m%position-1
        call Mii_sort(m)
    end subroutine Mii_commit

    !> return the key of the current element
    !>
    !> \param[in] m map object
    !>
    !> \returns current key
    integer function Mii_cur_key(m)
        type(Mii_type),intent(in) :: m

        Mii_cur_key = m%data(m%position)%key
    end function Mii_cur_key

    !> return the value of the current element
    !>
    !> \param[in] m map object
    !>
    !> \returns current value
    integer function Mii_cur_value(m)
        type(Mii_type),intent(in) :: m

        Mii_cur_value = m%data(m%position)%value
    end function Mii_cur_value

    !> frees the memory allocated for a map while leaving it in an
    !> undefined state
    !>
    !> \param[inout] m map object
    subroutine Mii_destroy(m)
        type(Mii_type),intent(inout) :: m

        if (allocated(m%data)) deallocate (m%data)
        if (allocated(m%remove_list)) deallocate (m%remove_list)
    end subroutine Mii_destroy

    !> checks for the existence of a given key in a map
    !>
    !> \param[in] m map object
    !>
    !> \param[in] k key
    !>
    !> \returns \c .true. if \c m contains at least one entry with key \c k,
    !> otherwise \c .false.
    logical function Mii_exist(m,k)
        type(Mii_type),intent(in) :: m
        integer,intent(in) :: k

        if (m%count==0) then
           Mii_exist = .false.
        else
           Mii_exist = m%data(Mii_find_key(m,k))%key==k
        end if
    end function Mii_exist

    !> checks for the existence of a given value in a map
    !>
    !> \param[in] m map object
    !>
    !> \param[in] v value
    !>
    !> \returns \c .true. if \c m contains at least one entry with value \c v,
    !> otherwise \c .false.
    logical function Mii_exist_value(m,v)
        type(Mii_type),intent(in) :: m
        integer,intent(in) :: v

        if (m%count==0) then
           Mii_exist_value = .false.
        else
           Mii_exist_value = m%data(Mii_find_value(m,v))%value==v
        end if
    end function Mii_exist_value

    !> returns the index in \c m%data to the entry with key \c k
    !>
    !> \param[in] m map object
    !>
    !> \param[in] k key to find
    !>
    !> \returns index of element with key \c k in \c m%data
    !>
    !> \warning \c m must contain at least one element!
    integer function Mii_find_key(m,k)
        type(Mii_type),intent(in) :: m
        integer,intent(in) :: k
        integer i,lo,hi

        lo = 1
        hi = m%count
        do while (lo<=hi)
           i = (lo+hi)/2
           if (k<m%data(i)%key) then
              hi = i-1
           else if (k>m%data(i)%key) then
              lo = i+1
           else
              exit
           end if
        end do
        Mii_find_key = i
    end function Mii_find_key

    !> returns the index in \c m%data to an entry with value \c v
    !>
    !> \param[in] m map object
    !>
    !> \param[in] v value to find
    !>
    !> \returns index of an element with value \c v in \c m%data. If
    !> \c v exists more then once, one of their indices is chosen
    !> arbitrarily. If \c v does not exist, an arbitrary yet valid
    !> index is returned.
    !>
    !> \warning \c m must contain at least one element!
    integer function Mii_find_value(m,v)
        type(Mii_type),intent(in) :: m
        integer,intent(in) :: v
        integer i

        do i=1,m%count
           if (m%data(i)%value==v) then
              Mii_find_value = i
              return
           end if
        end do
        Mii_find_value = 1
    end function Mii_find_value

    !> move the data pointer behind the last element
    !>
    !> \param[in,out] m map object
    subroutine Mii_forward(m)
        type(Mii_type),intent(inout) :: m

        m%position = m%count+1
    end subroutine Mii_forward

    !> basic initialization of a newly created map
    !>
    !> \note Call \c Mii_provide_capacity() before you actually insert data.
    !>
    !> \param[in,out] m map object
    subroutine Mii_init(m)
        type(Mii_type),intent(inout) :: m

        call Mii_clear(m)
    end subroutine Mii_init

    subroutine Mii_insert(m,k,v)
        type(Mii_type),intent(inout) :: m
        integer,intent(in) :: k,v

        call Mii_preinsert(m,k,v)
        call Mii_commit(m)
    end subroutine Mii_insert

    !> creates an inverse map where the old values are keys to the old keys
    !>
    !> \param[in] m existing map object
    !>
    !> \param[in] m inverse map object to create
    !>
    !> \returns error code as \c stat argument of \c allocate
    !>
    !> This does not care about possible duplicate keys or values. Any
    !> existing data in \c im is lost.
    integer function Mii_invert(m,im)
        type(Mii_type),intent(in) :: m
        type(Mii_type),intent(inout) :: im
        integer :: i,stat

        call Mii_init(im)

        stat = Mii_provide_capacity(im,Mii_count(m))
        no_error: if (stat==0) then

           ! avoid Mii_step() since this would break intent(in) for m
           do i=1,m%count
              call Mii_preinsert(im,m%data(i)%value,m%data(i)%key)
           end do

           call Mii_commit(im)
        end if no_error

        ! pass error code from Mii_provide_capacity()
        Mii_invert = stat
    end function Mii_invert

    !> return the value associated to a given key
    !>
    !> If \c k does not exist in \c m or is not unique, the behavior is
    !> undefined.
    !>
    !> \param[in] m map object
    !>
    !> \param[in] k key
    !>
    !> \returns associated value
    integer function Mii_map(m,k)
        type(Mii_type),intent(in) :: m
        integer,intent(in) :: k

        Mii_map = m%data(Mii_find_key(m,k))%value
    end function Mii_map

    !> intended to allow for efficient insertion of many data elements
    !>
    !> \c Mii_preinsert() may be called arbitrarily often as long as the
    !> capacity of \c m is sufficiently large. Before further usage,
    !> \c Mii_commit() needs to be called.
    !>
    !> \param[in,out] m map object
    !>
    !> \param[in] k key
    !>
    !> \param[in] v associated value
    subroutine Mii_preinsert(m,k,v)
        type(Mii_type),intent(inout) :: m
        integer,intent(in) :: k,v

        m%data(m%position) = Mii_data_type(k,v)
        m%position = m%position+1
    end subroutine Mii_preinsert

    !> mark elements for deletion by a consecutive call to \c Mii_purge()
    !>
    !> The benefit compared to \c Mii_remove() is that \c
    !> Mii_preremove() can be called safely while iterating through
    !> the map.
    !>
    !> \param[inout] m map object
    !>
    !> \param[in] k key to element to delete
    !>
    !> \warning After \c Mii_preremove(), no further calls except to
    !> \c Mii_preremove() should be issued until \c Mii_purge() was
    !> called.
    integer function Mii_preremove(m,k)
        type(Mii_type),intent(inout) :: m
        integer,intent(in) :: k
        integer,allocatable :: tmp(:)
        integer stat

        if (.not.allocated(m%remove_list)) then
           allocate(m%remove_list(1),stat=stat)
           if (stat/=0) then
              Mii_preremove = stat
              return
           end if
        end if
        if (size(m%remove_list)<m%n_remove+1) then
           allocate(tmp(m%n_remove),stat=stat)
           if (stat/=0) then
              Mii_preremove = stat
              return
           end if
           tmp = m%remove_list(1:m%n_remove)
           deallocate(m%remove_list)
           allocate(m%remove_list(2*m%n_remove),stat=stat)
           if (stat/=0) then
              Mii_preremove = stat
              return
           end if
           m%remove_list(1:m%n_remove) = tmp
           deallocate(tmp)
        end if

        m%n_remove = m%n_remove+1
        m%remove_list(m%n_remove) = k
        Mii_preremove = 0
    end function Mii_preremove

    !> make sure a map has a given capacity
    !>
    !> \param[in,out] m map object
    !>
    !> \param[in] c desired capacity
    !>
    !> \returns error code as \c stat argument of \c allocate
    integer function Mii_provide_capacity(m,c)
        type(Mii_type),intent(inout) :: m
        integer,intent(in) :: c
        integer stat
        type(Mii_data_type),allocatable :: tmp(:)

        if (.not.allocated(m%data)) then
           allocate (m%data(c),stat=stat)
        else if (size(m%data)<c) then
           allocate (tmp(m%count),stat=stat)
           if (stat/=0) then
              Mii_provide_capacity = stat
              return
           end if
           tmp = m%data(1:m%count)
           deallocate (m%data)
           allocate (m%data(c),stat=stat)
           m%data(1:m%count) = tmp
           deallocate (tmp)
        else
           stat = 0 ! nothing to do - also allocate would return 0 on success
        end if

        Mii_provide_capacity = stat
    end function Mii_provide_capacity

    subroutine Mii_purge(m)
        type(Mii_type),intent(inout) :: m
        integer i

        do i=1,m%n_remove
           call Mii_remove(m,m%remove_list(i))
        end do
        m%n_remove = 0
    end subroutine Mii_purge

    subroutine Mii_remove(m,k)
        type(Mii_type),intent(inout) :: m
        integer,intent(in) :: k

        call Mii_remove_element(m,Mii_find_key(m,k))
    end subroutine Mii_remove

    subroutine Mii_remove_element(m,e)
        type(Mii_type),intent(inout) :: m
        integer,intent(in) :: e
        integer :: i

        m%count = m%count-1
        do i=e,m%count
           m%data(i) = m%data(i+1)
        end do
    end subroutine Mii_remove_element

    !> rewind the data pointer back to the first element
    !>
    !> \param[in,out] m map object
    subroutine Mii_rewind(m)
        type(Mii_type),intent(inout) :: m

        m%position = 1
    end subroutine Mii_rewind

    !> assign a value to a specific key, if the key exists already,
    !> change the associated value instead of creating a 2nd entry
    !>
    !> \param[in,out] m map object
    !>
    !> \param[in] k key
    !>
    !> \param[in] k value
    subroutine Mii_safe_assign(m,k,v)
        type(Mii_type),intent(inout) :: m
        integer,intent(in) :: k,v
        integer :: pos

        if (m%count/=0) then
           pos = Mii_find_key(m,k)
           if (m%data(pos)%key==k) then
              m%data(pos)%value = v
           else
              call Mii_insert(m,k,v)
           end if
        else
           call Mii_insert(m,k,v)
        end if
    end subroutine Mii_safe_assign

    !> move the data pointer to the next element
    !>
    !> \param[in,out] m map object
    subroutine Mii_step(m)
        type(Mii_type),intent(inout) :: m

        m%position = m%position+1
    end subroutine Mii_step

    !> [private] sort elements according to their keys
    !>
    !> \param[in,out] m map object
    subroutine Mii_sort(m)
        type(Mii_type),intent(inout) :: m
        type(Mii_data_type) :: d
        integer i,j

        ! basic insertion sort without special optimizations
        do i=2,m%count
           d = m%data(i)
           ! assume everything <i is sorted, now insert d at the
           ! proper place <i while moving larger elements to the end,
           ! thus extending the range of sorted elements to everything
           ! <i+1.
           do j=i-1,1,-1
              if (m%data(j)%key<=d%key) exit
              m%data(j+1) = m%data(j)
           end do
           m%data(j+1) = d
        end do
    end subroutine Mii_sort

    !> access on a value based not on the respective key but of the
    !> internal position within \c data.
    !>
    !> \param[in] m map object
    !>
    !> \param[in] p index (starting with 1)
    !>
    !> \returns resulting value
    function Mi6r_at(m,p)
        real(kind=rk) :: Mi6r_at(6)
        type(Mi6r_type),intent(in) :: m
        integer,intent(in) :: p

        Mi6r_at = m%data(p)%value
    end function Mi6r_at

    !> empty a map
    !>
    !> \note This does not free the memory allocated for \c m.
    !> \param[in,out] m map object
    subroutine Mi6r_clear(m)
        type(Mi6r_type),intent(inout) :: m

        m%count = 0
        m%n_remove = 0
        call Mi6r_rewind(m)
    end subroutine Mi6r_clear

    !> return the number of elements in a map
    !>
    !> \param[in] m map object
    !>
    !> \returns number of elements
    integer function Mi6r_count(m)
        type(Mi6r_type),intent(in) :: m

        Mi6r_count = m%count
    end function Mi6r_count

    !> prepare a map for usage after elements have been added with
    !> \c Mi6r_preinsert()
    !>
    !> \warning This needs to be called before \c m is accessed by any other
    !> means than a further call to \c Mi6r_preinsert().
    !>
    !> \param[in] m map object
    subroutine Mi6r_commit(m)
        type(Mi6r_type),intent(inout) :: m

        m%count = m%position-1
        call Mi6r_sort(m)
    end subroutine Mi6r_commit

    !> return the key of the current element
    !>
    !> \param[in] m map object
    !>
    !> \returns current key
    integer function Mi6r_cur_key(m)
        type(Mi6r_type),intent(in) :: m

        Mi6r_cur_key = m%data(m%position)%key
    end function Mi6r_cur_key

    !> return the value of the current element
    !>
    !> \param[in] m map object
    !>
    !> \returns current value
    function Mi6r_cur_value(m)
        real(kind=rk) :: Mi6r_cur_value(6)
        type(Mi6r_type),intent(in) :: m

        Mi6r_cur_value = m%data(m%position)%value
    end function Mi6r_cur_value

    !> frees the memory allocated for a map while leaving it in an
    !> undefined state
    !>
    !> \param[inout] m map object
    subroutine Mi6r_destroy(m)
        type(Mi6r_type),intent(inout) :: m

        if (allocated(m%data)) deallocate (m%data)
        if (allocated(m%remove_list)) deallocate (m%remove_list)
    end subroutine Mi6r_destroy

    !> checks for the existence of a given key in a map
    !>
    !> \param[in] m map object
    !>
    !> \param[in] k key
    !>
    !> \returns \c .true. if \c m contains at least one entry with key \c k,
    !> otherwise \c .false.
    logical function Mi6r_exist(m,k)
        type(Mi6r_type),intent(in) :: m
        integer,intent(in) :: k

        if (m%count==0) then
           Mi6r_exist = .false.
        else
           Mi6r_exist = m%data(Mi6r_find_key(m,k))%key==k
        end if
    end function Mi6r_exist

    !> checks for the existence of a given value in a map
    !>
    !> \param[in] m map object
    !>
    !> \param[in] v value
    !>
    !> \returns \c .true. if \c m contains at least one entry with value \c v,
    !> otherwise \c .false.
    logical function Mi6r_exist_value(m,v)
        type(Mi6r_type),intent(in) :: m
        real(kind=rk),dimension(6),intent(in) :: v

        if (m%count==0) then
           Mi6r_exist_value = .false.
        else
           Mi6r_exist_value = all(m%data(Mi6r_find_value(m,v))%value==v)
        end if
    end function Mi6r_exist_value

    !> returns the index in \c m%data to the entry with key \c k
    !>
    !> \param[in] m map object
    !>
    !> \param[in] k key to find
    !>
    !> \returns index of element with key \c k in \c m%data
    !>
    !> \warning \c m must contain at least one element!
    integer function Mi6r_find_key(m,k)
        type(Mi6r_type),intent(in) :: m
        integer,intent(in) :: k
        integer i,lo,hi

        lo = 1
        hi = m%count
        do while (lo<=hi)
           i = (lo+hi)/2
           if (k<m%data(i)%key) then
              hi = i-1
           else if (k>m%data(i)%key) then
              lo = i+1
           else
              exit
           end if
        end do
        Mi6r_find_key = i
    end function Mi6r_find_key

    !> returns the index in \c m%data to an entry with value \c v
    !>
    !> \param[in] m map object
    !>
    !> \param[in] v value to find
    !>
    !> \returns index of an element with value \c v in \c m%data. If
    !> \c v exists more then once, one of their indices is chosen
    !> arbitrarily. If \c v does not exist, an arbitrary yet valid
    !> index is returned.
    !>
    !> \warning \c m must contain at least one element!
    integer function Mi6r_find_value(m,v)
        type(Mi6r_type),intent(in) :: m
        real(kind=rk),dimension(6),intent(in) :: v
        integer i

        do i=1,m%count
           if (all(m%data(i)%value==v)) then
              Mi6r_find_value = i
              return
           end if
        end do
        Mi6r_find_value = 1
    end function Mi6r_find_value

    !> move the data pointer behind the last element
    !>
    !> \param[in,out] m map object
    subroutine Mi6r_forward(m)
        type(Mi6r_type),intent(inout) :: m

        m%position = m%count+1
    end subroutine Mi6r_forward

    !> basic initialization of a newly created map
    !>
    !> \note Call \c Mi6r_provide_capacity() before you actually insert data.
    !>
    !> \param[in,out] m map object
    subroutine Mi6r_init(m)
        type(Mi6r_type),intent(inout) :: m

        call Mi6r_clear(m)
    end subroutine Mi6r_init

    subroutine Mi6r_insert(m,k,v)
        type(Mi6r_type),intent(inout) :: m
        integer,intent(in) :: k
        real(kind=rk),dimension(6),intent(in) :: v

        call Mi6r_preinsert(m,k,v)
        call Mi6r_commit(m)
    end subroutine Mi6r_insert

    !> creates an inverse map where the old values are keys to the old keys
    !>
    !> \param[in] m existing map object
    !>
    !> \param[in] m inverse map object to create
    !>
    !> \returns error code as \c stat argument of \c allocate
    !>
    !> This does not care about possible duplicate keys or values. Any
    !> existing data in \c im is lost.
    integer function Mi6r_invert(m,im)
        type(Mi6r_type),intent(in) :: m
        type(Mi6r_type),intent(inout) :: im
        integer :: i,stat

        ! there will probably never be a map
        ! real(kind=rk),dimension(6) -> integer. This just simulates a
        ! memory allocation error.
        Mi6r_invert = -1
    end function Mi6r_invert

    !> return the value associated to a given key
    !>
    !> If \c k does not exist in \c m or is not unique, the behavior is
    !> undefined.
    !>
    !> \param[in] m map object
    !>
    !> \param[in] k key
    !>
    !> \returns associated value
    function Mi6r_map(m,k)
        real(kind=rk) :: Mi6r_map(6)
        type(Mi6r_type),intent(in) :: m
        integer,intent(in) :: k

        Mi6r_map = m%data(Mi6r_find_key(m,k))%value
    end function Mi6r_map

    !> intended to allow for efficient insertion of many data elements
    !>
    !> \c Mi6r_preinsert() may be called arbitrarily often as long as the
    !> capacity of \c m is sufficiently large. Before further usage,
    !> \c Mi6r_commit() needs to be called.
    !>
    !> \param[in,out] m map object
    !>
    !> \param[in] k key
    !>
    !> \param[in] v associated value
    subroutine Mi6r_preinsert(m,k,v)
        type(Mi6r_type),intent(inout) :: m
        integer,intent(in) :: k
        real(kind=rk),dimension(6),intent(in) :: v

        m%data(m%position) = Mi6r_data_type(k,v)
        m%position = m%position+1
    end subroutine Mi6r_preinsert

    !> mark elements for deletion by a consecutive call to \c Mi6r_purge()
    !>
    !> The benefit compared to \c Mi6r_remove() is that \c
    !> Mi6r_preremove() can be called safely while iterating through
    !> the map.
    !>
    !> \param[inout] m map object
    !>
    !> \param[in] k key to element to delete
    !>
    !> \warning After \c Mi6r_preremove(), no further calls except to
    !> \c Mi6r_preremove() should be issued until \c Mi6r_purge() was
    !> called.
    integer function Mi6r_preremove(m,k)
        type(Mi6r_type),intent(inout) :: m
        integer,intent(in) :: k
        integer,allocatable :: tmp(:)
        integer stat

        if (.not.allocated(m%remove_list)) then
           allocate(m%remove_list(1),stat=stat)
           if (stat/=0) then
              Mi6r_preremove = stat
              return
           end if
        end if
        if (size(m%remove_list)<m%n_remove+1) then
           allocate(tmp(m%n_remove),stat=stat)
           if (stat/=0) then
              Mi6r_preremove = stat
              return
           end if
           tmp = m%remove_list(1:m%n_remove)
           deallocate(m%remove_list)
           allocate(m%remove_list(2*m%n_remove),stat=stat)
           if (stat/=0) then
              Mi6r_preremove = stat
              return
           end if
           m%remove_list(1:m%n_remove) = tmp
           deallocate(tmp)
        end if

        m%n_remove = m%n_remove+1
        m%remove_list(m%n_remove) = k
        Mi6r_preremove = 0
    end function Mi6r_preremove

    !> make sure a map has a given capacity
    !>
    !> \param[in,out] m map object
    !>
    !> \param[in] c desired capacity
    !>
    !> \returns error code as \c stat argument of \c allocate
    integer function Mi6r_provide_capacity(m,c)
        type(Mi6r_type),intent(inout) :: m
        integer,intent(in) :: c
        integer stat
        type(Mi6r_data_type),allocatable :: tmp(:)

        if (.not.allocated(m%data)) then
           allocate (m%data(c),stat=stat)
        else if (size(m%data)<c) then
           allocate (tmp(m%count),stat=stat)
           if (stat/=0) then
              Mi6r_provide_capacity = stat
              return
           end if
           tmp = m%data(1:m%count)
           deallocate (m%data)
           allocate (m%data(c),stat=stat)
           m%data(1:m%count) = tmp
           deallocate (tmp)
        else
           stat = 0 ! nothing to do - also allocate would return 0 on success
        end if

        Mi6r_provide_capacity = stat
    end function Mi6r_provide_capacity

    subroutine Mi6r_purge(m)
        type(Mi6r_type),intent(inout) :: m
        integer i

        do i=1,m%n_remove
           call Mi6r_remove(m,m%remove_list(i))
        end do
        m%n_remove = 0
    end subroutine Mi6r_purge

    subroutine Mi6r_remove(m,k)
        type(Mi6r_type),intent(inout) :: m
        integer,intent(in) :: k

        call Mi6r_remove_element(m,Mi6r_find_key(m,k))
    end subroutine Mi6r_remove

    subroutine Mi6r_remove_element(m,e)
        type(Mi6r_type),intent(inout) :: m
        integer,intent(in) :: e
        integer :: i

        m%count = m%count-1
        do i=e,m%count
           m%data(i) = m%data(i+1)
        end do
    end subroutine Mi6r_remove_element

    !> rewind the data pointer back to the first element
    !>
    !> \param[in,out] m map object
    subroutine Mi6r_rewind(m)
        type(Mi6r_type),intent(inout) :: m

        m%position = 1
    end subroutine Mi6r_rewind

    !> assign a value to a specific key, if the key exists already,
    !> change the associated value instead of creating a 2nd entry
    !>
    !> \param[in,out] m map object
    !>
    !> \param[in] k key
    !>
    !> \param[in] k value
    subroutine Mi6r_safe_assign(m,k,v)
        type(Mi6r_type),intent(inout) :: m
        integer,intent(in) :: k
        real(kind=rk),dimension(6),intent(in) :: v
        integer :: pos

        if (m%count/=0) then
           pos = Mi6r_find_key(m,k)
           if (m%data(pos)%key==k) then
              m%data(pos)%value = v
           else
              call Mi6r_insert(m,k,v)
           end if
        else
           call Mi6r_insert(m,k,v)
        end if
    end subroutine Mi6r_safe_assign

    !> move the data pointer to the next element
    !>
    !> \param[in,out] m map object
    subroutine Mi6r_step(m)
        type(Mi6r_type),intent(inout) :: m

        m%position = m%position+1
    end subroutine Mi6r_step

    !> [private] sort elements according to their keys
    !>
    !> \param[in,out] m map object
    subroutine Mi6r_sort(m)
        type(Mi6r_type),intent(inout) :: m
        type(Mi6r_data_type) :: d
        integer i,j

        ! basic insertion sort without special optimizations
        do i=2,m%count
           d = m%data(i)
           ! assume everything <i is sorted, now insert d at the
           ! proper place <i while moving larger elements to the end,
           ! thus extending the range of sorted elements to everything
           ! <i+1.
           do j=i-1,1,-1
              if (m%data(j)%key<=d%key) exit
              m%data(j+1) = m%data(j)
           end do
           m%data(j+1) = d
        end do
    end subroutine Mi6r_sort

end module map_module
