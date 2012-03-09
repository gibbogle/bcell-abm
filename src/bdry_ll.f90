module bdry_linked_list
use global
implicit none

type (boundary_type), pointer :: bdrylist
type (boundary_type), pointer :: checklist

contains

!***************************************************************************** 
! LIST_INSERT inserts ITEM into the linked list pointed to by HEAD.
!
!  Discussion:
!    In the original code the items in the linked list were maintained ascending sorted by 
!    the VALUE component of the individual boundary_types.  The ordering has been dropped.
!    (see linkedlist\linkedlist.f90 for the sorted version)
!  Licensing: This code is distributed under the GNU LGPL license.
!  Original author: John Burkardt
!  Parameters:
!    Input, type ( boundary_type ) pointer :: ITEM, a pointer to the boundary_type
!    to be inserted into the linked list.
!    Input/output, type ( boundary_type ) pointer :: HEAD, a pointer to the first boundary_type
!    in the linked list.
!*****************************************************************************
subroutine bdrylist_insert ( item, head )
type ( boundary_type ), pointer :: head
type ( boundary_type ), pointer :: item
!
!  In the case of an empty list.  
!
if ( .not. associated ( head ) ) then
    head => item
!    nullify ( item%previous )
    nullify ( item%next )
    return
endif
!
! Just add ITEM after the head of the list 
!
!write(*,*) 'head: ',head%site
!write(*,*) 'item: ',item%site
item%next => head%next
head%next => item
!item%previous => head
!write(*,*) 'added'
!write(*,*) 'head: ',head%site
!write(*,*) 'item: ',head%next%site

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine bdrylist_print ( head )
type ( boundary_type ), pointer :: head
integer i
type ( boundary_type ), pointer :: item

i = 0
item => head
do while ( associated ( item )) 
    i = i + 1
    write ( *, * ) i,item%site ! boundary_type components
    item => item%next
enddo
return
end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
logical function bdrylist_present ( site, head )
integer :: site(3)
type ( boundary_type ), pointer :: head, item

item => head
do while ( associated ( item )) 
    if (item%site(1)==site(1) .and. item%site(2)==site(2) .and. item%site(3)==site(3)) then
        bdrylist_present = .true.
        return
    endif
    item => item%next
enddo
bdrylist_present = .false.
end function

!------------------------------------------------------------------------------
! Delete from the list an entry with the given site.
!------------------------------------------------------------------------------
subroutine bdrylist_delete(site, head)
integer :: site(3)
type ( boundary_type ), pointer :: head
type ( boundary_type ), pointer :: item
type ( boundary_type ), pointer :: item_prev
type ( boundary_type ), pointer :: item_next

if ( .not. associated ( head ) ) then
    return
endif
item => head
if (item%site(1) == site(1) .and. item%site(2) == site(2) .and. item%site(3) == site(3)) then
    ! removing head
    head => item%next
!    nullify(head%previous)
!    nullify(item)
!    write(*,*) 'remove head: ',item%site
    deallocate(item)
    if ( associated ( head ) ) then
!        write(*,*) 'new head: ',head%site
    else
        nullify(head)
    endif
    return
endif
nullify(item_prev)
do while ( associated ( item ) )
    if (item%site(1) == site(1) .and. item%site(2) == site(2) .and. item%site(3) == site(3)) then
!        item_prev => item%previous
        item_next => item%next
        if ( associated ( item_prev ) ) then    ! not first in the list
            item_prev%next => item_next
        else
            head => item_next                   ! first in list
        endif
!        if ( associated ( item_next ) ) then    ! not last in the list
!            item_next%previous => item_prev
!        else
        if ( .not.associated ( item_next ) ) then   ! last in list
!            write(*,*) 'Last in list'
            nullify(item_prev%next)             
        endif
!        nullify(item)
        deallocate(item)
        return
    endif
    item_prev => item
    item => item%next
end do
end subroutine

end module

!program main
!
!!*****************************************************************************80
!!
!!! MAIN is the main program for the linked list example.
!!
!!  Licensing:
!!
!!    This code is distributed under the GNU LGPL license.
!!
!!  Modified:
!!
!!    30 December 2007
!!
!!  Author:
!!
!!    John Burkardt
!!
!  implicit none
!
!  integer data_num
!
!  write ( *, '(a)' ) ' '
!  call timestamp ( )
!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) 'LINKED_LIST:'
!  write ( *, '(a)' ) '  FORTRAN90 version.'
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '  Demonstrate how pointers can be used to define'
!  write ( *, '(a)' ) '  and manipulate a linked list.'
!
!  data_num = 10
!  call test01 ( data_num )
!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) 'LINKED_LIST:'
!  write ( *, '(a)' ) '  Normal end of execution.'
!
!  write ( *, '(a)' ) ' '
!  call timestamp ( )
!
!  stop
!end

!subroutine test01 ( data_num )
!
!!*****************************************************************************80
!!
!!! TEST01 uses a linked list to store and sort random data.
!!
!!  Discussion:
!!
!!    This routine requires a user defined linked-list library.
!!
!!  Licensing:
!!
!!    This code is distributed under the GNU LGPL license.
!!
!!  Modified:
!!
!!    30 December 2007
!!
!!  Author:
!!
!!    John Burkardt
!!
!  use linked_list_library
!
!  implicit none
!
!  integer data_num
!  type ( boundary_type ), pointer :: head
!  integer i, idel
!  type ( boundary_type ), pointer :: item
!  real r
!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) 'TEST01'
!  write ( *, '(a)' ) '  Create, one at a time, a sequence of integers.'
!  write ( *, '(a)' ) '  As each integer is created, insert it into a sorted'
!  write ( *, '(a)' ) '  list.  Print the list at the end.'
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '  We carry out this task using a linked list.'
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '  Initial data generation:'
!  write ( *, '(a)' ) ' '
!
!  nullify ( head )
!
!  do i = 1, data_num
!!
!!  Generate a new item.
!!
!    allocate ( item )
!
!    item%generation = i
!
!    call random_number ( harvest = r )
!
!    item%value = int ( 1000.0 * r )
!    nullify ( item%next )
!    nullify ( item%previous )
!
!    write ( *, '(2x,i8,2x,i8)' ) i, item%value
!!
!!  Insert the new item into the linked list.
!!
!    call list_insert ( item, head )
!
!  end do
!!
!!  Print the linked list.
!!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '  Contents of sorted linked list:'
!  write ( *, '(a)' ) ' '
!
!  call list_print ( head )
!  
!  idel = 5
!  call list_delete(idel, head)
!  write(*,*) 'deleted: ',idel
!  call list_print ( head )
!
!  return
!end

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
subroutine timestamp ( )
  implicit none

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end