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
    nullify ( item%next )
    return
endif
!
! Just add ITEM after the head of the list 
item%next => head%next
head%next => item
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
        item_next => item%next
        if ( associated ( item_prev ) ) then    ! not first in the list
            item_prev%next => item_next
        else
            head => item_next                   ! first in list
        endif
        if ( .not.associated ( item_next ) ) then   ! last in list
            nullify(item_prev%next)             
        endif
        deallocate(item)
        return
    endif
    item_prev => item
    item => item%next
end do
end subroutine

end module
