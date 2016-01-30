module matrix_crack_module
!
!  Purpose:
!   this module defines a matrix crack object
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    15/01/16  B. Y. Chen            Original code
!
use parameter_module, only : NDIM, DP, ZERO, PI, HALFCIRC, MIN_ELEM_SIZE, MAX_ELEM_SIZE

implicit none

private

type :: matrix_crack
  private ! hide components from external operation
  ! list of type components:
  ! - coords     : planar coordinates of a node passed by the crack
  real(DP) :: coords(2) = ZERO
end type matrix_crack

type :: ply_crack_list
  private
  type(matrix_crack), allocatable :: cracks(:)
  integer                         :: num_cracks = 0
end type ply_crack_list

interface set
  module procedure set_ply_crack_list
end interface

public :: ply_crack_list, set, add_newcrack, newcrack_ok


contains


  pure subroutine set_ply_crack_list(crack_list, maxncracks)
    type(ply_crack_list), intent(inout) :: crack_list
    integer,              intent(in)    :: maxncracks
    
    allocate(crack_list%cracks(maxncracks))
  
  end subroutine set_ply_crack_list


  pure function cracklist_set(crack_list) result(is_set)
    type(ply_crack_list), intent(in) :: crack_list
    
    logical :: is_set
    
    is_set = .false.
    
    if (allocated(crack_list%cracks)) is_set = .true.
    
  end function cracklist_set


  pure function cracklist_full(crack_list) result(is_full)
    type(ply_crack_list), intent(in) :: crack_list
    
    logical :: is_full
    
    is_full = .false.
    
    if (crack_list%num_cracks == size(crack_list%cracks)) is_full = .true.
  
  end function cracklist_full
  
  
  pure function newcrack_ok (crack_list, point, theta) result(ok)
    type(ply_crack_list), intent(in) :: crack_list
    real(DP),             intent(in) :: point(2)
    real(DP),             intent(in) :: theta
    
    logical            :: ok
    type(matrix_crack) :: crackj
    real(DP)           :: unit_vect(2), distance
    integer            :: j
    
    if (.not. cracklist_set(crack_list)) then
      ok = .false.
      return
    end if
    
    if (cracklist_full(crack_list)) then
      ok = .false.
      return
    end if
    
    if (crack_list%num_cracks == 0) then
      ok = .true.
      return
    end if
    
    ok = .true.
    unit_vect = ZERO
    distance  = ZERO
    j = 0
    
    unit_vect(1) = -sin(theta/HALFCIRC * PI)
    unit_vect(2) =  cos(theta/HALFCIRC * PI)
    
    do j = 1, crack_list%num_cracks
      crackj   = crack_list%cracks(j)
      distance = abs(dot_product(unit_vect, crackj%coords-point))
      !if (distance < MIN_CRACK_SPACING) then
      !  ok = .false.
      !  exit
      !end if
      if ( 0.01*MIN_ELEM_SIZE < distance .and. distance < 2.0*MAX_ELEM_SIZE ) then
        ok = .false.
        exit
      end if
    end do
    
  end function newcrack_ok


  pure subroutine add_newcrack(crack_list, coords)
    type(ply_crack_list),     intent(inout) :: crack_list
    real(DP),                 intent(in)    :: coords(2)
    
    integer :: i
    i = 0
    
    ! make sure that it is OK to add new crack before calling this subroutine
    
    crack_list%num_cracks = crack_list%num_cracks + 1
    
    i = crack_list%num_cracks
    
    crack_list%cracks(i)%coords = coords
  
  end subroutine add_newcrack
  

end module matrix_crack_module