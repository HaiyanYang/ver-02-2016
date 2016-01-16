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
use parameter_module, only : NDIM, DP, ZERO

implicit none

private

type, public :: matrix_crack
  private ! hide components from external operation
  ! list of type components:
  ! - coords     : planar coordinates of a node passed by the crack
  ! - active     : true if this matrix crack is present in the mesh
  real(DP) :: coords(NDIM-1) = ZERO
  logical  :: active = .false.
end type matrix_crack

interface update
  module procedure update_matrix_crack
end interface

interface extract
  module procedure extract_matrix_crack
end interface


public :: update, extract




contains




  pure subroutine update_matrix_crack (this_matrix_crack, coords, active)
  ! Purpose:
  ! to update the components of this matrix_crack; it is used both before and during
  ! analysis to set the initial component values and update the runtime 
  ! component values, respectively.
  
    type(matrix_crack), intent(inout) :: this_matrix_crack
    real(DP),           intent(in)    :: coords(NDIM-1)
    logical,            intent(in)    :: active
    
    this_matrix_crack%coords = coords
    this_matrix_crack%active = active

  end subroutine update_matrix_crack 
  
  

  pure subroutine extract_matrix_crack (this_matrix_crack, coords, active)
  ! Purpose:
  ! to extract all the components of this matrix_crack
  
    type(matrix_crack), intent(in)  :: this_matrix_crack
    real(DP),           intent(out) :: coords(NDIM-1)
    logical,            intent(out) :: active
    
    coords = this_matrix_crack%coords
    active = this_matrix_crack%active

  end subroutine extract_matrix_crack


end module matrix_crack_module