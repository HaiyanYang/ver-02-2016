module fnode_module
!
!  Purpose:
!   this module defines a node object, which contains coordinates and displacement
!
!  Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    15/07/15  B. Y. Chen            Original code
!
use parameter_module, only : INTACT, NDIM, DP, ZERO

implicit none

private

type, public :: fnode ! a node with enriched d.o.f
  private ! hide components from external operation
  ! list of type components:
  ! - x         : coordinates
  ! - u         : displacement vector
  real(DP) :: x(NDIM) = ZERO
  real(DP) :: u(NDIM) = ZERO
end type fnode

interface empty
  module procedure empty_fnode
end interface

interface update
  module procedure update_fnode
end interface

interface extract
  module procedure extract_fnode
end interface

interface operator(+)
  module procedure plus_fnode
end interface

interface operator(-)
  module procedure minus_fnode
end interface

interface operator(*)
  module procedure ratio_fnode
end interface


public :: empty, update, extract, &
        & operator(+), operator(-), operator(*)




contains



  
  elemental subroutine empty_fnode(this_fnode)
  ! Purpose:
  ! to deallocate all the components of this object
  
    type(fnode), intent(inout) :: this_fnode
    
    this_fnode%x     = ZERO
    this_fnode%u     = ZERO

  end subroutine empty_fnode


  
  elemental function plus_fnode(fnode1, fnode2) result(fnode3)
  ! Purpose:
  ! to add two fnode objects, component by component
  
    type(fnode), intent(in) :: fnode1, fnode2
    type(fnode)             :: fnode3
    
    fnode3%x = fnode1%x + fnode2%x
    fnode3%u = fnode1%u + fnode2%u

  end function plus_fnode


  
  elemental function minus_fnode(fnode1, fnode2) result(fnode3)
  ! Purpose:
  ! to subtract two fnode objects (second from first), component by component
  
    type(fnode), intent(in) :: fnode1, fnode2
    type(fnode)             :: fnode3
    
    fnode3%x = fnode1%x - fnode2%x
    fnode3%u = fnode1%u - fnode2%u

  end function minus_fnode


  
  elemental function ratio_fnode(r, fnode1) result(fnode3)
  ! Purpose:
  ! to multiply all the components of an fnode object by a constant
  
    type(fnode), intent(in) :: fnode1
    real(DP),    intent(in) :: r
    type(fnode)             :: fnode3
    
    fnode3%x = r * fnode1%x
    fnode3%u = r * fnode1%u
 
  end function ratio_fnode



  pure subroutine update_fnode (this_fnode, x, u)
  ! Purpose:
  ! to update the components of this fnode; it is used both before and during
  ! analysis to set the initial component values and update the runtime 
  ! component values, respectively. status and error messages are needed to 
  ! flag an error when inputs are not valid
  
    type(fnode),        intent(inout) :: this_fnode
    real(DP), optional, intent(in)    :: x(NDIM), u(NDIM)
     
    if (present(x))     this_fnode%x     = x           
    if (present(u))     this_fnode%u     = u           

  end subroutine update_fnode 
  
  

  pure subroutine extract_fnode (this_fnode, x, u)
  ! Purpose:
  ! to extract all the components of this fnode
  
    type(fnode),          intent(in)  :: this_fnode
    real(DP),   optional, intent(out) :: x(NDIM),  u(NDIM)
    
    if(present(x))     x     = this_fnode%x 
    if(present(u))     u     = this_fnode%u

  end subroutine extract_fnode


end module fnode_module