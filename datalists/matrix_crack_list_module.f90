!***************************************!
!   the global list of matrix_cracks            !
!***************************************!

module matrix_crack_list_module
use matrix_crack_module, only : ply_crack_list

implicit none
save

! list of all matrix cracks in all plies
type(ply_crack_list), allocatable :: lam_crack_list(:)

contains

  subroutine empty_matrix_crack_list()

    if(allocated(lam_crack_list)) deallocate(lam_crack_list)

  end subroutine empty_matrix_crack_list

end module matrix_crack_list_module
