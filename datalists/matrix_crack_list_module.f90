!***************************************!
!   the global list of matrix_cracks            !
!***************************************!

module matrix_crack_list_module
use matrix_crack_module, only: matrix_crack

implicit none
save

! list of all matrix cracks in all plies
type(matrix_crack), allocatable :: matrix_crack_list(:,:)

! mininum crack spacing, a parameter defined at input
real(DP) :: min_crack_spacing


contains

  subroutine empty_matrix_crack_list()

    if(allocated(matrix_crack_list)) deallocate(matrix_crack_list)

  end subroutine empty_matrix_crack_list

end module matrix_crack_list_module
