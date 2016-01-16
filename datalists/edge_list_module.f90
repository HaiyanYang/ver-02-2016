!***************************************!
!   the global list of edges            !
!***************************************!

module edge_list_module
use fedge_module, only: fedge

implicit none
save

type(fedge), allocatable :: edge_list(:)


contains

  subroutine empty_edge_list()

    if(allocated(edge_list)) deallocate(edge_list)

  end subroutine empty_edge_list

end module edge_list_module
