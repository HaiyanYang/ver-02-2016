!***************************************!                             
!   the list of elements                                 
!***************************************!  

module predelam_list_module

implicit none
save

integer, allocatable :: predelam_elems(:)
integer, allocatable :: predelam_interf

contains

  subroutine empty_predelam_list()
  
    if (allocated(predelam_elems))  deallocate(predelam_elems)
    if (allocated(predelam_interf)) deallocate(predelam_interf)
  
  end subroutine empty_predelam_list

end module predelam_list_module