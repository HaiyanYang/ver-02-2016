subroutine set_fnm_predelam()         
use predelam_list_module, only: predelam_elems, predelam_interf 
                                      
  integer :: npdelem                  
                                      
  npdelem=1              

  allocate(predelam_elems(npdelem))     
  allocate(predelam_interf)             

  predelam_elems(1)=1

  predelam_interf=1

end subroutine set_fnm_predelam
