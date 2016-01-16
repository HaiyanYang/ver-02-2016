subroutine set_fnm_edges()            
use fedge_module,     only: update    
use edge_list_module, only: edge_list 
                                      
  integer :: nedge=0                  
                                      
  nedge=96              
  allocate(edge_list(nedge))          
  call update(edge_list(72), constrained=.true.) 
  call update(edge_list(4), constrained=.true.) 

end subroutine set_fnm_edges
