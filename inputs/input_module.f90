module input_module
use parameter_module, only: DIRLENGTH

implicit none

public

! define global variable for output directory
character(len=DIRLENGTH), save :: indir

contains

  ! input materials
  include 'inputs/fnm_materials.f90'

  subroutine set_fnm_nodes()            
    use parameter_module, only: DP, NDIM, ZERO  
    use node_list_module, only: node_list 
    use fnode_module,     only: update    
                                      
    integer  :: nnode
    integer  :: i
    real(DP) :: x(NDIM), u(NDIM)
    
    nnode = 0
    i = 0  
    x = ZERO
    u = ZERO
  
    open (unit=111, file=trim(indir)//'fnm_nodes.txt', status='old', action='read')
    
    read(111, *) nnode
    allocate(node_list(nnode)) 
    
    do i = 1, nnode
      read(111, *) x
      call update(node_list(i), x=x, u=u)
    end do
    
    close(111)
  
  end subroutine set_fnm_nodes


  subroutine set_fnm_edges()            
    use fedge_module,     only: update    
    use edge_list_module, only: edge_list 
                                      
    integer :: nedge, iedge                 
    
    nedge = 0
    iedge = 0
    
    open (unit=112, file=trim(indir)//'fnm_edges.txt', status='old', action='read')
    
    read(112, *) nedge    
    allocate(edge_list(nedge))
    
    do
      read(112, *) iedge
      if (iedge <= 0) exit
      call update(edge_list(iedge), tie_bcd=.true.)
    end do
    
    close(112)

  end subroutine set_fnm_edges


  subroutine set_fnm_elems()                                 
    use parameter_module,      only: DP, ZERO                        
    use elem_list_module,      only: layup, elem_list,&        
                          & elem_node_connec, elem_edge_connec 
    use fBrickLam_elem_module, only: plyblock_layup, set       
                                                           
    integer :: nelem                                   
    integer :: elnnode                                   
    integer :: elnedge                                   
    integer :: nplyblk     
    real(DP):: angle
    integer :: nplies                                   
    integer, allocatable :: nodecnc(:), edgecnc(:)     
    integer :: i   
  
                                                           
    nelem   = 0                        
    elnnode = 0                              
    elnedge = 0                              
    nplyblk = 0
    angle   = ZERO
    nplies  = 0
    i = 0
    
    open (unit=113, file=trim(indir)//'fnm_elems.txt', status='old', action='read')
    
    read(113, *) nelem, elnnode, elnedge          
    
    allocate(elem_list(nelem))
                   
    allocate(elem_node_connec(elnnode,nelem))  
    allocate(nodecnc(elnnode))
    nodecnc = 0 
        
    allocate(elem_edge_connec(elnedge,nelem))
    allocate(edgecnc(elnedge)) 
    edgecnc = 0 
    
    read(113, *) nplyblk
    allocate(layup(nplyblk))   
                                  
    do i = 1, nplyblk
      read(113, *) angle, nplies
      layup(i)=plyblock_layup(angle=angle,nplies=nplies)
    end do                                              
                                                 
    do i = 1, nelem
      call set(elem_list(i), NPLYBLKS=nplyblk)
      read(113, *) nodecnc
      read(113, *) edgecnc
      elem_node_connec(:,i)=nodecnc(:)
      elem_edge_connec(:,i)=edgecnc(:)
    end do 
    
    close(113)
    
  end subroutine set_fnm_elems
  
  
  ! predelam
  subroutine set_fnm_predelam()         
    use predelam_list_module, only: predelam_elems, predelam_interf 
                                      
    integer :: npdelem, i, iel                  
                                      
    npdelem = 0
    i = 0
    iel = 0
    
    open (unit=114, file=trim(indir)//'fnm_predelam.txt', status='old', action='read')              

    read(114, *) npdelem
    
    if (npdelem > 0) then
    
      allocate(predelam_elems(npdelem))     
      allocate(predelam_interf)             

      read(114, *) predelam_interf
    
      do i = 1, npdelem
        read(114, *) iel
        predelam_elems(i)=iel  
      end do
      
    end if
    
    close(114)
    
  end subroutine set_fnm_predelam
  
  
  ! matrix cracks
  subroutine set_fnm_matrix_crack()         
    use parameter_module,         only: MIN_CRACK_SPACING
    use matrix_crack_module,      only: set
    use matrix_crack_list_module, only: lam_crack_list 
                                      
    integer :: nplyblock, maxncrack, i                 
    
    nplyblock = 0
    maxncrack = 0
    i = 0
    
    open (unit=115, file=trim(indir)//'fnm_matrix_crack.txt', status='old', action='read')              

    ! read user-input min crack spacing parameter
    read(115, *) MIN_CRACK_SPACING
    
    ! read no. of plyblocks in the laminate, and maximum no. of cracks possible for each plyblock
    read(115, *) nplyblock
    read(115, *) maxncrack
    
    if (maxncrack > 0) then
      ! allocate laminate crack list
      allocate(lam_crack_list(nplyblock))
      ! allocate each ply crack list
      do i = 1, nplyblock
        call set(lam_crack_list(i),maxncrack)
      end do
    end if
    
    close(115)
    
  end subroutine set_fnm_matrix_crack
  

end module input_module