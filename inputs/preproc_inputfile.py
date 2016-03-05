################################################################
################# Preprocessing for the FNM ####################
################################################################
##
##
################################################################
##  so far, only works on models with a single fnm part, meshed
##  with 3D brick-type elements (C3D8(R) or SC8(R)); it prepares
##  a new input file and fnm modules with fBricklam type only.
################################################################
##
##
################################################################
##  Applicable abaqus input file:
##  - multiple normal Abaqus parts allowed
##  - single fnm part (named 'fnm') meshed with either C3D8(R) 
##    or SC8(R) elements
##  - fnm part definition can include *Node, *Element, *Nset and
##    *Elset, elset=predelam
##  - Nsets in the fnm part can have keyword 'tie' in its name
##    so that all the edges involved in this nset have a tie bcd
##  - assembly definition can be relatively generic, except that
##    new nsets of the nodes on the fnm part cannot be defined.
##  - boundary conditions must be applied on the nsets
################################################################



# glb objects defined for FNM
from preproc_classes import*

import math
# operating system module
import os, sys
# shutil module has functions for different kinds of copying
import shutil

# define get integer
def GetInteger2():
    nnn = None
    while nnn != 1 and nnn !=2:
        try:
            s=sys.stdin.readline()
            nnn = int(s)
        except ValueError:
            print ('Invalid Number. Enter integer. ')
    return nnn


#***************************************************************
#   define parameters
#***************************************************************
# ndim    : dimension; currently only support 3D
# nprops  : no. of input material properties used in uel code (min 1)
# nsvars  : no. of sol. dpdnt var. used and ouput by uel code (min 1)
# uelcode : code for the user element fBrickLam, used in uel input file
# uellinecount  : max no. of entries in a data line in uel input file (16)
# uellinelength : line length in writing uel input file (max. =80 in *Element)
# fnmlinelength : line length in writing fnm input file (max. =132)
ndim    = 3
nprops  = 1
nsvars  = 1
uelcode = 308
uellinecount  = 14
uellinelength = 70
fnmlinelength = 120

#***************************************************************
#   fetch original input file and open uel files
#***************************************************************
# ask for job name
jobname      = raw_input('abaqus job name:')
# abaqus original input file
abqinputfile = jobname+'.inp'    
# uel input file     
uelinputfile = 'uel_'+jobname+'.inp'  
# uel nodes and elems files, to be included in uel input file
uelnodesfile = 'uel_nodes.inp'
uelelemsfile = 'uel_elems.inp'
# open files
abq_input = open(abqinputfile,'r')
uel_input = open(uelinputfile,'w')
uel_nodes = open(uelnodesfile,'w')
uel_elems = open(uelelemsfile,'w')

#***************************************************************
#   Open fnm input files to be written during pre-processing
#***************************************************************
fnm_nodes = open('fnm_nodes.txt','w')  # list of all nodes
fnm_edges = open('fnm_edges.txt','w')  # list of all nodes
fnm_elems = open('fnm_elems.txt','w')  # list of all elems
fnm_predelam = open('fnm_predelam.txt','w') # list of all predelam elems and the predelam interf
fnm_matrix_crack = open('fnm_matrix_crack.txt','w') # matrix crack info



#***************************************************************
#       Read Abaqus input file
#***************************************************************

# Store all lines first
All_lines = abq_input.readlines()

# Find the length of all lines
lenAll = len(All_lines)


#==================================================
# READ HEADER SECTION:
#==================================================
header  = []    # list of header lines
header.extend(All_lines[0:5])


#==================================================
# READ PARTS SECTION: 
# stores all lines of ordinary parts
# reads info of the fnm part and preprocess the following:
# - store part's real nodes, elems & its real nodes
# - find all breakable edges & create fl. nodes for each edge 
# - add fl. nodes to part's node list & elem's node list
# - add edges to part's edge list & elem's edge list
# - store nsets and the associated edges
# - store elsets
# ** NOTE **: 
# - only ONE fnm part is read
# - only a SINGLE-layer mesh is supported for the fnm part
# - only ONE type of elem is supported in the fnm part (all elems will be FNM-ized)
# - only *Node, *Element, *Nset, *Elset sections are supported in the fnm part
#==================================================
fnmparts = []    # list of all fnm parts in the model
ordparts = []    # list of all ordinary parts in the model
predelam = []    # predelam elset (only one) in the model

# find the line no. of *Part and store in jparts
jparts = [j for j,line in enumerate(All_lines) if '*Part' in line]

for jp in jparts:

    # read Part name
    pname = All_lines[jp][12:].rstrip()
    
    # if fnm is not a keyword, then this part is an ordinary part, store everything
    if not ('fnm' in pname):
        ordparts.append( opart(lines=[]) )
        for j in range(jp,lenAll):
            jline = All_lines[j]
            ordparts[-1].lines.append(jline)
            # break out of for loop if end of part is reached
            if('*End Part' in jline):
                break        
    
    # proceed to preprocessing if it is a fnm part
    else:
        # check if fnmparts is not more than one
        if (len(fnmparts) > 1):
            print("ERROR: more than one fnm part is not yet supported!")
            sys.exit()
            
        # ask for layup from user for this fpart
        # isSymLam: yes if the laminate layup is symmetric 
        # rawlayup: layup of the laminate in a list of all plies
        # blklayup: layup of the laminate in a list of plyblocks
        rawlayup = []
        blklayup = []

        # ask if the lamiante layup is symmetric
        print (" Is the model a half-laminate with symmetric mid-plane:  1=yes  2=no ")
        isSymLam = GetInteger2()

        # ask for a list of ply angles
        rawlayup = \
        input('Enter fibre angles (int/float) of all plies in the model, \
        from bottom ply to top or mid ply (if the model is half-laminate); \
        separate each ply by a comma; in case of a single ply, end with a comma:')
        # check if it is a list and if all elements in the list are numbers
        while ( (not isinstance(rawlayup, (list,tuple))) or \
         (not all(isinstance(item, (int,float)) for item in rawlayup)) ):
            rawlayup = \
            input('Enter fibre angles (int/float) of all plies in the model, \
            from bottom ply to top or mid ply (if the model is half-laminate); \
            separate each ply by a comma; in case of a single ply, end with a comma:')  

        # ask for the thickness of a single ply
        plythick = \
        input('Enter the thickness of a single ply (positive real number):')
        while ( not ( isinstance(plythick, float) and plythick > 0 ) ):
            plythick = \
            input('Enter the thickness of a single ply (positive real number):')
           
        # find blocked plies and update blklayup
        # initiate blklayup
        blklayup.append( plyblk(angle=rawlayup[0], nplies=0, thickness=0.0) )
        for plyangle in rawlayup:
            if (plyangle == blklayup[-1].angle):
                blklayup[-1].nplies += 1
                blklayup[-1].thickness += plythick
            else:
                blklayup.append( plyblk(angle=plyangle, nplies=1, thickness=plythick) )

        # if the laminate is symmetric, change the nplies for the last plyblk
        if (isSymLam == 1):
            blklayup[-1].nplies = 2 * blklayup[-1].nplies
    
        # create a new part in parts list
        fnmparts.append(fpart(name=pname, layup=blklayup, nodes=[], NtN=[], \
        edges=[], elems=[], toprnds=[], botrnds=[], nsets=[], elsets=[]))
        
        # find the line of *Node, *Element, *Nset, *Elset
        # ** NOTE: 
        # only ONE Node section (default)
        # only ONE type of elem is supported; all elems in this part will be FNM-ized
        # nset can be multiple
        # only ONE elset, which is the predelam elset, is supported
        jpend  = next( j for j in range(jp,lenAll) if '*End Part'in All_lines[j] )
        jnode  = next( j for j in range(jp,jpend)  if '*Node'    in All_lines[j] )
        jelems =     [ j for j in range(jp,jpend)  if '*Element' in All_lines[j] ]
        jnsets =     [ j for j in range(jp,jpend)  if '*Nset'    in All_lines[j] and 'internal' not in All_lines[j] ]
        jlsets =     [ j for j in range(jp,jpend)  if '*Elset'   in All_lines[j] and 'internal' not in All_lines[j] ]
        
        # report error if multiple elem defs exist in this fnm part
        if (len(jelems) != 1):
            print("ERROR: exactly ONE type of elem is supported in fnm part!")
            sys.exit()
        
        # report error if multiple elset defs exist in this fnm part
        if (len(jlsets) > 1):
            print("ERROR: exactly ONE elset, predelam, is supported in fnm part!")
            sys.exit()
        
        # report error if the elset is not named predelam
        if (len(jlsets) == 1):
            if not ('predelam' in All_lines[jlsets[0]]):
                print("ERROR: only the elset named predelam is supported in fnm part!")
                sys.exit()
        
        # read real(original) nodes of this part 
        for jn in range(jnode+1,jpend):
            nline = All_lines[jn]
            # break out of for loop if end of node section is reached
            if('*' in nline):
                break
            # read the coords of this node into a list of float numbers
            coords = []
            for t in nline.split(','):
                try:
                    coords.append(float(t))
                except ValueError:
                    pass
            # store the coords in nodes list of this part
            fnmparts[-1].nodes.append(node(x=coords[1], y=coords[2], z=coords[3]))
            # update the node-to-node (NtN) matrix row length
            fnmparts[-1].NtN.append([])
        ##check node correctness
        #for nd in fnmparts[-1].nodes:
        #    print(str(nd.x)+','+str(nd.y)+','+str(nd.z))
        
        # append the NtN matrix column to the correct length (nnode)
        for r in range(len(fnmparts[-1].NtN)):
            fnmparts[-1].NtN[r] = [0]*len(fnmparts[-1].NtN)
            
        # read elems of this part
        jelem = jelems[0]
        for je in range(jelem+1,jpend):
            eline = All_lines[je]
            # break out of for loop if end of elem section is reached
            if('*' in eline):
                break
            # read the index and real nodes of this elem into a list of int numbers
            el = []
            for t in eline.split(','):
                try:
                    el.append(int(t))
                except ValueError:
                    pass
            id  = el[0]  # elem index
            nds = el[1:] # elem real nodes
            # store the index and real nodes in elem list of this part
            fnmparts[-1].elems.append(element(index=id, nodes=nds, edges=[]))
            # update the top and bot surf real nodes lists
            halfnds = len(nds)/2
            for j in range(halfnds):
                # bot surf real nodes
                if not(nds[j] in fnmparts[-1].botrnds):
                    fnmparts[-1].botrnds.append(nds[j])
                # top surf real nodes
                if not(nds[halfnds+j] in fnmparts[-1].toprnds):
                    fnmparts[-1].toprnds.append(nds[halfnds+j])
            # form edges
            # in 3D FNM for composites, only edges parallel to shell plane are breakable
            # j = 0: bottom edges; 
            # j = 1: top edges
            for j in range(2):
                # j = 0: loop over real nodes on bot surf
                # j = 1: loop over real nodes on top surf
                for i in range(halfnds):
                    # ith node on surf j
                    row = nds[ j*halfnds + i ] - 1
                    # (i+1)th node on surf j; 
                    # if i is the last node, i+1 is the first node
                    if i == halfnds-1:
                        col = nds[ j*halfnds ] - 1
                    else:
                        col = nds[ j*halfnds + i + 1 ] - 1
                    # fill in the NtN matrix:
                    # fill in the edge index composed of the two end nodes
                    if fnmparts[-1].NtN[row][col]==0:
                    # this pair of nodes hasn't formed an edge
                    # a new edge will be formed
                        # indices of 2 fl. nodes on this new edge
                        # they are two new nodes to be added in the part's node list
                        fn1 = len(fnmparts[-1].nodes)
                        fn2 = fn1 + 1
                        fnmparts[-1].nodes.append( node(0.0,0.0,0.0) )
                        fnmparts[-1].nodes.append( node(0.0,0.0,0.0) )
                        # form a new edge and append to existing list of edges 
                        # the new edge has 4 nodes: 2 real, 2 floating
                        fnmparts[-1].edges.append(edge(nodes=[row+1,col+1,fn1+1,fn2+1]))
                        # fill the new edge index in the NtN matrix
                        # nodes in rev. order makes the same edge in rev. dir.
                        fnmparts[-1].NtN[row][col] =  len(fnmparts[-1].edges)
                        fnmparts[-1].NtN[col][row] = -len(fnmparts[-1].edges)
                        # append this edge no. in this elem
                        fnmparts[-1].elems[-1].edges.append(fnmparts[-1].NtN[row][col]) 
                        # append the fl. nodes in this elem
                        fnmparts[-1].elems[-1].nodes.extend([fn1+1,fn2+1])
                    else:
                    # this pair of nodes has already formed an edge
                        eg = fnmparts[-1].NtN[row][col]
                        # append this edge no. in this elem 
                        fnmparts[-1].elems[-1].edges.append(eg)
                        # find the two fl. nodes on this edge
                        fn1 = fnmparts[-1].edges[abs(eg)-1].nodes[2]
                        fn2 = fnmparts[-1].edges[abs(eg)-1].nodes[3]
                        # append the fl. nodes in this elem
                        # edge is in the same order as saved
                        if eg > 0:
                            fnmparts[-1].elems[-1].nodes.extend([fn1,fn2])
                        # edge is in the rev. order as saved
                        else:
                            fnmparts[-1].elems[-1].nodes.extend([fn2,fn1])
        
        ##check elem correctness
        #for el in fnmparts[-1].elems:
        #    print(str(el.index)+','+str(el.nodes)+','+str(el.edges))
        ##check NtN
        #print(str(fnmparts[-1].NtN))
        
        # read nsets of this part
        for jns in jnsets:
            nline = All_lines[jns].rstrip()
            # remove 'generate' in the line if present
            if ('generate' in All_lines[jns]):
                nline = All_lines[jns].replace(', generate','').rstrip()
            # add this nset in the list of nsets in this fpart
            fnmparts[-1].nsets.append( nset( name=nline, rnodes=[], edges=[] ) )
            # read nodes in the nset
            # if generate is used, then calculate all nodes;
            # otherwise, read all nodes directly
            if ('generate' in All_lines[jns]):
                nline = All_lines[jns+1]
                nl = []
                for t in nline.split(','):
                    try:
                        nl.append(int(t))
                    except ValueError:
                        pass
                nds = nl[0] # start node
                ndf = nl[1] # final node
                try:
                    itv = nl[2] # interval
                except IndexError:
                    itv = 1
                for n in range(nds,ndf+1,itv):
                    fnmparts[-1].nsets[-1].rnodes.append(n)
            else:
                # read the lines of nodes in this nset
                nl = [] # list of node to be filled
                for n in range(jns+1,jpend):
                    nline = All_lines[n]
                    # break out of loop if end of section encountered
                    if ('*' in nline):
                        break
                    for t in nline.split(','):
                        try:
                            nl.append(int(t))
                        except ValueError:
                            pass
                fnmparts[-1].nsets[-1].rnodes.extend(nl)
            # find the edges involved in this nset, 
            # and include the fl. nodes in the nset
            # extract this nset from fpart
            nst = fnmparts[-1].nsets[-1]
            # loop over all node pairs in this nset
            for n1 in range(len(nst.rnodes)-1):
                for n2 in range(n1+1,len(nst.rnodes)):
                    rnd = nst.rnodes[n1]-1
                    cnd = nst.rnodes[n2]-1
                    # if this node pair forms an edge
                    if (fnmparts[-1].NtN[rnd][cnd]!=0):
                        # get this edge number
                        eg = abs(fnmparts[-1].NtN[rnd][cnd])
                        #print(' node '+str(rnd)+' node '+str(cnd)+' forms edge '+str(eg))
                        # store this edge in the nset
                        nst.edges.append(eg)
            # update this nset in fpart
            fnmparts[-1].nsets[-1] = nst 
        
        # read the predelam elset of the fnm part
        for jels in jlsets:
            eline = All_lines[jels].rstrip()
            # remove 'generate' in the line if present
            if ('generate' in All_lines[jels]):
                eline = All_lines[jels].replace(', generate','').rstrip()
            # add this elset in the list of elsets in this fpart
            fnmparts[-1].elsets.append( elset( name=eline, elems=[] ) )
            # read elems in the elset
            # if generate is used, then calculate all elems;
            # otherwise, read all elems directly
            if ('generate' in All_lines[jels]):
                eline = All_lines[jels+1]
                el = []
                for t in eline.split(','):
                    try:
                        el.append(int(t))
                    except ValueError:
                        pass
                els = el[0] # start elem
                elf = el[1] # final elem
                try:
                    itv = el[2] # interval
                except IndexError:
                    itv = 1
                for e in range(els,elf+1,itv):
                    fnmparts[-1].elsets[-1].elems.append(e)
            else:
                # read the lines of nodes in this nset
                el = [] # list of elems to be filled
                for e in range(jels+1,jpend):
                    eline = All_lines[e]
                    # break out of loop if end of section encountered
                    if ('*' in eline):
                        break
                    for t in eline.split(','):
                        try:
                            el.append(int(t))
                        except ValueError:
                            pass
                fnmparts[-1].elsets[-1].elems.extend(el)
            # store the elset in the predelam list
            predelam.append(fnmparts[-1].elsets[-1])

# check if fnmparts is one
if (len(fnmparts) != 1):
    print("ERROR: exactly one fnm part is supported!")
    sys.exit()        

#==================================================
# read assembly section:
# - only a single assembly is supported
#==================================================
assembly = []

# find the line no. of *Assembly and store in jassemblies
jassemblies = [j for j,line in enumerate(All_lines) if '*Assembly' in line]
if (len(jassemblies)!=1):
    print("ERROR: exactly ONE assembly is supported!")
    sys.exit()

# copy everything of the assembly
ja = jassemblies[0]
for j in range(ja,lenAll):
    jline = All_lines[j]
    assembly.append(jline)
    if '*End Assembly' in jline:
        break
#print(assembly)

#==================================================
# read material section:
#==================================================
materials = []

# find the line no. of *Material
jmaterials = [j for j,line in enumerate(All_lines) if '*Material' in line]

# copy everything of the material after the 1st line starting with '*Material'
for ij,jmat in enumerate(jmaterials):
    # break out if it goes to the 2nd *Material line.
    if ij == 1:
        break
    for j in range(jmat,lenAll):
        jline = All_lines[j]
        materials.append(jline)
        if '**' in jline:
            break

#==================================================
# read interaction property section:
#==================================================
interaction_props = []

# find the line no. of *Surface Interaction and store
jinteraction_props = [j for j,line in enumerate(All_lines) if '*Surface Interaction' in line]

# copy everything of the section
for jintp in jinteraction_props:
    for j in range(jintp,lenAll):
        jline = All_lines[j]
        interaction_props.append(jline)
        if '**' in jline:
            break

#==================================================
# read initial boundary conditions and interaction:
#==================================================
bcds  = []       # list of all initial boundary conditions
interaction = [] # the initial interaction

# find the first separation line '** -------------' which separates initial step with the 
# other steps
jdash = next(j for j,line in enumerate(All_lines) \
if '** ----------------------------------------------------------------' in line)

# find the lines with *boundary
jbcds = [j for j in range(0,jdash) if '*Boundary' in All_lines[j]]

# find the line with ** Interaction:
jinteraction = [j for j, line in enumerate(All_lines) if '** Interaction:' in line]

# loop over all bcds, store them without modification
for jb in jbcds:
    for k in range(jb,jdash):
        bline = All_lines[k]
        if ('**' in bline):
            break
        bcds.append(bline)
        
# store interaction without modification
if (len(jinteraction) > 0):
    for k in range(jinteraction[0]+1, jdash):
        iline = All_lines[k]
        if ('**' in iline):
            break
        interaction.append(iline)
    
#print(bcds)


#==================================================
# read steps and
# add control parameters
#==================================================
steps  = []
jsteps = [j for j, line in enumerate(All_lines) if '*Step' in line]
for jstep in jsteps:
    for j in range(jstep,lenAll):
        jline = All_lines[j]
        if ('*End Step' in jline):
            # add control parameters
            steps.append('*Controls, reset\n')
            steps.append('*Controls, parameters=time incrementation\n')
            steps.append('3200, 4000, , 6000, 4800, 50, , 125, , , \n')
            steps.append('**\n')
            steps.append(jline)
            break
        steps.append(jline)
#print(step)

   


#***************************************************************
#       write nodes of the fnm part
#***************************************************************  
# extract layup of the fnm part
blklayup = fnmparts[0].layup
# find no. of plyblocks
nplyblk = len(blklayup)
# find no. of nodes and edges in a ply of this mesh
nnode_p = len(fnmparts[0].nodes)
nedge_p = len(fnmparts[0].edges)
# find internal nodes for an interface of this mesh
nnodein = nedge_p
# find the total no. of nodes in this mesh
nnodett = nplyblk * nnode_p + (nplyblk-1) * nnodein

# write fnm_nodes 1st line: nnode
fnm_nodes.write(str(nnodett)+' \n')

for ipb in range(nplyblk):
    # calculate the bot and top real node z-coordinate
    # zbot = 0 if this is the 1st plyblk
    if (ipb == 0):
        zbot = 0.0
    # zbot = thickness of last plyblk
    else:
        zbot = zbot + blklayup[ipb-1].thickness
    # ztop = zbot + thickness of this plyblk
    ztop = zbot + blklayup[ipb].thickness
    # loop over all nodes in the single ply mesh
    for cntr0, nd in enumerate(fnmparts[0].nodes):
        # current node id of the node on the ith plyblk
        cntr = ipb * nnode_p + (cntr0+1)
        # check if this node is a real node on the bot/top surf, or a fl. node
        # must use cntr0+1 as bot/topnds lists are for nodes in 1st plyblk
        if ((cntr0+1) in fnmparts[0].botrnds):
            zz = zbot
        elif ((cntr0+1) in fnmparts[0].toprnds):
            zz = ztop
        else:
            zz = 0.0
        # write this node coords in uel_nodes.inp
        uel_nodes.write\
        (str(cntr)+', '+str(nd.x)+', '+str(nd.y)+', '+str(zz)+'\n')
        # write this node coords in fnm node_list
        fnm_nodes.write\
        (str(nd.x)+' '+str(nd.y)+' '+str(zz)+' \n')

# write the additional nodes of interfaces if they're present
if (nplyblk > 1):
    for jintf in range(nplyblk-1):
        # loop over all edges in the base ply mesh
        # each edge has one additional node 
        for cntr0 in range(nedge_p):
            cntr = nplyblk * nnode_p + jintf * nedge_p + cntr0 + 1
            # write this node coords in uel_nodes.inp
            uel_nodes.write(str(cntr)+', 0.0, 0.0, 0.0\n')
            # write this node coords in fnm node_list
            fnm_nodes.write('0.0 0.0 0.0 \n')
            

#***************************************************************
#       write edges
#*************************************************************** 
# find the total no. of edges in this mesh
nedgett = nplyblk * nedge_p  
fnm_edges.write(str(nedgett)+' \n')

# update bcd edges
for nst in fnmparts[0].nsets:
    # only constrain the edges in nst with keyword 'tie'
    if ('tie' in nst.name):
        # if all the real nodes in this nst are on the bot surface, then
        # only store the bot plyblk nodes, DO NOT include the other plies
        if all(n in fnmparts[0].botrnds for n in nst.rnodes):
            pstart = 0
            pend   = 1
        # if all the real nodes in this nst are on the top surface, then
        # only store the top plyblk nodes, DO NOT include the other plies
        elif all(n in fnmparts[0].toprnds for n in nst.rnodes):
            pstart = nplyblk-1
            pend   = nplyblk
        # otherwise, store corresponding nodes of all plyblks
        else:
            pstart = 0
            pend   = nplyblk
        # find the bcd edges
        for jpb in range(pstart,pend):
            for jeg in nst.edges:
                jedge = jeg + jpb * nedge_p 
                fnm_edges.write(str(jedge)+' \n')
# mark end of file with -1
fnm_edges.write('-1 \n')

#***************************************************************
#       write elems
#***************************************************************   
# find no. of elems in a ply mesh
nelem_p = len(fnmparts[0].elems) 
# find total no. of elems in the laminate
# it is the same as nelem_p, as a fnm elem contains all plies&interfs
nelemtt = nelem_p
# find the no. of r+f nodes in an elem of a single plyblk
elnndrf_p = len(fnmparts[0].elems[0].nodes)
# find the no. of edges in an elem of a single plyblk
elnedge_p = len(fnmparts[0].elems[0].edges)
# find the no. of r+f nodes in an elem of the laminate
elnndrf_l = elnndrf_p * nplyblk
# find the no. of interface internal nodes in an elem of the laminate
elnndin_l = elnedge_p * (nplyblk-1)
# find the total no. of nodes in an elem of the laminate
elnndtt_l = elnndrf_l + elnndin_l
# find the no. of edges in an elem of the laminate
elnedge_l = elnedge_p * nplyblk

fnm_elems.write(str(nelemtt)+' '+str(elnndtt_l)+' '+str(elnedge_l)+' \n')

fnm_elems.write(str(nplyblk)+' \n')
# write layup array
for jpb in range(nplyblk):
    angle  = str(blklayup[jpb].angle)
    nplies = str(blklayup[jpb].nplies)
    fnm_elems.write(angle+' '+nplies+' \n')

for jel in range(nelemtt):
    elnds_p = []
    elegs_p = []
    elnds_l = []
    elegs_l = []
    # find the r+f nodes in this elem ply
    elnds_p = fnmparts[0].elems[jel].nodes
    # find the edges in this elem ply
    elegs_p = fnmparts[0].elems[jel].edges
    # find the r+f nodes & edges in this elem laminate; 
    # abs is needed to remove the negative sign on some of the edge no.
    for jpb in range(nplyblk):
        elnds_l.extend( [     x  + nnode_p * jpb  for x in elnds_p ] )
        elegs_l.extend( [ abs(x) + nedge_p * jpb  for x in elegs_p ] )
    # find the internal nodes of interfs in this elem laminate
    # intern. nodes are listed after r and f nodes in the node list
    # they are assigned to each edges of the interfaces
    # so the elem's edge connec is used for assignment of intern. nodes
    if nplyblk > 1:
        for jit in range(nplyblk-1):
            elnds_l.extend( [ abs(x) + nnode_p * nplyblk + nedge_p * jit for x in elegs_p ] )
            
    #**** write elem's nodal and edge connec to uel&fnm_elems ****
    
    #** node cnc 
    # start the line with elem index jel+1
    eline = [str(jel+1)+',']  # node cnc dataline for uel_elems
    fline = ['']              # node cnc dataline for fnm_elems
    cntr  = 0                 # line entry count  for uel_elems
    # add the node no. to the line one by one
    for k in elnds_l:
        # if the uel line gets too long, continue on next line
        if (len(eline[-1]+str(k)) >= uellinelength) or \
        (cntr >= uellinecount):
            eline.append('')
            cntr = 0
        # add the node no. to the line and update line count
        eline[-1] = eline[-1]+str(k)+','
        cntr      = cntr + 1
        # add the node no. to the fnm line 
        fline[0] = fline[0]+str(k)+' '
    # remove the last comma from the eline
    eline[-1] = eline[-1][:-1]
    
    #** edge cnc
    gline = ['']  # edge cnc dataline for fnm_elems
    # add the node no. to the line one by one
    for k in elegs_l:
        # add the node no. to the line
        gline[0] = gline[0]+str(k)+' '
    
    # write the line of elem node connec
    for l in eline:
        uel_elems.write(l+'\n')
        
    # write this elem in fnm_elems
    # write nodecnc array
    fnm_elems.write(fline[0]+' \n')
    # write edgecnc array
    fnm_elems.write(gline[0]+' \n')



#***************************************************************
#       write predelam
#*************************************************************** 
# if no predelam, write 0 in fnm_predelam
if (len(predelam) == 0):
    fnm_predelam.write('0 \n')
    
# check if there is only ONE predelam
if (len(predelam) > 1):
    print("ERROR: more than one predelam is not yet supported!")
    sys.exit()

# write the elem indices in the predelam elset
for pd in predelam:
    npdelem = len(pd.elems)
    fnm_predelam.write(str(npdelem)+' \n')
    # ask for the predelam interface no.
    pdinterf = \
    input('Enter the pre-delamination interface, \
    1 means the first interface from the bottom (positive integer number):')
    while ( not ( isinstance(pdinterf, int) and pdinterf > 0 ) ):
        pdinterf = \
        input('Enter the pre-delamination interface, \
        1 means the first interface from the bottom (positive integer number):')
    fnm_predelam.write(str(pdinterf)+' \n')
    # write all elems in the predelam elset
    for jel in pd.elems:
        fnm_predelam.write(str(jel)+' \n')


#***************************************************************
#       write matrix_crack
#*************************************************************** 
## ask for minimum crack spacing
#minspacing = \
#input('Enter the minimum spacing between two matrix cracks (positive real number):')
#while ( not ( isinstance(minspacing, float) and minspacing > 0.0 ) ):
#    minspacing = \
#    input('Enter the minimum spacing between two matrix cracks (positive real number):')

# ask for min and max elem sizes
minl = 0
while ( not ( isinstance(minl, float) and minl > 0.0 ) ):
    minl = \
    input('Enter the minimum element size in the mesh (positive real number):')
maxl = 0
while ( not ( isinstance(maxl, float) and maxl >= minl ) ):
    maxl = \
    input('Enter the maximum element size in the mesh (real number, >= min_elem_size):')
    
maxncrack = int(math.sqrt(nelem_p))

#fnm_matrix_crack.write(str(minspacing)+' \n')
fnm_matrix_crack.write(str(minl)+' '+str(maxl)+' \n')
fnm_matrix_crack.write(str(nplyblk)+' \n')
fnm_matrix_crack.write(str(maxncrack)+' \n')


#***************************************************************
#       write uel input file
#*************************************************************** 
#**** write HEADER ****
for hline in header:
    uel_input.write(str(hline[0:]))

#**** write ordinary part ****
for op in ordparts:
    for line in op.lines:
        uel_input.write(str(line[0:]))

#**** write FNM PART   ****
# part name 
uel_input.write('*Part, name='+fnmparts[0].name+'\n')
# part nodes
uel_input.write('*NODE,INPUT=uel_nodes.inp    \n')
# user element definition
uel_input.write('*USER ELEMENT, TYPE=U'+str(uelcode)+\
', NODES='+str(elnndtt_l)+', COORDINATES='+str(ndim)+\
', PROPERTIES='+str(nprops)+', VARIABLES='+str(nsvars)+'\n')
uel_input.write('1,2,3\n')
# elements and fnm elset
uel_input.write('*ELEMENT, TYPE=U'+str(uelcode)+', ELSET=fnm, INPUT=uel_elems.inp \n')
# write the mandatory uel property line (not needed for calculation)
uel_input.write('*UEL PROPERTY, ELSET=fnm\n')
uel_input.write('1\n')
# elset (predelam) does not need to be written here in uel input    
# nsets
nsets  = fnmparts[0].nsets
for nst in nsets:
    # write the nset name
    uel_input.write(nst.name+'\n')
    # nst dataline for uel_input, to be filled, and line count initiated
    nstline = ['']
    cntr    = 0
    # if all the real nodes in this nst are on the bot surface, then
    # only store the bot plyblk nodes, DO NOT include the other plies
    if all(n in fnmparts[0].botrnds for n in nst.rnodes):
        pstart = 0
        pend   = 1
    # if all the real nodes in this nst are on the top surface, then
    # only store the top plyblk nodes, DO NOT include the other plies
    elif all(n in fnmparts[0].toprnds for n in nst.rnodes):
        pstart = nplyblk-1
        pend   = nplyblk
    # otherwise, store corresponding nodes of all plyblks
    else:
        pstart = 0
        pend   = nplyblk
    # find nst nodes in plyblks
    for jpb in range(pstart,pend):
        # add the real nodes to the list one by one
        for n in nst.rnodes:
            # find the corresponding node on the jpb-th plyblk
            k = n + jpb * nnode_p
            # if the uel line gets too long, continue on the next line
            if (len(nstline[-1]+str(k)) >= uellinelength) or \
               (cntr >= uellinecount):
                nstline.append('')
                cntr = 0
            # add the node no. to the line and update line count
            nstline[-1] = nstline[-1]+str(k)+','
            cntr        = cntr + 1
        # add the fl. nodes to the list one by one, if it is not a tie nst
        if not ('tie' in nst.name):
            for eg in nst.edges:
                k1 = fnmparts[0].edges[eg-1].nodes[2] + jpb * nnode_p
                k2 = fnmparts[0].edges[eg-1].nodes[3] + jpb * nnode_p
                # if the uel line gets too long, continue on the next line
                if (len(nstline[-1]+str(k1)+str(k2)) >= uellinelength) or \
                   (cntr >= uellinecount):
                    nstline.append('')
                    cntr = 0
                # add the nodes to the line
                nstline[-1] = nstline[-1]+str(k1)+','+str(k2)+','
                cntr        = cntr + 2
    # remove the last comma from the line
    nstline[-1] = nstline[-1][:-1]
    # write all original nodes of the nset
    for nstl in nstline:
        uel_input.write(nstl+'\n')
# end writing fnm part in uel input
uel_input.write('*End Part\n')

#**** write ASSEMBLY ****
for aline in assembly:
    uel_input.write(str(aline[0:]))
    
#**** write Material ****
for mline in materials:
    uel_input.write(str(mline[0:]))
    
#**** write interaction property ****
for iline in interaction_props:
    uel_input.write(str(iline[0:]))

#**** write initial BCDs ****
for bline in bcds:
    uel_input.write(str(bline[0:]))
    
#**** write initial interaction ****
for iline in interaction:
    uel_input.write(str(iline[0:]))

#**** write steps ****
for sline in steps:
    uel_input.write(str(sline[0:]))



#*************************************************************** 
# close all open files
#*************************************************************** 
#   close input files
abq_input.close()
uel_input.close()
#   close nodes files
fnm_nodes.close()
uel_nodes.close()
#   close edges file
fnm_edges.close()
#   close elems files
fnm_elems.close()
uel_elems.close()
#   close predelam file
fnm_predelam.close()
#   close matrix crack file
fnm_matrix_crack.close()


#*************************************************************** 
# copy fnm input file to main directory
#*************************************************************** 

# get current working directory
cwd = os.getcwd()
# parent directory of cwd
pwd = os.path.dirname(cwd)
# copy fnm input file to parent directory of preprocessing directory (which is assumed to be the working directory)
shutil.copy (uelinputfile,pwd)
shutil.copy (uelnodesfile,pwd)
shutil.copy (uelelemsfile,pwd)
