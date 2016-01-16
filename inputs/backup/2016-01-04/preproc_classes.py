#***************************************************************
#****************** Define material classes ********************
#***************************************************************

class lamina_modulus:

    def __init__(self, E1, E2, G12, G23, nu12, nu23):
        self.E1   = E1
        self.E2   = E2
        self.G12  = G12
        self.G23  = G23
        self.nu12 = nu12
        self.nu23 = nu23

class lamina_strength:

    def __init__(self, Xt, Xc, Yt, Yc, Sl, St):
        self.Xt = Xt
        self.Xc = Xc
        self.Yt = Yt
        self.Yc = Yc
        self.Sl = Sl
        self.St = St

class lamina_fibreToughness:

    def __init__(self, GfcT, GfcC):
        self.GfcT = GfcT
        self.GfcC = GfcC

class lamina:

    def __init__(self, modulus, strength, fibreToughness):
        self.modulus  = modulus
        self.strength = strength
        self.fibreToughness  = fibreToughness


class cohesive_modulus:

    def __init__(self, Dnn, Dtt, Dll):
        self.Dnn = Dnn
        self.Dtt = Dtt
        self.Dll = Dll

class cohesive_strength:

    def __init__(self, tau_nc, tau_tc, tau_lc):
        self.tau_nc = tau_nc
        self.tau_tc = tau_tc
        self.tau_lc = tau_lc

class cohesive_toughness:

    def __init__(self, Gnc, Gtc, Glc, alpha):
        self.Gnc = Gnc
        self.Gtc = Gtc
        self.Glc = Glc
        self.alpha = alpha

class cohesive:

    def __init__(self, modulus, strength, toughness):
        self.modulus   = modulus
        self.strength  = strength
        self.toughness = toughness




#***************************************************************
#****************** Define geometry classes ********************
#***************************************************************
class plyblk:

    def __init__(self, angle, nplies, thickness):
        self.angle  = angle
        self.nplies = nplies
        self.thickness = thickness


class node:

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class edge:

    def __init__(self, nodes):
        self.nodes = nodes


class element:

    def __init__(self, index, nodes, edges):
        self.index = index
        self.nodes = nodes
        self.edges = edges

   
class nset:

    def __init__(self, name, rnodes, edges):
        self.name     = name
        self.rnodes   = rnodes
        self.edges    = edges


class elset:

    def __init__(self, name, elems):
        self.name     = name
        self.elems    = elems
    

class opart:

    def __init__(self, lines):
        self.lines = lines
    
        
class fpart:

    def __init__(self, name, layup, nodes, NtN, edges, elems, toprnds, botrnds, nsets, elsets):
        self.name    = name
        self.layup   = layup
        self.nodes   = nodes
        self.NtN     = NtN
        self.edges   = edges
        self.elems   = elems
        self.toprnds = toprnds
        self.botrnds = botrnds
        self.nsets   = nsets
        self.elsets  = elsets


#class bcd:
#
#    def __init__(self, name, type, nsets, firstdof=0, lastdof=0, value=0.):
#        self.name     = name
#        self.type     = type
#        self.nsets    = nsets
#        self.firstdof = firstdof
#        self.lastdof  = lastdof
#        self.value    = value
    
class instance:

    def __init__(self, lines):
        self.lines = lines  
#class instance:
#
#    def __init__(self, name, part, translation=[0.,0.,0.]):
#        self.name = name
#        self.part = part
#        self.translation = translation

    
class assembly:

    def __init__(self, name, instances, nsets, elsets, surfaces, constraints):
        self.name      = name
        self.instances = instances
        self.nsets     = nsets
        self.elsets    = elsets
        self.surfaces  = surfaces
        self.constraints = constraints
        

