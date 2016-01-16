from numpy import*

# set lambda
l = 0.75

# initialize T, Ta, Tb
Tmatrix  = array([[0. for x in range(7)] for x in range(7)])
TmatrixA = array([[0. for x in range(5)] for x in range(5)])
TmatrixB = array([[0. for x in range(5)] for x in range(5)])

# initialize K, Ka, Kb
Kmatrix  = array([[0. for x in range(7)] for x in range(7)])
KmatrixA = array([[0. for x in range(5)] for x in range(5)])
KmatrixB = array([[0. for x in range(5)] for x in range(5)])

# Form T
for i in range(6):
  Tmatrix[i][i] = 1.
Tmatrix[6][1] = 1.-l
Tmatrix[6][4] = l
#print(Tmatrix)

# Form Ta and Tb
for i in range(4):
  TmatrixA[i][i] = 1.
  TmatrixB[i][i] = 1.
TmatrixA[4][1] = 1.-l
TmatrixA[4][2] = l
TmatrixB[4][0] = 1.-l
TmatrixB[4][3] = l
#print(TmatrixA)
#print(TmatrixB)

# Form Ka and Kb
for i in range(5):
  KmatrixA[i][i] = 1.
  KmatrixB[i][i] = 2.
KmatrixA[0][1] = 0.5
KmatrixA[1][0] = 0.5
KmatrixA[0][4] = 0.5
KmatrixA[4][0] = 0.5
KmatrixB[0][1] = 0.5
KmatrixB[1][0] = 0.5
KmatrixB[0][4] = 0.5
KmatrixB[4][0] = 0.5

# connectivities of el A and el B
elA = [1,2,5,4,7]
elB = [2,3,6,5,7]

# Assemble Ka and Kb into K
for i, ni in enumerate(elA):
  for j, nj in enumerate(elA):
    Kmatrix[ni-1][nj-1] = Kmatrix[ni-1][nj-1] + KmatrixA[i][j]  
for i, ni in enumerate(elB):
  for j, nj in enumerate(elB):
    Kmatrix[ni-1][nj-1] = Kmatrix[ni-1][nj-1] + KmatrixB[i][j]
    
# Kmatrix after constrain
#cKmatrix = dot(dot(Tmatrix.transpose(),Kmatrix),Tmatrix)
cKmatrix = dot(Tmatrix.transpose(),dot(Kmatrix,Tmatrix))
print(cKmatrix)

# KmatrixA after constrain
cKmatrixA = dot(dot(TmatrixA.transpose(),KmatrixA),TmatrixA)
cKmatrixB = dot(dot(TmatrixB.transpose(),KmatrixB),TmatrixB)

# initialize cK2
cKmatrix2  = array([[0. for x in range(7)] for x in range(7)])

# Assemble cKa and cKb into cK2
for i, ni in enumerate(elA):
  for j, nj in enumerate(elA):
    cKmatrix2[ni-1][nj-1] = cKmatrix2[ni-1][nj-1] + cKmatrixA[i][j] 
for i, ni in enumerate(elB):
  for j, nj in enumerate(elB):
    cKmatrix2[ni-1][nj-1] = cKmatrix2[ni-1][nj-1] + cKmatrixB[i][j]
print(cKmatrix2)

#print(Kmatrix.transpose())
#print(TmatrixA)
#print(TmatrixB)