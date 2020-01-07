import os
import sys

# sys.path.append(r"C:\Program Files\Cubit 15.4\bin\")
# import cubit
# cubit.init(['cubit','-nojournal','-nogui','nographics'])

f = open("inputData.csv","r")
fLines = f.readlines()
x = []
y = []
for i in range(0,len(fLines)):
  cLine = fLines[i].strip().split(",")
  x.append(float(cLine[0]))
  y.append(float(cLine[1]))

f.close()

cubit.cmd("reset")
cubit.cmd("create surface circle radius 1 zplane ")
for i in range(0,len(x)):
  cubit.cmd("create vertex " + str(x[i]) + " " + str(y[i]) + " 0 on surface 1 ")
  if i == 0:
    newVertex = "2"
  else:
    newVertex = newVertex + " " + str(i+2)
  
cubit.cmd("imprint volume all with vertex " + newVertex)
cubit.cmd("delete free vertex all")
cubit.cmd("compress ids")

cubit.cmd("block 1 surf 1")
cubit.cmd("nodeset 1 vertex " + newVertex)
cubit.cmd("nodeset 2 curve all")
cubit.cmd("block 1 element type SHELL4")
cubit.cmd("surface all scheme pave")

maxNodes = 1000
validMesh = False
m0 = 1e-2
m2 = 1.0
while validMesh == False:
  cubit.cmd("delete mesh")
  m1 = (m0 + m2) / 2.
  if ((abs(m1 - m0) < 1e-5) or (abs(m1 - m2) < 1e-5)):
    validMesh = True
  cubit.cmd("surf all size " + str(m1))
  cubit.cmd("mesh surface all")
  nNodes = cubit.get_node_count()
  if nNodes > maxNodes:
    m0 = m1
  elif nNodes < maxNodes:
    m2 = m1
  else:
    validMesh = True
  if nNodes > maxNodes:
    validMesh = False

# cubit.cmd("block 1 element type shell4")
cubit.cmd("export abaqus 'test.inp'  overwrite  everything ")