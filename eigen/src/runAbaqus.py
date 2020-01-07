# -*- coding: mbcs -*-
import os
import sys

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

def main(x,y):
	makeGeom(x,y)
	makeMesh()
	makeMaterial()
	makeAssembly()
	makeAssembly()
	makeStep()
	makeBoundaryConditions()
	makeJob()
	submitJob()

def makeGeom(x,y):
	mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=10.0)
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0), point1=(1.0, 0.0))
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])
	del mdb.models['Model-1'].sketches['__profile__']
	
	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.14, name='__profile__', 
		sheetSize=5.64, transform=
		mdb.models['Model-1'].parts['Part-1'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['Part-1'].faces[0], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['Part-1'].edges[0], 
		sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
	
	mdb.models['Model-1'].parts['Part-1'].projectReferencesOntoSketch(filter=
		COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(x, y+0.01), point2=(x, y))
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(x, y), point2=(x+0.01, y))
	
	mdb.models['Model-1'].parts['Part-1'].PartitionFaceBySketch(faces=
		mdb.models['Model-1'].parts['Part-1'].faces[0:1], sketch=
		mdb.models['Model-1'].sketches['__profile__'], sketchUpEdge=
		mdb.models['Model-1'].parts['Part-1'].edges[0])
	
	del mdb.models['Model-1'].sketches['__profile__']
	mdb.models['Model-1'].parts['Part-1'].ignoreEntity(entities=(
		mdb.models['Model-1'].parts['Part-1'].vertices[1:3], 
		mdb.models['Model-1'].parts['Part-1'].edges[0:2]))

def makeMesh():
	mdb.models['Model-1'].parts['Part-1'].setMeshControls(elemShape=QUAD, regions=mdb.models['Model-1'].parts['Part-1'].faces[0:1])
	mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.14)
	mdb.models['Model-1'].parts['Part-1'].generateMesh()
	
	mdb.models['Model-1'].parts['Part-1'].setElementType(elemTypes=(ElemType(
		elemCode=M3D4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
		hourglassControl=DEFAULT), ElemType(elemCode=M3D3, elemLibrary=STANDARD)), 
		regions=(mdb.models['Model-1'].parts['Part-1'].faces[0:1], ))
	
	mdb.models['Model-1'].parts['Part-1'].setElementType(elemTypes=(ElemType(
		elemCode=S8R, elemLibrary=STANDARD), ElemType(elemCode=STRI65, 
		elemLibrary=STANDARD)), regions=(
		mdb.models['Model-1'].parts['Part-1'].faces[0:1], ))

def makeMaterial():
	mdb.models['Model-1'].Material(name='Material-1')
	mdb.models['Model-1'].materials['Material-1'].Density(table=((0.0001, ), ))
	mdb.models['Model-1'].materials['Material-1'].Elastic(table=((100000.0, 0.3), ))
	
	mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
		integrationRule=GAUSS, material='Material-1', name='Section-1', 
		nodalThicknessField='', numIntPts=3, poissonDefinition=DEFAULT, 
		preIntegrate=OFF, temperature=GRADIENT, thickness=0.01, thicknessField='', 
		thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
		
	mdb.models['Model-1'].parts['Part-1'].Set(faces=
		mdb.models['Model-1'].parts['Part-1'].faces[0:1], name='Set-2')
	mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
		offsetField='', offsetType=MIDDLE_SURFACE, region=
		mdb.models['Model-1'].parts['Part-1'].sets['Set-2'], sectionName=
		'Section-1', thicknessAssignment=FROM_SECTION)

def makeAssembly():
	mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', part=mdb.models['Model-1'].parts['Part-1'])

def makeStep():
	mdb.models['Model-1'].FrequencyStep(limitSavedEigenvectorRegion=None, name='Step-1', numEigen=1, previous='Initial')

def makeBoundaryConditions():
	mdb.models['Model-1'].rootAssembly.Set(edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges[0:1], name='Set-1')
	mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
		distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
		'BC-1', region=mdb.models['Model-1'].rootAssembly.sets['Set-1'], u1=0.0, 
		u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
	
	mdb.models['Model-1'].rootAssembly.Set(name='Set-2', vertices=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].vertices[1:2])
	mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
		distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
		'BC-2', region=mdb.models['Model-1'].rootAssembly.sets['Set-2'], u1=0.0, 
		u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

def makeJob():
	mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
		explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
		memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, name=
		'drumEigen', nodalOutputPrecision=SINGLE, queue=None, resultsFormat=ODB, 
		scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

def submitJob():
	mdb.jobs['drumEigen'].submit(consistencyChecking=OFF)

if __name__ == "__main__":
	x,y = sys.argv[-2:]
	main(float(x),float(y))
	