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

def main(simID):
	makeGeom(simID)
	makeMesh()
	makeMaterial()
	makeAssembly()
	makeAssembly()
	makeStep()
	makeBoundaryConditions()
	makeJob(simID)
	submitJob()

def makeGeom(simID):
	mdb.ModelFromInputFile(inputFileName='test_' + str(simID) + '.inp', name='temp')
	del mdb.models['Model-1']
	mdb.models.changeKey(fromName='temp', toName='Model-1')
	del mdb.models['Model-1'].steps['DefaultSet']
	# mdb.models['Model-1'].parts.changeKey(fromName='PART-DEFAULT', toName='Part-1')

	mdb.models['Model-1'].parts['PART-DEFAULT'].DatumPlaneByPrincipalPlane(offset=0.0, 
		principalPlane=XYPLANE)
	mdb.models['Model-1'].parts['PART-DEFAULT'].DatumPlaneByPrincipalPlane(offset=1.0, 
		principalPlane=YZPLANE)
	mdb.models['Model-1'].parts['PART-DEFAULT'].DatumAxisByTwoPlane(plane1=
		mdb.models['Model-1'].parts['PART-DEFAULT'].datums[6], plane2=
		mdb.models['Model-1'].parts['PART-DEFAULT'].datums[7])
	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.19, name='__profile__', 
		sheetSize=7.85, transform=
		mdb.models['Model-1'].parts['PART-DEFAULT'].MakeSketchTransform(
		sketchPlane=mdb.models['Model-1'].parts['PART-DEFAULT'].datums[6], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=mdb.models['Model-1'].parts['PART-DEFAULT'].datums[8], 
		sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
	mdb.models['Model-1'].parts['PART-DEFAULT'].projectReferencesOntoSketch(filter=
		COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
	mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
		0.0, 0.0), point1=(1.0, 0.0))
	mdb.models['Model-1'].parts['PART-DEFAULT'].Shell(sketch=
		mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=RIGHT, 
		sketchPlane=mdb.models['Model-1'].parts['PART-DEFAULT'].datums[6], 
		sketchPlaneSide=SIDE1, sketchUpEdge=
		mdb.models['Model-1'].parts['PART-DEFAULT'].datums[8])
	del mdb.models['Model-1'].sketches['__profile__']

def makeMesh():
	mdb.meshEditOptions.setValues(enableUndo=True, maxUndoCacheElements=0.5)
	mdb.models['Model-1'].parts['PART-DEFAULT'].associateMeshWithGeometry(elemFaces=
		mdb.models['Model-1'].parts['PART-DEFAULT'].elementFaces, geometricEntity=
		mdb.models['Model-1'].parts['PART-DEFAULT'].faces[0])

	mdb.models['Model-1'].parts['PART-DEFAULT'].setElementType(elemTypes=(ElemType(
		elemCode=S8R, elemLibrary=STANDARD), ElemType(elemCode=S4R, 
		elemLibrary=STANDARD)), regions=(
		mdb.models['Model-1'].parts['PART-DEFAULT'].faces[0:1], ))

def makeMaterial():
	mdb.models['Model-1'].Material(name='Material-1')
	mdb.models['Model-1'].materials['Material-1'].Density(table=((0.0001, ), ))
	mdb.models['Model-1'].materials['Material-1'].Elastic(table=((100000.0, 0.3), ))
	
	mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
		integrationRule=GAUSS, material='Material-1', name='Section-1', 
		nodalThicknessField='', numIntPts=3, poissonDefinition=DEFAULT, 
		preIntegrate=OFF, temperature=GRADIENT, thickness=0.01, thicknessField='', 
		thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
		
	mdb.models['Model-1'].parts['PART-DEFAULT'].Set(faces=
		mdb.models['Model-1'].parts['PART-DEFAULT'].faces[0:1], name='Set-2')
	mdb.models['Model-1'].parts['PART-DEFAULT'].SectionAssignment(offset=0.0, 
		offsetField='', offsetType=MIDDLE_SURFACE, region=
		mdb.models['Model-1'].parts['PART-DEFAULT'].sets['Set-2'], sectionName=
		'Section-1', thicknessAssignment=FROM_SECTION)

def makeAssembly():
	mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='PART-DEFAULT_1', part=mdb.models['Model-1'].parts['PART-DEFAULT'])

def makeStep():
	mdb.models['Model-1'].FrequencyStep(limitSavedEigenvectorRegion=None, name='Step-1', numEigen=1, previous='Initial')

def makeBoundaryConditions():
	mdb.models['Model-1'].rootAssembly.Set(edges=mdb.models['Model-1'].rootAssembly.instances['PART-DEFAULT_1'].edges[0:1], name='Set-1')
	mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
		distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
		'BC-1', region=mdb.models['Model-1'].rootAssembly.instances['PART-DEFAULT_1'].sets['NS2'], u1=0.0, 
		u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
	
	# mdb.models['Model-1'].rootAssembly.Set(name='Set-2', vertices=mdb.models['Model-1'].rootAssembly.instances['PART-DEFAULT_1'].vertices[1:2])
	mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
		distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
		'BC-2', region=mdb.models['Model-1'].rootAssembly.instances['PART-DEFAULT_1'].sets['NS1'], u1=0.0, 
		u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

def makeJob(simID):
	mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
		explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
		memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, name=
		'drumEigen_' + str(simID), nodalOutputPrecision=SINGLE, queue=None, resultsFormat=ODB, 
		scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

def submitJob():
	mdb.jobs['drumEigen_' + str(simID)].submit(consistencyChecking=OFF)

if __name__ == "__main__":
	simID = sys.argv[-1]
	main(simID)
	