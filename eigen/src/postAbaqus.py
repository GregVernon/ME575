import os
import sys

def main():
	f = evaluateObjectiveFunction()
	writeObjectiveValue(f)
	
def evaluateObjectiveFunction():
	odb = session.openOdb('C:/Users/gregj/Documents/GitHub/ME575/eigen/src/drumEigen.odb')
	odb = session.odbs['C:/Users/gregj/Documents/GitHub/ME575/eigen/src/drumEigen.odb']
	f = odb.steps['Step-1'].frames[1].frequency
	odb.close()
	return f
	
def writeObjectiveValue(f):
	fRes = open("C:/Users/gregj/Documents/Abaqus/Temp/objectiveFunction.csv","w+")
	fRes.write(str(f))

if __name__ == "__main__":
	f = main()