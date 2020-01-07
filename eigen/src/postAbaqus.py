import os
import sys

def main():
	f = evaluateObjectiveFunction()
	writeObjectiveValue(f)
	
def evaluateObjectiveFunction():
	odb = session.openOdb('drumEigen.odb')
	odb = session.odbs['drumEigen.odb']
	f = odb.steps['Step-1'].frames[1].frequency
	odb.close()
	return f
	
def writeObjectiveValue(f):
	fRes = open("objectiveFunction.csv","w+")
	fRes.write(str(f))

if __name__ == "__main__":
	f = main()