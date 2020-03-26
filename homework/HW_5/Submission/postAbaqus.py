import os
import sys

def main(simID):
	f = evaluateObjectiveFunction(simID)
	writeObjectiveValue(f,simID)
	
def evaluateObjectiveFunction(simID):
	odb = session.openOdb('drumEigen_' + str(simID) +'.odb')
	odb = session.odbs['drumEigen_' + str(simID) + '.odb']
	f = odb.steps['Step-1'].frames[1].frequency
	odb.close()
	return f
	
def writeObjectiveValue(f,simID):
	fRes = open("objectiveFunction_" + str(simID) + ".csv","w+")
	fRes.write(str(f))

if __name__ == "__main__":
	simID = sys.argv[-1]
	f = main(int(simID))