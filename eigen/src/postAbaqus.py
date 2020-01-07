def main():
	session.openOdb('C:/Users/gregj/Documents/Abaqus/Temp/drumEigen.odb')
	odb = session.odbs['C:/Users/gregj/Documents/Abaqus/Temp/drumEigen.odb']
	f = odb.steps['Step-1'].frames[1].frequency
	odb.close()
	return f

if __name__ == "__main__"
	f = main()