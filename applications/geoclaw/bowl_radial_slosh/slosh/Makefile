PYTHON=python
EXE=slosh

topo:
	${PYTHON} maketopo.py

data:
	${PYTHON} setrun.py

output:
	@echo "Usage : "
	@echo "To run example in serial : "
	@echo "       $$ ${EXE}"
	@echo " "
	@echo "To run example in parallel : "
	@echo "       $$ mpirun -n 4 ${EXE}"

plots: 
	${PYTHON} setplot.py

clean:
	rm -f *.data pyclaw.log gauges.kml bowl.topotype2 hump.xyz

usage:
	@echo "Usage : "
	@echo "To create data files needed by GeoClaw : "
	@echo "       $$ make data"
	@echo " "
	@echo "To produce output"
	@echo "       $$ make output"
	@echo " "
	@echo "To create plots using VisClaw"
	@echo "       $$ make plots"
	@echo " "
	@echo "To create KMZ file for plotting in Google Earth"
	@echo "       $$ make kml"
	@echo " "
	@echo "To clean data and topo files"
	@echo "       $$ make clean"
