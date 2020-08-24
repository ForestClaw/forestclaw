PYTHON=python
EXE=tohoku

topo:
	${PYTHON} maketopo.py

data:
	${PYTHON} setrun.py

output:
	@echo "Usage (output) : "
	@echo "To run example in serial : "
	@echo "       $$ ${EXE}"
	@echo " "
	@echo "To run example in parallel : "
	@echo "       $$ mpirun -n 4 ${EXE}"

plots: 
	${PYTHON} setplot.py

kml: 
	${PYTHON} setplot_kml.py

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

