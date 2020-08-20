PYTHON=python

EXE=chile2010

topo:
	${PYTHON} maketopo.py  # Essentially maketopo.py

data:
	${PYTHON} setrun.py  # Essentially setrun.py

output:
	@echo "Usage : "
	@echo "    To run example in serial : "
	@echo "           $$ ${EXE}"
	@echo "     "
	@echo "    To run example in parallel : "
	@echo "           $$ mpirun -n 4 ${EXE}"
	@echo "     "

	# ${EXE}

plots: 
	${PYTHON} setplot.py   # Essentially setplot.py

kml: 
	${PYTHON} setplot_kml.py     # Version of setplot.py for creating KML files.
