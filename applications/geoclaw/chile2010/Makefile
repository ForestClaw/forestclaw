PYTHON=python

EXE=chile2010

topo:
	${PYTHON} maketopo.py  # Essentially maketopo.py

data:
	${PYTHON} make_data.py  # Essentially setrun.py

output:
	@echo "Usage : "
	@echo "To run example in serial : "
	@echo "       $$ ${EXE}"
	@echo " "
	@echo "To run example in parallel : "
	@echo "       $$ mpirun -n 4 ${EXE}"
	@echo " "
	@echo "Running the serial example using options in fclaw_options.ini"
	${EXE}

plots: 
	${PYTHON} make_plots.py   # Essentially setplot.py

kml: 
	${PYTHON} make_kml.py     # Version of setplot.py for creating KML files.
