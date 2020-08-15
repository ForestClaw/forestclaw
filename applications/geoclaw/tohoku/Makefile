PYTHON=python

topo:
	${PYTHON} make_topo.py

data:
	${PYTHON} make_data.py

output:
	@echo "Usage : "
	@echo "To run example in serial : "
	@echo "       $$ tohoku"
	@echo " "
	@echo "To run example in parallel : "
	@echo "       $$ mpirun -n 4 tohoku"

plots: 
	${PYTHON} make_plots.py

kml: 
	${PYTHON} make_kml.py
