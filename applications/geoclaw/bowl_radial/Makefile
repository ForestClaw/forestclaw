PYTHON=python
EXE=bowl

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

