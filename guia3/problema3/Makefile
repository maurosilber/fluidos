all: parameters run

SOLVER: src/*
	cp -r src/ ../../GHOST/3D
	cd ../../GHOST/3D/src && make
	mv ../../GHOST/3D/bin/README_output.txt .
	mv ../../GHOST/3D/bin/BOUSS SOLVER

parameters: parameter.inp SOLVER
	python parameters.py

run: output/*
	for d in $^; do $(MAKE) $$d/vx.0001.out; done;

clean/%:
	cd $* && rm *.txt *.out || true

output/%/vx.0001.out: output/%/parameter.inp output/%/SOLVER
	$(MAKE) clean/output/$*
	cd output/$* && time ./SOLVER > /dev/null
