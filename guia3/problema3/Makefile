.PHONY: all
all: parameters run

SOLVER: src/*
	cp -r src/ ../../GHOST/3D
	cd ../../GHOST/3D/src && make
	mv ../../GHOST/3D/bin/README_output.txt .
	mv ../../GHOST/3D/bin/BOUSS SOLVER

.PHONY: parameters
parameters: parameter.inp parameters.py
	python parameters.py

.PHONY: run
run: output/*
	for d in $^; do $(MAKE) $$d/vx.0001.out; done;

.PHONY: clean
clean/%:
	cd $* && rm *.txt *.out

output/%/vx.0001.out: output/%/parameter.inp
	-$(MAKE) clean/output/$*
	-rm output/$*/SOLVER
	ln SOLVER output/$*/SOLVER
	cd output/$* && time ./SOLVER > /dev/null
