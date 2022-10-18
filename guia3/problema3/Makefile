OUTPUT := $(wildcard output/*)

.PHONY: all
all: parameters $(OUTPUT)

SOLVER: src/*
	cp -r src/ ../../GHOST/3D
	cd ../../GHOST/3D/src && make
	mv ../../GHOST/3D/bin/README_output.txt .
	mv ../../GHOST/3D/bin/BOUSS SOLVER

.PHONY: parameters
parameters: parameter.inp parameters.py
	python parameters.py

.PHONY: $(OUTPUT)
$(OUTPUT):
	$(MAKE) $@/vx.0001.out

clean/%:
	cd $* && rm *.txt *.out

output/%/vx.0001.out: output/%/parameter.inp
	-$(MAKE) clean/output/$*
	-rm output/$*/SOLVER
	ln SOLVER output/$*/SOLVER
	cd output/$* && time ./SOLVER > /dev/null
