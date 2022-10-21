ROOT := ../..
CONDA_BIN := mamba
OUTPUT := $(wildcard output/*)

.PHONY: all
all: parameters $(OUTPUT)

.PHONY: install
install: $(ROOT)/conda/

$(ROOT)/conda/: $(ROOT)/environment.yaml
	if [! -d $(ROOT)/conda]; \
	then \
		$(CONDA_BIN) env create --prefix $(ROOT)/conda --file $(ROOT)/environment.yaml; \
	else \
		$(CONDA_BIN) env update  --prefix $(ROOT)/conda --file $(ROOT)/environment.yaml; \
	fi 

SOLVER: src/*
	cp -r src/ $(ROOT)/GHOST/3D
	cd $(ROOT)/GHOST/3D/src && make
	mv $(ROOT)/GHOST/3D/bin/README_output.txt .
	mv $(ROOT)/GHOST/3D/bin/BOUSS SOLVER
	cd $(ROOT)/GHOST/3D/src && make clean

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
