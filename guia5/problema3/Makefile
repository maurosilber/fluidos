ROOT := ../..
CONDA_BIN := micromamba
CONDA_ENV := $(ROOT)/conda
SOLVER := MHDB

.PHONY: install
install: $(CONDA_ENV)

$(CONDA_ENV): $(ROOT)/environment.GHOST.yaml
	if [ ! -d $(CONDA_ENV) ]; \
	then \
		$(CONDA_BIN) create --prefix $(CONDA_ENV) --file $(ROOT)/environment.GHOST.yaml; \
	else \
		$(CONDA_BIN) update --prefix $(CONDA_ENV) --file $(ROOT)/environment.GHOST.yaml; \
	fi
	touch $(CONDA_ENV)

SOLVER: src/*
	cp -r src/ $(ROOT)/GHOST/3D
	cd $(ROOT)/GHOST/3D/src && make
	mv $(ROOT)/GHOST/3D/bin/README_output.txt .
	mv $(ROOT)/GHOST/3D/bin/$(SOLVER) SOLVER
	cd $(ROOT)/GHOST/3D/src && make clean

.PHONY: clean
clean:
	cd output/ && rm *.txt *.out

output/vx.0001.out: output/parameter.inp
	-$(MAKE) clean
	-rm output/SOLVER
	ln SOLVER output/SOLVER
	cd output/ && time ./SOLVER > /dev/null
