all: parameters run

parameters:
	python parameters.py

run: output/*
	for d in $^; do $(MAKE) $$d/vx.0001.out; done;

clean/%:
	cd $* && rm *.txt *.out || true

output/%/vx.0001.out: output/%/parameter.inp
	$(MAKE) clean/output/$*
	cd output/$* && time ./SOLVER > /dev/null
