all: output/*
	for d in $^; do $(MAKE) $$d/vx.0001.out; done;

clean/%:
	cd $* && rm *.txt *.out

output/%/vx.0001.out: output/%/parameter.inp
	cd output/$* && time ./SOLVER > /dev/null
