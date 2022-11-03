CONDA_BIN := micromamba
CONDA_ENV := ./env

.PHONY: install
install: $(CONDA_ENV)

$(CONDA_ENV): environment.yaml
	if [ ! -d $(CONDA_ENV) ]; \
	then \
		echo hola; \
		$(CONDA_BIN) create --prefix $(CONDA_ENV) --file environment.yaml; \
	else \
		echo chau; \
		$(CONDA_BIN) update --prefix $(CONDA_ENV) --file environment.yaml; \
	fi
	touch $(CONDA_ENV)