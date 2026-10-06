CONDA_ENV_DIR := $(shell dirname $(CONDA_EXE))
CONDA_ENV_NAME := literature_harmonization

LIFTOVER_KEY ?= default

TARGETS := dependencies add_studies clean_table liftover_bcftools liftover_test harmonization gwasstudio_files

all:
	@echo "Available targets: $(TARGETS)"

dependencies:
	conda env update --file environment.yml

dependencies-force:
	conda env update --file environment.yml --prune
	@echo "Forcing reinstall of gwaspipe from main..."
	conda run --no-capture-output -n "$(CONDA_ENV_NAME)" \
		python -m pip install --force-reinstall --no-cache-dir --no-deps \
		"git+https://github.com/ht-diva/gwaspipe.git@main"

add_studies:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/add_studies.py

liftover_bcftools:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/literature_liftover_bcftools.py

liftover_test:
	$(MAKE) harmonization LIFTOVER_KEY=bcftools
	$(MAKE) liftover_bcftools
	$(MAKE) harmonization LIFTOVER_KEY=gwaslab_standard
	$(MAKE) harmonization LIFTOVER_KEY=gwaslab_bridge
	python liftover_comparison.py

clean_table:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/clean_literature_table.py

harmonization:
	source "$(CONDA_ENV_DIR)/activate" "$(CONDA_ENV_NAME)" && \
	LIFTOVER_KEY="$(LIFTOVER_KEY)" python scripts/literature_harmonization.py

gwasstudio_files:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/make_gwasstudio_files.py

combine_output:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/combine_output.py

.PHONY: all dependencies add_studies clean_table liftover_bcftools liftover_test harmonization gwasstudio_files
