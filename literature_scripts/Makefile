CONDA_ENV_DIR := $(shell dirname $(CONDA_EXE))
CONDA_ENV_NAME := literature_harmonization

TARGETS := dependencies add_studies liftover clean_table harmonization gwasstudio_files

all:
	@echo "Available targets: $(TARGETS)"

dependencies:
	conda env update --file environment.yml

add_studies:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/add_studies.py

liftover:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/literature_liftover.py

clean_table:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/clean_literature_table.py

harmonization:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/literature_harmonization.py

gwasstudio_files:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/make_gwasstudio_files.py

combine_output:
	source $(CONDA_ENV_DIR)/activate $(CONDA_ENV_NAME) && \
	python scripts/combine_output.py

.PHONY: all dependencies add_studies liftover clean_table harmonization gwasstudio_files
