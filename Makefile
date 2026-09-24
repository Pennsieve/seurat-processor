.PHONY: help test build clean run-local

SERVICE_NAME  ?= "seurat-processor"

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make build            - build the Docker image"
	@echo "make test             - run the processor locally via docker-compose"
	@echo "make run-local        - run the processor directly with Rscript (no Docker)"
	@echo "make clean            - remove output files and stop containers"

build:
	docker compose -f docker-compose.yml build

test:
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor

run-local:
	INPUT_DIR=./data/input OUTPUT_DIR=./data/output R_MAX_VSIZE=64Gb Rscript processor/main.R

clean:
	docker compose -f docker-compose.yml down --remove-orphans
	rm -rf data/output/*
