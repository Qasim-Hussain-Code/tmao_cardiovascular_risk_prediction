.PHONY: help setup test analysis clean

help:
	@echo "setup     install the dependencies"
	@echo "test      run the test suite"
	@echo "analysis  run the full analysis and write tables and figures"
	@echo "clean     remove generated outputs and caches"

setup:
	python -m pip install -r requirements.txt

test:
	python -m pytest tests/

analysis:
	python scripts/run_analysis.py

clean:
	rm -rf results/figures/*.png results/tables/*.csv results/tables/*.json
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -rf .pytest_cache
