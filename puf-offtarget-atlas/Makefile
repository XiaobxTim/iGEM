.PHONY: install test lint typecheck check app example

install:
	python -m pip install -e ".[dev]"

test:
	pytest -q

lint:
	ruff check .

typecheck:
	mypy src

check: test lint typecheck

app:
	streamlit run app/streamlit_app.py

example:
	pufscan scan --query AACGUCUAUA --max-mismatches 1 --gencode-fasta tests/data/synthetic.fa --gencode-gtf tests/data/synthetic.gtf --gtex-expression tests/data/expression.tsv --target-tissue Liver --mode binding_only --no-structure --output results

