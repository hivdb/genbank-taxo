start:
	python -m venv env

all:
	./env/bin/python work.py

.PHONY:
	start all