start:
	python -m venv env

all:
	./env/bin/python work.py

.PHONY:
	all start

.DEFAULT_GOAL := all
