init:
	python3 -m venv env
	@echo "env" >> .gitignore
	./env/bin/pip3 install -r requirements.txt

all:
	./env/bin/python work.py

.PHONY:
	all start

.DEFAULT_GOAL := all
