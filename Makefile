.PHONY: help prepare-dev clean

VENV_NAME?=venv
VENV_ACTIVATE=$(VENV_NAME)/bin/activate
PYTHON=${VENV_NAME}/bin/python3

.DEFAULT: help
help:
	@echo "make prepare-dev"
	@echo "	   preparing the development environment (use only once with sudo)"
	@echo "	    - install python3, pip and virtualenv"
	@echo "	    - create virtual environment inside the current directory"
	@echo "make venv"
	@echo "	    creating the virtual environment"
	@echo "	     - install required packages"
	@echo "make clean"
	@echo "	   remove virtualenv and .egg-info"

prepare-dev:
	sudo apt-get -y install python3 python3-pip virtualenv
	make venv

venv: $(VENV_ACTIVATE)
$(VENV_ACTIVATE): requirements.txt
	test -d $(VENV_NAME) || virtualenv -p python3 $(VENV_NAME)
	${PYTHON} -m pip install -U pip
	. $(VENV_ACTIVATE); pip install -Ur requirements.txt
	touch $(VENV_ACTIVATE)

clean:
	test -r $(VENV_NAME) && rm -r $(VENV_NAME)