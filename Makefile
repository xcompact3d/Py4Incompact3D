#
#        FILE: Makefile
#      AUTHOR: Paul Bartholomew <ptb08@ic.ac.uk>
# DESCRIPTION: Makefile for Py4Incompact3D: runs tests and builds documentation.
#

all: test

test:
	python3 -m unittest discover

.PHONY: doc
doc:
	make -C doc/ latexpdf

