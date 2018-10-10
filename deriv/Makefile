#
# 
#

all: build

build:
	f2py3 tdma.f90 -m tdma -h tdma.pyf
	f2py3 -c tdma.pyf tdma.f90

test:
	make -C ../ test

clean:
	rm *.so tdma.pyf
