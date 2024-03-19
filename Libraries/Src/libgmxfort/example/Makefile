.PHONY: all clean

all:
	@gfortran main.f90 `pkg-config --cflags --libs libgmxfort` -Wall

clean:
	@rm a.out
