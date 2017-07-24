.DEFAULT_GOAL := all

MPICC = mpiCC
MPICFLAGS = -std=c++11
MPICOPTFLAGS = -O3 -g -lpng
MPILDFLAGS =

TARGETS = mandelbrot_serial$(EXEEXT) mandelbrot_joe$(EXEEXT) mandelbrot_susie$(EXEEXT) mandelbrot_ms$(EXEEXT)

all: $(TARGETS)

serial: mandelbrot_serial$(EXEEXT)

susie: mandelbrot_susie$(EXEEXT) mandelbrot_serial$(EXEEXT)

joe: mandelbrot_joe$(EXEEXT)

ms: mandelbrot_ms$(EXEEXT)

SRCS_COMMON = render.cc 

DISTFILES += $(SRCS_COMMON) $(DEPS_COMMON)

mandelbrot_serial$(EXEEXT): mandelbrot_serial.cc $(SRCS_COMMON) $(DEPS_COMMON)
	$(MPICC) $(MPICFLAGS) $(MPICOPTFLAGS) -I/data/apps/boost/1.57/include \
	    -o $@ mandelbrot_serial.cc $(SRCS_COMMON) $(MPILDFLAGS)

mandelbrot_joe$(EXEEXT): mandelbrot_joe.cc $(SRCS_COMMON) $(DEPS_COMMON)
	$(MPICC) $(MPICFLAGS) $(MPICOPTFLAGS) -I/data/apps/boost/1.57/include \
	    -o $@ mandelbrot_joe.cc $(SRCS_COMMON) $(MPILDFLAGS)

mandelbrot_susie$(EXEEXT): mandelbrot_susie.cc $(SRCS_COMMON) $(DEPS_COMMON)
	$(MPICC) $(MPICFLAGS) $(MPICOPTFLAGS) -I/data/apps/boost/1.57/include \
	    -o $@ mandelbrot_susie.cc $(SRCS_COMMON) $(MPILDFLAGS)

mandelbrot_ms$(EXEEXT): mandelbrot_ms.cc $(SRCS_COMMON) $(DEPS_COMMON)
	$(MPICC) $(MPICFLAGS) $(MPICOPTFLAGS) -I/data/apps/boost/1.57/include \
	    -o $@ mandelbrot_ms.cc $(SRCS_COMMON) $(MPILDFLAGS)

clean:
	rm -f $(TARGETS) *.png Mandelbrot.o* Mandelbrot.e* core.*

# eof
