L = FM_files
OBJS1 = $(L)/fmsave.o $(L)/fmzm90.o $(L)/fm.o
OBJS2 = $(L)/TestFM.o $(L)/SampleFM.o
OBJS = $(OBJS1) $(OBJS2)
SRCS = $(L)/fmsave.f95 $(L)/fmzm90.f95 $(L)/fm.f95 $(L)/TestFM.f95 $(L)/SampleFM.f95
FMLIB = libfmlib.a
PROGS = TestFM SampleFM
FC = gfortran-4.7
FFLAGS = -O

.SUFFIXES: .f90 .o
	.f90.o:
	$(FC) $(FFLAGS) -c $<

all: $(PROGS)

fmlib.a: $(OBJS1)
	ar cr $@ $(OBJS1)
	ranlib $@

TestFM: $(L)/TestFM.o fmlib.a
	$(FC) $(FFLAGS) -o $@ $(L)/TestFM.o fmlib.a

SampleFM: $(L)/SampleFM.o $(L)/fmlib.a
	$(FC) $(FFLAGS) -o $@ $(L)/SampleFM.o fmlib.a

test: TestFM
	$(L)//TestFM