FC = mpif90
OLEVEL = -O3 -r8
#OLEVEL = -O2 -fbounds-check -g -fbacktrace -fdump-core -ffpe-trap=zero,invalid,overflow -Wuninitialized -Wall
RM = rm -f

#----------------------------------------------------------------------
# Base code
#----------------------------------------------------------------------
MAIN = main.o
OBJ_FILES = alloc.o calc.o computemetrics.o constants.o postprocess.o readfiles.o setup.o startmpi.o utility.o writefiles.o

OBJ_MODS = mod_post.o 
OBJS = $(OBJ_MODS) $(OBJ_FILES) $(MAIN)
EXEC = blpost.exe
LDFLAGS =
#----------------------------------------------------------------------

$(EXEC): $(OBJS) Makefile 
	$(FC) $(PRECFLAGS) $(OLEVEL) -o $@ $(OBJS) $(FCLIBS) $(LDFLAGS)

%.o:%.F90 Makefile 
	$(FC) $(PRECFLAGS) $(INCLUDE) $(OLEVEL) -c $< -o $@

%.o:%.cpp Makefile
	$(CXX) $(PRECFLAGS) $(INCLUDE) -O2 -c $< -o $@

%.o:%.cxx Makefile
	$(CXX) $(PRECFLAGS) $(INCLUDE) -O2 -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(EXEC) $(MAIN) $(OBJS) $(OBJ_MODS) $(GRID) *.mod
