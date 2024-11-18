.SUFFIXES:
EXEC = test_pzgemm test_pcgemm test_pdgemm test_psgemm 
SRC = scalapack_interface.F90 test_pvgemm.F90
obj_path = _obj_d

ifdef MKLROOT
# Intel ifx using MKL ScaLapack
ifdef STATICSCALAPACK
LIB = -L$(HOME)/src/scalapack -lscalapack
else
ifdef ILP
LIB = -lmkl_scalapack_ilp64  -lmkl_blacs_intelmpi_ilp64
else
LIB = -lmkl_scalapack_lp64  -lmkl_blacs_intelmpi_lp64
endif
endif
#FPP = ifort -E
#FF  = mpiifort -c
#LD  = mpiifort
FPP = ifx -E
FF  = mpiifx -c
LD  = mpiifx
FFLAGS = -warn all -module $(obj_path) -I./$(obj_path)
ifdef ILP
FFLAGS := -i8 -qmkl-ilp64=parallel $(FFLAGS)
LDFLAGS = -qmkl-ilp64=parallel
else
FFLAGS := -qmkl=parallel $(FFLAGS)
LDFLAGS = -qmkl=parallel
endif
else
# GNU gfortran
# On Debian: 
# apt install gfortran 
# apt install libopenmpi-dev 
# apt install libscalapack-openmpi2.1 libscalapack-openmpi2.1-dev
LIB = -lscalapack-openmpi 
FPP = gfortran -E
FF = mpifort -c -J $(obj_path)
LD = mpifort
FFLAGS = -Wall -I./$(obj_path) $(INCLUDE)
LDFLAGS = 
endif


all: $(EXEC)

test_pdgemm: PPFLAGS = -DfD
test_psgemm: PPFLAGS = -DfS
test_pzgemm: PPFLAGS = -DfZ
test_pcgemm: PPFLAGS = -DfC

$(EXEC): % : $(SRC) 
	-mkdir -p $(obj_path)
	for i in $(subst .F90,,$^) ; do $(FPP) $(PPFLAGS) $$i.F90 > $(obj_path)/$$i.f90 ; done 
	for i in $(subst .F90,,$^) ; do $(FF) $(FFLAGS) -o $(obj_path)/$$i.o $(obj_path)/$$i.f90 ; done
	$(LD) $(LDFLAGS) -o $@ $(addprefix $(obj_path)/, $(subst .F90,.o,$+)) $(LIB) 


test_d: test_pdgemm
test_s: test_psgemm
test_z: test_pzgemm
test_c: test_pcgemm

test_s test_d test_c test_z: 
	@for nodes in 1 2 3 4 5 6 ; do         \
	echo '-------------------'           ; \
	echo "Test with np = $$nodes"        ; \
	echo '-------------------'           ; \
	echo mpirun -n $$nodes ./$<          ; \
	mpirun -n $$nodes ./$<               ; \
	echo ""                              ; \
	done

.PHONY: clean realclean

clean:
	-rm -f $(obj_path)/*.o 
	-rm -f $(obj_path)/*.mod 

realclean: clean
	-rm -f $(EXEC)

