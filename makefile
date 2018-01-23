SRC = src
OBJ = obj
MOD = $(SRC)/mod
OUT = Xall
TARGET = $(HOME)/bin/$(OUT)

#===============================================================================
# Compilatore e opzioni di compilazione
FORT = mpifort

# di default usa le opzioni per il PGI
LDFLAGS = -fast -mp -Mcache_align -Mdalign -Mdclchk -Kieee -Ktrap=fp -i4

FCFLAGS = $(LDFLAGS) -module $(OBJ) -c

ifdef gnu
   FORT = mpifort # gfortran-6
   LDFLAGS = -O2 -Wall -J$(OBJ) -I$(OBJ) -fopenmp 
      ifdef deb
      LDFLAGS = -g -Wall -J$(OBJ) -I$(OBJ) -fopenmp -ffpe-trap=invalid,zero,overflow,underflow \
                -fbounds-check -fbacktrace -Warray-temporaries  -fcheck-array-temporaries
      endif
   FCFLAGS = $(LDFLAGS) -fimplicit-none -c 
endif

ifdef intel
FORT = mpif90
LDFLAGS = -O2 -openmp -module $(OBJ) -I/opt/mpi/openmpi-1.6.5-intel/lib/
   ifdef deb
       LDFLAGS = -g -check all -traceback -check bounds -debug all -openmp \
                 -module $(OBJ) -I/opt/mpi/openmpi-1.6.5-intel/lib/
   endif
FCFLAGS = $(LDFLAGS) -warn all -c
endif

#===============================================================================
# file sorgenti
SRC1 = $(MOD)/prec.f03
SRC2 = $(MOD)/var_type.f03 $(MOD)/cell_type.f03 $(MOD)/block_type.f03 $(MOD)/physical_par.f03
SRC3 = $(MOD)/proc_pointers.f03 $(MOD)/solution_par.f03 $(MOD)/motion_par.f03 $(MOD)/work_var.f03 $(MOD)/concon.f03
SRC4 = $(filter-out $(SRC1) $(SRC2) $(SRC3), $(wildcard $(MOD)/*.f03))
SRC_MAIN = $(wildcard $(SRC)/*.f03)

OBJ1 = $(patsubst $(MOD)/%.f03,$(OBJ)/%.o,$(SRC1))
OBJ2 = $(patsubst $(MOD)/%.f03,$(OBJ)/%.o,$(SRC2))
OBJ3 = $(patsubst $(MOD)/%.f03,$(OBJ)/%.o,$(SRC3))
OBJ4 = $(patsubst $(MOD)/%.f03,$(OBJ)/%.o,$(SRC4))
OBJ_MAIN = $(patsubst $(SRC)/%.f03,$(OBJ)/%.o,$(SRC_MAIN))

#===============================================================================
# Rules

xall : $(OBJ_MAIN)
	@rm -f $(TARGET)
	@echo Assemblaggio $(TARGET)
	@$(FORT) $(OBJ_MAIN) $(OBJ4) $(OBJ3) $(OBJ2) $(OBJ1) $(LDFLAGS) -o $(TARGET)
	@echo

$(OBJ_MAIN) : $(OBJ)/%.o: $(SRC)/%.f03 $(OBJ4) 
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(SRC)/$(<F) -o $@

$(OBJ4) : $(OBJ)/%.o: $(MOD)/%.f03 $(OBJ3) 
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(MOD)/$(<F) -o $@

$(OBJ3) : $(OBJ)/%.o: $(MOD)/%.f03 $(OBJ2)
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(MOD)/$(<F) -o $@

$(OBJ2) : $(OBJ)/%.o: $(MOD)/%.f03 $(OBJ1)
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(MOD)/$(<F) -o $@

$(OBJ1) : $(OBJ)/%.o: $(MOD)/%.f03 /bin/ls
	@echo Compilazione $(<F)
	@$(FORT) $(FCFLAGS) $(MOD)/$(<F) -o $@

#
# /bin/ls è, in teoria, un qualunque file che non viene mai aggiornato
#
# la sequenza così costruita serve a eseguire "opzioni" prima di cominciare la
# compilazione
#
# se non viene inserito il passaggio per /bin/ls e si fa dipendere OBJ1
# direttamente da "opzioni", ogni volta che si esegue il make viene ricompilato
# tutto il codice
#
/bin/ls : opzioni

opzioni :
	@echo
	@echo ++++++++++++++++++++
	@echo "Sorgenti in $(SRC)"
	@echo
	@echo "Compilatore:"
	@echo "   " `type $(FORT)` " Wrapper di " ${FC}
	@echo
	@echo "Opzioni di compilazione:"
	@echo "   $(FCFLAGS)"
	@echo
	@echo "Opzioni di linking:"
	@echo "   $(LDFLAGS)"
	@echo ++++++++++++++++++++
	@echo

clean :
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	rm -f $(TARGET)

cleanall :
	find $(OBJ) -type f -exec rm {} \;
	find $(HOME)/bin -name "${OUT}*" -exec rm {} \;
