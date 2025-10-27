PETSC_DIR ?= petsc
PETSC_ARCH ?= arch-linux-opt

CXX = mpicxx
SRC = main.cpp
OUT = main.exe

.PHONY: deps build run

INCLUDES = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIBS     = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -Wl,-rpath=$(PETSC_DIR)/$(PETSC_ARCH)/lib

build: deps
	@echo ">> Installazione completata"

$(OUT): $(SRC)
	$(CXX) $(SRC) -o $(OUT) $(INCLUDES) $(LIBS)

run: $(OUT)
	mpirun -n 1 ./$(OUT)

clean:
	rm -f $(OUT)

deps:
	@echo ">> Installazione della libreria PETSc"
	git clone -b release https://gitlab.com/petsc/petsc.git petsc
	cd petsc && ./configure \
  	--with-cc=mpicc \
  	--with-cxx=mpicxx \
  	--with-fc=mpif90 \
  	--download-fblaslapack \
  	--with-debugging=0
	cd petsc && make all all -j$(nproc)
	cd petsc && make check

