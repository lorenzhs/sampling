CXX ?= g++

MKLROOT ?= /opt/intel/compilers_and_libraries/linux/mkl
MKL ?= ${MKLROOT}/lib/intel64
MKLFLAGS = -L${MKL} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
MKLLD = -Wl,-Bstatic ${MKLFLAGS} -Wl,-Bdynamic -ldl

LDFLAGS=-lpthread libstocc.a
CFLAGS=-std=c++14 -Wall -Wextra -Werror -g
OPT=-Ofast -DNDEBUG -march=native
DEBUG=-O0 -march=native -fsanitize=address
CFLAGS+=-IDistributedSampling/lib -IDistributedSampling/lib/tools -IDistributedSampling/extern/stocc -IDistributedSampling/extern/dSFMT

PGO=0
ifeq ("clang","$(findstring clang,$(CXX))")
        OPT += -flto
else ifeq ("icc","$(findstring icc,$(CXX))")
        OPT += -ipo
else
        OPT += -flto=8
        PGO=1
endif

flags ?= # runtime flags

# use 64-bit integers
# explicit 32-bit integer suffix
ifeq ($(B64),0) # we need MKL but remove std_gen
SUFF:=32
CFLAGS+=-DNOSTD -I${MKLROOT}/include
LDFLAGS+=${MKLLD}
else
ifneq ($(B64),) #64-bit mode
CFLAGS+=-DUSE64BIT
SUFF:=64
else # we need the MKL
CFLAGS+=-I${MKLROOT}/include
LDFLAGS+=${MKLLD}
endif
endif

# use stable fixer
ifneq ($(STABLE),)
CFLAGS+=-DFIX_STABLE
SUFF:=${SUFF}S
endif

GIT_VERSION := $(shell git describe --abbrev=7 --dirty --always --tags)
CFLAGS+=-DVERSION=\"$(GIT_VERSION)\"

.PHONY: rand pgo

rand: rand.cpp include/*.h
ifneq ($(SUFF),)
	@echo "Output name is rand${SUFF}"
endif
	${CXX} ${CFLAGS} ${OPT} -o rand${SUFF} rand.cpp ${LDFLAGS}

pgo: rand.cpp include/*.h
ifneq ($(SUFF),)
	@echo "Output name is rand${SUFF}-pgo"
endif
	rm -f rand.gcda
	${CXX} ${CFLAGS} ${OPT} -fprofile-generate -o rand${SUFF}-pgo rand.cpp ${LDFLAGS}
	./rand${SUFF}-pgo -i 1000 -k 65536
	./rand${SUFF}-pgo -i 100 -k 1048576
	./rand${SUFF}-pgo -i 10 -k 33554432
	${CXX} ${CFLAGS} ${OPT} -fprofile-use -o rand${SUFF}-pgo rand.cpp ${LDFLAGS}

debug:
ifneq ($(SUFF),)
	@echo "Output name is rand${SUFF}-dbg"
endif
	${CXX} ${CFLAGS} ${DEBUG} -o rand${SUFF}-dbg rand.cpp ${LDFLAGS}

run:
	@LD_LIBRARY_PATH=${MKL} ./rand ${flags}

buildall:
	make rand B64=0 &
	make rand B64=0 STABLE=1 &
	make rand B64=1 &
	make rand B64=1 STABLE=1 &

buildall-pgo:
ifeq ($(PGO),1)
	B64=0 make pgo
	B64=0 STABLE=1 make pgo
	B64=1 make pgo
	B64=1 STABLE=1 make pgo
endif

everything: buildall buildall-pgo
	make rand &
	STABLE=1 make rand &
ifeq ($(PGO),1)
	make pgo
	STABLE=1 make pgo
endif
