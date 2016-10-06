COMPLR = $(shell which clang >/dev/null; echo $$?)
ifeq "$(COMPLR)" "0"
CC=clang
else
CC=gcc
endif

#CC=gcc
ifeq ($(CC), clang)
    CFLAGS=-O3 -pedantic -Wall -Wextra -Weverything
else
    CFLAGS=-O3 -pedantic -Wall -Wextra -std=c99
endif

TGT=diead

$(TGT): main.o RNGTest.o RNGUtil.o
	@echo [$(CC)] compiling $@
	@$(CC) $(CFLAGS) $^ -o $@
	@echo [$(CC)] done.

main.o: main.c RNGTest.h RNGUtil.h
	@echo [$(CC)] compiling $@
	@$(CC) $(CFLAGS) -c -o $@ $<

RNGTest.o: RNGTest.c RNGTest.h
	@echo [$(CC)] compiling $@
	@$(CC) $(CFLAGS) -c -o $@ $<

RNGUtil.o: RNGUtil.c RNGUtil.h
	@echo [$(CC)] compiling $@
	@$(CC) $(CFLAGS) -c -o $@ $<

all: $(TGT)

.PHONY: clean

clean:
	@echo cleaning everything
	rm -f $(TGT) *.o
	@echo done.
