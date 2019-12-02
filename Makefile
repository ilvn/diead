ifeq "$(shell which clang >/dev/null; echo $$?)" "0"
	CC := clang
else
	CC := gcc
endif

CFLAGS := -O3 -pedantic -Wall -Wextra -std=c99
ifeq ($(CC), clang)
    CFLAGS += -Weverything
endif

TGT := diead

.PHONY: clean all

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

clean:
	@echo cleaning everything
	rm -f $(TGT) *.o
	@echo done.
