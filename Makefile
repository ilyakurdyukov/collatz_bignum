
APPNAME := collatz_test
SRCNAME := main.c

MFLAGS := -march=native
CFLAGS := -Wall -Wextra -O3 $(MFLAGS)
SFLAGS := -fno-asynchronous-unwind-tables -masm=intel
PROF_ARGS := --ones 100000
CLANG := 0

.PHONY: clean all profile check

all: $(APPNAME)

clean:
	rm -f $(APPNAME)

$(APPNAME): $(SRCNAME)
	$(CC) -DAPPNAME=$(APPNAME) $(CFLAGS) -s -o $@ $<

$(APPNAME).s: $(SRCNAME)
	$(CC) -DAPPNAME=$(APPNAME) $(CFLAGS) $(SFLAGS) -S -o $@ $<

profile:
	@$(MAKE) --no-print-directory MFLAGS="$(MFLAGS) -fprofile-generate" $(APPNAME)
	./$(APPNAME) $(PROF_ARGS)
	rm $(APPNAME)
ifneq ($(CLANG),0)
	llvm-profdata merge -output=default.profdata *.profraw
	rm -f *.profraw
endif
	@$(MAKE) --no-print-directory MFLAGS="$(MFLAGS) -fprofile-use" $(APPNAME)

check: $(APPNAME)
	./$(APPNAME) --num 989345275647 | \
		grep -q "mul3 = 506, div2 = 842, total = 1348$$"
	./$(APPNAME) --lut 20 --ones 100000 | \
		grep -q "mul3 = 481603, div2 = 863323, total = 1344926$$"


