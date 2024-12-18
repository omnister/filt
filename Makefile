OBJECTS=filter.o
CCFLAGS=-O -lm
DEMOS=tests/filterdemo tests/demobrf tests/demolpf tests/demolpt \
      tests/demobpf1 tests/demolpt2 tests/demobpf2 tests/demos_sdw
ALL=filter.c bpf brf hpf lpf.1 bessel.h Makefile $(DEMOS)

SYSBIN=/usr/local/bin/
MANDIR=/usr/local/man/man1/
CATDIR=/usr/local/man/cat1/


all: lpf bpf hpf brf

lpf: $(OBJECTS) bessel.h
	cc $(OBJECTS) -o lpf $(CCFLAGS) 
hpf: lpf 
	ln -s lpf hpf
bpf: lpf
	ln -s lpf bpf
brf: lpf
	ln -s lpf brf


install: lpf bpf hpf brf
	-cp lpf $(SYSBIN)
	-ln -s $(SYSBIN)/lpf $(SYSBIN)/bpf
	-ln -s $(SYSBIN)/lpf $(SYSBIN)/hpf
	-ln -s $(SYSBIN)/lpf $(SYSBIN)/brf

	-cp lpf.1 $(MANDIR)
	-ln -s $(MANDIR)/lpf.1 $(MANDIR)/bpf.1
	-ln -s $(MANDIR)/lpf.1 $(MANDIR)/hpf.1
	-ln -s $(MANDIR)/lpf.1 $(MANDIR)/brf.1

	/bin/rm -f $(CATDIR)/lpf.1
	/bin/rm -f $(CATDIR)/bpf.1
	/bin/rm -f $(CATDIR)/hpf.1
	/bin/rm -f $(CATDIR)/brf.1

tar: $(ALL) 
	tar -cvf filt.tar $(ALL)

.c.o:
	cc -c $(CCFLAGS) $<
