

DIRS = T D src

all: stamp-h Makefile Makefile.config
	$(RM) make.all
	for i in $(DIRS); do $(MAKE) -C $$i all 2>&1 | tee -a make.all; done
clean: Makefile
	$(RM) make.clean
	for i in $(DIRS); do $(MAKE) -C $$i clean 2>&1 | tee -a make.clean; done
distclean: Makefile
	$(RM) make.distclean
	for i in $(DIRS) tools doc; do $(MAKE) -C $$i distclean 2>&1 | tee -a make.distclean; done
	$(RM) config.cache config.log config.h make.all make.install make.clean make.distclean make.doc stamp-h gmon.out

./configure : ./configure.in
	cd . && autoconf

./config.h.in : ./stamp-h.in
./stamp-h.in: ./configure.in
#	cd . && autoheader
	date > ./stamp-h.in

config.h: stamp-h
stamp-h: ./config.h.in config.status
	./config.status

Makefile Makefile.config: ./stamp-h.in ./config.h.in ./Makefile.in ./Makefile.config.in config.status
	./config.status

config.status: ./configure
	./config.status --recheck
