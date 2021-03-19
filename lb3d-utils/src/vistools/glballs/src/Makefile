BASEPATH=..


include ../Makefile.config

TARGETS=glballs generator
TARGETS_SRC=./glballs.cc ./generator.cc

SRCS=$(TARGETS_SRC)\

OBJS=glballs.o ../D/*.o ../T/*.o

.PHONY : all install clean distclean test

all : .depend $(TARGETS) ../Makefile.config Makefile

glballs: $(OBJS) 
	$(CXX) $(LFLAGS) -L$(QTLIBDIR) -lqt-mt -L$(GLULIBDIR) -lGLU -L$(GLLIBDIR) -lGL -o glballs $(OBJS)

clean:
	-rm -f *~ core .depend
	-rm -f generator glballs
	-rm -f QuadPlane QuadPlane.o

distclean: clean

generator: generator.o
	$(CXX) -lm generator.o -o generator 

QuadPlane: ./QuadPlane.cc
	$(CXX) -lm -DWITH_QT $(CXXFLAGS) $(CPPFLAGS) ./QuadPlane.cc -o QuadPlane

QuadSurface: ./QuadSurface.cc
	$(CXX) -DWITH_QT $(CXXFLAGS) $(CPPFLAGS) -lm ./QuadSurface.cc -o QuadSurface

ifeq (.depend,$(wildcard .depend))
include .depend
endif
