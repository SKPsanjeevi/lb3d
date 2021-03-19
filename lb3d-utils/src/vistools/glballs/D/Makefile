BASEPATH=..


include ../Makefile.config

OBJMOC=	\
 D_CGLBox_moc.o\
 D_CDisplay_moc.o\
 D_CLoader_moc.o

SRCMOC=$(patsubst %.o,%.cc,$(OBJMOC))

OBJS= \
 D_CGLBox.o\
 D_CGL_Circle.o\
 D_CGL_Cube.o\
 D_CGL_Sphere.o\
 D_CGL_Texture.o\
 D_CGL_Blip.o\
 D_CDisplay.o\
 D_CLoader.o\
 D_CGL.o\
 D_CGL_Borders.o

SRCS=$(patsubst %.o,./%.cc,$(OBJS))

####### Build rules

all: .depend $(OBJS) $(SRCMOC) $(OBJMOC)

clean:
	-rm -f $(OBJS) $(OBJMOC) $(SRCMOC)
	-rm -f *~ .depend

distclean: clean

ifeq (.depend,$(wildcard .depend))
include .depend
endif
