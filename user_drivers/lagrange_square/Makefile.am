noinst_PROGRAMS= \
square0 \
square1 \
square2 \
square3 \
square4


#------------------------------------------------------------------------
square0_SOURCES = square0.cc

square0_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------
square1_SOURCES = square1.cc

square1_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------
square2_SOURCES = square2.cc

square2_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------
square3_SOURCES = square3.cc

square3_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------
square4_SOURCES = square4.cc

square4_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)



#
#


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES = -I@includedir@  



#EXTRA_DIST = \
#fluid_iliac.1.ele   fluid_iliac.1.node  \
#solid_iliac.1.ele   solid_iliac.1.node \
#solid_iliac.1.face  solid_iliac.poly \
#fluid_iliac.1.face  fluid_iliac.poly  \
#fluid_iliac_with_extensions.1.ele   fluid_iliac_with_extensions.1.node  \
#solid_iliac_with_extensions.1.ele   solid_iliac_with_extensions.1.node \
#solid_iliac_with_extensions.1.face  solid_iliac_with_extensions.poly \
#fluid_iliac_with_extensions.1.face  fluid_iliac_with_extensions.poly  \
#iliac.xda    iliac_with_extensions.xda

# Extra actions for clean and distclean
distclean-local:
	touch Validation
	rm -r -f Validation
	rm -f *~
clean-local:
	touch Validation
	rm -r -f Validation
	rm -f *~
