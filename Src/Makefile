#
# 				Makefile
#                              Mark NOBLE
#                       Ecole des Mines de Paris
#
# created first 04.09.97 M.N.    (Solaris 2.5, Linux RedHat)
######################################################################
# Instructions: (in Makefile.config)
#	Make changes where necessary to reflect the needs of your
#	system, by commenting and uncommenting necessary lines
######################################################################
#                     Include options for your machine
include ./Makefile.config
#
######################################################################
# Instructions: Add in programs to be compiled
#
# Place to put Binaries
BINDIR=.
#
# Post-stack Depth 1D depth migration
#
OBJ1 = invdfe.o Time_2d.o
invdfe:$(OBJ1)
	$(FC) -o $(BINDIR)/$@ $(OBJ1) $(LDFLAGS)
