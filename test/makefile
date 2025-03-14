##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for core SMS++)                                        #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
# NAME can be both main and test, the corresponding .cpp is compiled
# but the final executable is always called "test"
NAME = test

# basic directory
DIR = ./

# here is where the specific MCFClass to be tested is decided
# debug switches
SW = -g -std=c++17 -ferror-limit=1 -lz -lbz2
# production switches
#SW = -O3 -DNDEBUG -std=c++17

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

all default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# LagrangianDualSolver
LDSSLVSDR = ../
include $(LDSSLVSDR)makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(LDSSLVOBJ)

# libreries
MLIB =  $(LDSSLVLIB)

$(NAME): $(MOBJ) $(DIR)$(NAME).o
	$(CC) -o test $(DIR)$(NAME).o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(LDSSLVINC)

# includes
MH =    $(BNDSLVH)

# compile command

$(DIR)$(NAME).o: $(DIR)$(NAME).cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
