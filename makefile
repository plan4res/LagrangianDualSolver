##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of LagrangianDualSolver                                         #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by LagrangianDualSolver, i.e., core SMS++.              #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by LagrangianDualSolver.                                            #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(LgDSLVSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(LgDSLVOBJ)   = the final object(s) / library                   #
#           $(LgDSLVH)     = the .h files to include                         #
#           $(LgDSLVINC)   = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone                               #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

LgDSLVOBJ = $(LgDSLVSDR)LagrangianDualSolver.o 

LgDSLVINC = -I$(LgDSLVSDR)

LgDSLVH   = $(LgDSLVSDR)LagrangianDualSolver.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(LgDSLVOBJ) $(LgDSLVSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(LgDSLVSDR)LagrangianDualSolver.o: $(LgDSLVSDR)LagrangianDualSolver.cpp \
	$(LgDSLVH) $(SMS++OBJ)
	$(CC) -c $*.cpp -o $@ $(LgDSLVINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
