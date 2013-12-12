#=================================================================================================#
# Makefile template for mixing language programming (C, C++, Fortran)
# Author: Zhihua Ma
# Email:  jackycfd@gmail.com
#=================================================================================================#
# A typical project's structure
#  project/
#  |-- bin
#  |-- doc
#  |-- include
#  |-- lib
#  |-- obj
#  '-- src
#      |-- C
#      |-- CPP
#      '-- F90
#=================================================================================================#
# The common compilers and flags (C, C++, Fortran)
# usually, you can keep the same
CC	= gcc
CXX     = g++
FC	= gfortran
CFLAGS 	= -c #-g 
FFLAGS	= -c #-g
OLEVEL  = -O3
#==================================================================================================#
# the directories
# you can keep these variables almost the same
DIR_BIN := bin
DIR_DOC := doc
DIR_INC := include
DIR_LIB := lib
DIR_OBJ := obj
DIR_SRC := src
DIR_SRC_C  := $(DIR_SRC)/C
DIR_SRC_CPP:= $(DIR_SRC)/CPP
DIR_SRC_F90:= $(DIR_SRC)/F90


# the include path
INCLUDE = -I$(DIR_INC)
#=================================================================================================#
# FILEs' name (executable, header, source, object and library)
# You must modify the variables according to your case

# the executable file name
# must not be empty
FILE_EXE := shock

# the header files' name
FILE_HEADER := Riemann.hh

# the source files' name
# at least one source file
FILE_C := 
FILE_CPP:= Riemann.cpp ExactShock.cpp NumericalShock.cpp main.cpp 
FILE_F90:= 

# the library file's name
# must not be empty
FILE_LIB:= libShock.a


#***************************************************************************************************#
# FILEs with path and name
# In very few cases you need to modify the variables
# the source files
SRC_C   := $(addprefix $(DIR_SRC_C)/,$(FILE_C))
SRC_CPP := $(addprefix $(DIR_SRC_CPP),$(FILE_CPP))
SRC_F90 := $(addprefix $(DIR_SRC_F90)/,$(FILE_F90))

#the header files
HEADER := $(addprefix $(DIR_INC)/,$(FILE_HEADER))

# the object files
OBJ_C  := $(addprefix $(DIR_OBJ)/,$(FILE_C:.c=.o))
OBJ_CPP:= $(addprefix $(DIR_OBJ)/,$(FILE_CPP:.cpp=.o))
OBJ_F90:= $(addprefix $(DIR_OBJ)/,$(FILE_F90:.f90=.o)) 
OBJ_ALL:= $(OBJ_C) $(OBJ_CPP) $(OBJ_F90)

# the executable file
EXE    := $(DIR_BIN)/$(FILE_EXE)

# the OpenGL library
ifeq ($(shell uname),CYGWIN_NT-5.2-WOW64)
	LIBGL= -lglut -lglu -lgl
endif
ifeq ($(shell uname),CYGWIN_NT-6.1)
	LIBGL= -lglut -lglu -lgl
endif
ifeq ($(shell uname),Linux)
	LIBGL= -lglut -lGLU -lGL
endif
ifeq ($(shell uname),Darwin)
	LIBGL = -framework OpenGL -framework GLUT
endif

# the library
LIB     := $(addprefix $(DIR_LIB)/,$(FILE_LIB))
LIBLINK := $(patsubst %.a,%,$(patsubst lib%,%,$(FILE_LIB)))
LIBSYS  := -lstdc++ $(LIBGL)
#***************************************************************************************************#
# In very few cases you need to modify the targes, pre and rules
# default target by convention is ''all''
.PHONY: ALL clean

ALL:  $(LIB) $(EXE) 

$(EXE):  $(OBJ_F90) $(LIB) 
#	$(FC)  $(OBJ_F90) -L./$(DIR_LIB)/ -l$(LIBLINK) -lstdc++  -o $@ $(OLEVEL)
#	$(FC) $^ $(LIBSYS) -g -o $@ 
	$(CXX) $^ $(LIBSYS) $(OLEVEL) -o $@

$(LIB): $(OBJ_C) $(OBJ_CPP)
	ar rv $@ $?

$(OBJ_C): $(HEADER)

$(OBJ_CPP): $(HEADER)

$(DIR_OBJ)/%.o: $(DIR_SRC_C)/%.c
	$(CC) $(INCLUDE) $(CFLAGS) $(OLEVEL) $< -o $@ 

$(DIR_OBJ)/%.o: $(DIR_SRC_CPP)/%.cpp
	$(CXX) $(INCLUDE) $(CFLAGS) $(OLEVEL) $< -o $@ 

$(DIR_OBJ)/%.o: $(DIR_SRC_F90)/%.f90
	$(FC) $(FFLAGS) $(OLEVEL) $< -o $@

ECHO:
	@echo FILE_C: $(FILE_C) 
	@echo FILE_CPP: $(FILE_CPP)  
	@echo FILE_F90: $(FILE_F90)
	@echo EXE: $(EXE) 
	@echo LIB: $(LIB)	
	@echo OBJ_F90: $(OBJ_F90)

clean:
	$(RM) $(OBJ_C) $(OBJ_CPP) $(OBJ_F90) $(LIB) 
