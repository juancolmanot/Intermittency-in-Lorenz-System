# ====================================================================
# Compiler info and paths
# ====================================================================
CC = gcc
OMP = -fopenmp
GDB = -g

# Directories
SRCDIR = src
INCDIR = include
BUILDDIR = build
COMMON_SRCDIR = ../../common/src
COMMON_INCDIR = ../../common/include

# Main program (to be set via command line or manually)
MAIN =

# Libraries
LDLIBS = -lgsl -lgslcblas -lopenblas -lm -fopenmp

# Compile flags
# ====================================================================
# Warning flags
flag_w01 = -Wall -Wextra -Wconversion -pedantic -Wno-unused-parameter
# Debugging flags
flag_d01 = -g
# Optimization flags
flag_o01 = -O3 -ftree-vectorize -ftree-loop-vectorize -funroll-loops
flag_o02 = -march=native -Ofast -ffast-math
flag_o03 = -fopt-info-vec -ftree-vectorizer-verbose=2

CFLAGS = ${flag_o01} ${flag_o02} ${flag_w01} ${flag_d01} ${flag_o03} ${OMP}

# Header files (local and common)
LOCAL_HEADERS_LIST = $(LOCAL_HEADERS)
COMMON_HEADERS_LIST = $(COMMON_HEADERS)

# Transform header list into .o and .h file lists
LOCAL_HEADERS_O := $(addsuffix .o,$(LOCAL_HEADERS_LIST))
LOCAL_HEADERS_H := $(addsuffix .h,$(LOCAL_HEADERS_LIST))
COMMON_HEADERS_O := $(addsuffix .o,$(COMMON_HEADERS_LIST))
COMMON_HEADERS_H := $(addsuffix .h,$(COMMON_HEADERS_LIST))

# Object files (for dependencies)
LOCAL_OBJS := $(LOCAL_HEADERS_O:%=$(BUILDDIR)/%)
COMMON_OBJS := $(COMMON_HEADERS_O:%=$(BUILDDIR)/%)

# All object files
OBJS = $(LOCAL_OBJS) $(COMMON_OBJS)

# Rules
.PHONY: all clean

all: $(BUILDDIR)/$(MAIN)

# Rule to build the main program
$(BUILDDIR)/$(MAIN): $(SRCDIR)/$(MAIN).o $(OBJS)
	@mkdir -p $(BUILDDIR)
	@$(CC) $(CFLAGS) -I$(INCDIR) -I$(COMMON_INCDIR) $^ -o $@ $(LDLIBS)

# Rule to compile local headers into object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(LOCAL_HEADERS_H:%=$(INCDIR)/%)
	@mkdir -p $(BUILDDIR)
	@$(CC) $(CFLAGS) -I$(INCDIR) -I$(COMMON_INCDIR) -c $< -o $@

# Rule to compile common headers into object files
$(BUILDDIR)/%.o: $(COMMON_SRCDIR)/%.c $(COMMON_HEADERS_H:%=$(COMMON_INCDIR)/%)
	@mkdir -p $(BUILDDIR)
	@$(CC) $(CFLAGS) -I$(INCDIR) -I$(COMMON_INCDIR) -c $< -o $@

# Clean up
clean:
	@rm -f $(BUILDDIR)/*

.PHONY: clean all
