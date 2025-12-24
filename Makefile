# =============================================================================
# Makefile for POLRAD 2.0 (JLab ifarm version)
# Target: polrad20.exe
# Compiler: gfortran
# =============================================================================

# 1. Compiler Setup
FC = gfortran

# 2. CERNLIB Paths (From your ifarm info)
# We use the path you verified: /u/scigroup/cvmfs/scicomp/sw/el9/cernlib/2023
CERN_ROOT = /u/scigroup/cvmfs/scicomp/sw/el9/cernlib/2023

# 3. Compiler Flags
# -O2: Optimization
# -g: Debug info
# -std=legacy: Essential for 1990s Fortran code
# -fno-align-commons: Prevents memory misalignment in COMMON blocks
# -ffixed-line-length-none: (Double 'f') Prevents code cutoff at column 72
# NOTE: Removed -fdefault-real-8 to fix the REAL(16) mismatch errors.
FFLAGS = -O2 -g -std=legacy -fno-align-commons \
         -ffixed-line-length-none -Wno-unused-variable

# 4. Libraries
# We link statically against mathlib, kernlib, and packlib (contains Minuit)
LIBS = -L$(CERN_ROOT)/lib -lmathlib -lkernlib -lpacklib

# 5. Targets
TARGET = polrad20.exe
SRCS = polrad20.f
OBJS = $(SRCS:.f=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "Linking $(TARGET) with CERNLIB..."
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) $(LIBS)
	@echo "Build successful! Run./$(TARGET) to execute."

.f.o:
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
