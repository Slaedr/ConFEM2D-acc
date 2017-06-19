# Makefile for confem2d-acc
# Aditya Kashi
#
# Make sure environment variables CC and CXX have been set with the compilers to use.
# eg. make CXX=nvcc
# Set DEBUG while/before calling make for a debug build
# Set VIENNACL_DIR as the path where ViennaCL include files are present

NAME := confem2d
PREFIX := build

INCLUDES := ${VIENNACL_DIR}
CXXUSERDEFS := 

#Profile using gprof
PROFILE= #-pg

CXXFLAGS := -Winline -std=c++11

# if DEBUG is not defined, code is compiled in debug mode. Otherwise, optimizations are enabled
ifndef DEBUG

  $(info "Compiling with optimizations")
  CXXFLAGS := ${CFLAGS} -O3

else

  $(info "Compiling debug version")
  CXXFLAGS := ${CFLAGS} -Wall
  CXXUSERDEFS := ${CXXUSERDEFS} -DDEBUG

endif

#LFLAGS := $(LFLAGS) -fPIC
LIBS := 

libsrcs =$(wildcard *.cpp)
libobjst =$(libsrcs:.cpp=.o)
libobjs = $(foreach obj,$(libobjst),$(PREFIX)/$(obj))

$(NAME): $(libobjs)
	$(CXX) $(LFLAGS) -o $(PREFIX)/$(NAME) $(libobjs) $(LIBS) $(PROFILE)

$(PREFIX)/%.o: %.c
	${CXX}  ${CXXFLAGS} -c -o $@ $<  ${INCLUDES} ${LIBS} ${PROFILE} ${LFLAGS} ${CXXUSERDEFS}

#.PHONY : clean
clean:
	rm -f $(PREFIX)/$(NAME) $(libobjs)

