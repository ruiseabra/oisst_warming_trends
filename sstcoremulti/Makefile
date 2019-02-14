# Makefile for sstcoremulti

CC	      	= g++
GDB_FLAGS 	= -Wall -ggdb3
CC_OPTS 	  = -O2 -pedantic -std=c++11
LIBS		    = -lm

TARGET		= sstcoremulti
SOURCES		= main.cpp Scanner.cpp ExtractSST.cpp Run.cpp SSTfiles.cpp Str_split.cpp TargetCores.cpp
OBJECTS		= main.o Scanner.o ExtractSST.o Run.o SSTfiles.o Str_split.o TargetCores.o

# Implicit rules

#.SUFFIXES:
.SUFFIXES: .cpp

%.o:	%.cpp
  $(CC) -c $(GDB_FLAGS) $(CC_OPTS) $(INCDIR) $<
  
# Build rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
  $(CC) $(GDB_FLAGS) $(CC_OPTS) $(OBJECTS) $(INCDIR) $(LIBDIR) $(LIBS) -o $(TARGET)

depend:
  @makedepend $(SOURCES) 2> /dev/null

clean:	
  -rm -f *.o *~ *.bak sstcoremulti

# DO NOT DELETE THIS LINE -- make depend depends on it.
