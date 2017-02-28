NAME = qreg

VERSION_MAJOR = 1
VERSION_MINOR = 2.2
VERSION       = $(VERSION_MAJOR).$(VERSION_MINOR)


LIBRARY = lib$(NAME).so
LIB_VER = $(LIBRARY).$(VERSION)
STATIC  = $(LIBRARY:.so=.a)

EXE = $(NAME).exe


CPP = g++
ASM = as -alhnd

CPPFLAGS = -std=c++14 -fmessage-length=0 -fPIC -pthread -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -Wextra -Wconversion
LDFLAGS  = -shared
ASMFLAGS = -S -fverbose-asm
OBFLAGS  = -std=c++14 -fmessage-length=0 -fPIC -pthread -O3 -g3 -pedantic -Wall -I/usr/include/openbabel-2.0

SOURCES = atom.cpp box.cpp interface.cpp idconverter.cpp logger.cpp config.cpp ifbase.cpp iforc.cpp ifuff.cpp ifopengaff.cpp ifopenuff.cpp ifevb.cpp partitioner.cpp qreg.cpp

EXE_SRC = exe.cpp


HEADERS = $(SOURCES:.cpp=.hpp)

EXE_H   = $(EXE_SRC:.cpp=.hpp)


OBJECTS = $(SOURCES:.cpp=.o)
SS      = $(SOURCES:.cpp=.s)
LST     = $(SOURCES:.cpp=.lst)

EXE_OBJ = $(EXE_SRC:.cpp=.o)
EXE_SS  = $(EXE_SRC:.cpp=.s)
EXE_LST = $(EXE_SRC:.cpp=.lst)



%.o: %.cpp $(SOURCES) $(EXE_SRC) $(HEADERS) $(EXE_H)
	$(CPP) $(CPPFLAGS) -c $< -o $@
	$(CPP) $(CPPFLAGS) $(ASMFLAGS) $< -o $(@:.o=.s)
	$(ASM) $(@:.o=.s) > $(@:.o=.lst) 


all: clean libv static exe


libv: $(OBJECTS)
	$(CPP) $(LDFLAGS) -o $(LIB_VER) $(OBJECTS) $(CPPFLAGS)


static: $(OBJECTS)
	ar -cq $(STATIC) $(OBJECTS)


exe: $(OBJECTS) $(EXE_OBJ)
	$(CPP) $(CPPFLAGS) $(OBJECTS) $(EXE_OBJ) -o $(EXE)


clean:
	rm -rf $(LIBRARY)* $(STATIC) $(EXE) $(OBJECTS) $(EXE_OBJ) $(SS) $(EXE_SS) $(LST) $(EXE_LST) a.out opengaff.o openuff.o opengaff openuff

opengaff:
	rm -rf opengaff.o opengaff
	$(CPP) $(OBFLAGS) -c opengaff.cpp -o opengaff.o
	$(CPP) $(OBFLAGS) opengaff.o -o opengaff -lopenbabel -L/usr/lib/openbabel/2.3.2

openuff:
	rm -rf openuff.o openuff
	$(CPP) $(OBFLAGS) -c openuff.cpp -o openuff.o
	$(CPP) $(OBFLAGS) openuff.o -o openuff -lopenbabel -L/usr/lib/openbabel/2.3.2
