
INC = -I. -I./include
LIB = -lfl
FLG = 
LEX = flex
CXX = g++ -g -O
BIS = bison

SRC = b++.cc a++.cc cards.cc errrec.cc output.cc beam.cc marker.cc part.cc \
      joint.cc material.cc accgrav.cc equilibrium.cc jprim.cc ic.cc pointmass.cc \
      springdamper.cc matrix.cc helpers.cc storage.cc translate.cc \
      sforce.cc gforce.cc vforce.cc vtorque.cc mbdyn.cc mathem.cc \
      manipulate.cc debug.cc variable.cc
OBJ = b++.o  a++.o  cards.o  errrec.o  output.o  beam.o  marker.o  part.o \
      joint.o material.o accgrav.o equilibrium.o jprim.o ic.o pointmass.o \
      springdamper.o matrix.o helpers.o storage.o translate.o sforce.o gforce.o \
      vforce.o vtorque.o mbdyn.o mathem.o manipulate.o debug.o variable.o
      
all: a++

a++: $(OBG)$(OBJ)
	$(CXX) -o a++ $(OBJ) $(LIB)

./src/b++.cc ./src/a++.cc : a++.lex b++.y
	$(LEX) a++.lex; mv lex.yy.c ./src/a++.cc
	$(BIS) b++.y -d -k; mv b++.tab.c ./src/b++.cc

$(OBJ): %.o : ./src/%.cc
	$(CXX) -c $(INC) $(FLG) -o $@ -c $<

depend: a++.cc
	makedepend -- -c ./include/$(INC) -- ./src/$(SRC)

.CLEAN: all clean veryclean

clean:
	rm -f *.o ./src/a++.cc ./src/b++.cc *~ ./src/*~ *.output \
	Makefile.bak *.err *.msg lex.log *.mbd *.ref formula.txt \
	b++.tab.h


# FOR CODE EXCEPTIONS TYPE FLG = -DEXCEPTIONS
# DO NOT DELETE
