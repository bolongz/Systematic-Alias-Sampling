EXE = test
DEFINES = -g -ggdb
CC = g++ -std=c++11 -Wall
OPTFLAG= -O3 
CFLAGS = $(OPTFLAG)
COMPILE = $(CC) $(CFLAGS) -c
OBJ= test.o 

ALL: $(EXE)
$(EXE): $(OBJ) sas.h
	$(CC) $(CFLAGS) -o  $(EXE) $(OBJ)

test.o: test.cpp sas.h
	@$(COMPILE) -o  $@  $< 
clean:
	\rm -f *.o  test
