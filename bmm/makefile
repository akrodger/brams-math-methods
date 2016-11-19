#
#Makefile written by Abram Rodgers
#Originally written September 10 2016
#

CC = gcc
CFLAGS = -g -Wall

SRCPATH = src/

OBJPATH = obj/

NUM_SRCS = InterpolateTester.c MathMethods.c Stack.c

NUM_OBJS = $(NUM_SRCS:.c=.o)

LIN_SRCS = LinAlgTester.c MathMethods.c Stack.c

LIN_OBJS = $(LIN_SRCS:.c=.o)

#This is the name of the executable that showcases the linear algebra functs
LINEAR_ALGEBRA = LinAlg

#This is the name of the executable that showcases the numerical methods functs
NUMERICS = "Numer"

default: brams-math-methods

brams-math-methods: $(NUM_OBJS) $(LIN_OBJS)
	$(CC) $(CLFAGS) $(NUM_OBJS) -lm -o $(NUMERICS)
	$(CC) $(CLFAGS) $(LIN_OBJS) -lm -o $(LINEAR_ALGEBRA)
	mkdir -p -- "obj"
	mv *.o $(OBJPATH)
    

%.o : $(SRCPATH)%.c
	$(CC) $(CFLAGS) -lm -c $<

#
#make using the object files stored in the obj folder. crashes if obects not 
#present.
#

objects:
	mv obj/*.o .
	$(CC) $(CLFAGS) $(NUM_OBJS) -lm -o $(NUMERICS)
	$(CC) $(CLFAGS) $(LIN_OBJS) -lm -o $(LINEAR_ALGEBRA)
	mkdir -p -- "obj"
	mv *.o $(OBJPATH)

#
#clean for when you want to remove executables and object files
#
# if obj doesnt exit then make it, thats what the mkdir lines do

clean:
	mkdir -p -- "obj"
	rm obj/*
	rm $(LINEAR_ALGEBRA)
	rm $(NUMERICS)


    
