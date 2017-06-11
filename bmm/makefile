#
#Makefile written by Abram Rodgers
#Originally written September 10 2016
#

CC = gcc
CFLAGS = -g -Wall

SRCPATH = src/

OBJPATH = obj/

NUM_SRCS = Interpolator.c MathMethods.c Stack.c

NUM_OBJS = $(NUM_SRCS:.c=.o)

LIN_SRCS = LinAlgTester.c MathMethods.c Stack.c

LIN_OBJS = $(LIN_SRCS:.c=.o)

LP_SRCS = ODE_Solver_LpCircle.c MathMethods.c Stack.c gnuplot_i.c

LP_OBJS = $(LP_SRCS:.c=.o)

TEST_SRCS = testing_ground.c MathMethods.c Stack.c gnuplot_i.c

TEST_OBJS = $(TEST_SRCS:.c=.o)
#This is the name of the executable that showcases the linear algebra functs
LINEAR_ALGEBRA = "LinAlg"

#This is the name of the executable that showcases the numerical methods functs
INTERPOLATOR = "interpolator"

#This is the name of the executable that showcases AB3 ODE Solver and Lp Circles
LP_CIRCLE = "OdeSolver"

#This is the name of the executable that is used for testing new funcitons
TESTING = "testing"

default: brams-math-methods

brams-math-methods: $(NUM_OBJS) $(LIN_OBJS) $(LP_OBJS)
	$(CC) $(CLFAGS) $(NUM_OBJS) -lm -o $(INTERPOLATOR)
	$(CC) $(CLFAGS) $(LIN_OBJS) -lm -o $(LINEAR_ALGEBRA)
	$(CC) $(CLFAGS) $(LP_OBJS) -lm -o $(LP_CIRCLE)
	mkdir -p -- "obj"
	mv *.o $(OBJPATH)
    

%.o : $(SRCPATH)%.c
	$(CC) $(CFLAGS) -lm -c $<

#
#make using the object files stored in the obj folder. crashes if obects not 
#present.
#

tests: $(TEST_OBJS)
	$(CC) $(CLFAGS) $(TEST_OBJS) -lm -o $(TESTING)
	mkdir -p -- "obj"
	mv *.o $(OBJPATH)

objects:
	mv obj/*.o .
	$(CC) $(CLFAGS) $(NUM_OBJS) -lm -o $(INTERPOLATOR)
	$(CC) $(CLFAGS) $(LIN_OBJS) -lm -o $(LINEAR_ALGEBRA)
	mkdir -p -- "obj"
	mv *.o $(OBJPATH)

#
#clean for when you want to remove executables and object files
#
# if obj doesnt exit then make it, thats what the mkdir lines do

clean:
	mkdir -p -- "obj"
	rm -f obj/*
	rm -f $(LINEAR_ALGEBRA)
	rm -f $(INTERPOLATOR)
	rm -f $(LP_CIRCLE)
	rm -f $(TESTING)
	rm -f gnuplot*
	rm -f *.png

