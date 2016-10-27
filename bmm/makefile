#
#Makefile written by Abram Rodgers
#Originally written September 10 2016
#

CC = gcc
CFLAGS = -g -Wall

SRCPATH = src/

OBJPATH = obj/

SRCS = main.c MathMethods.c

OBJS = $(SRCS:.c=.o)

APP = bmm

default: brams-math-methods

brams-math-methods: $(OBJS)
	$(CC) $(CLFAGS) $(OBJS) -lm -o $(APP)
	mkdir -p -- "obj"
	mv *.o $(OBJPATH)

%.o : $(SRCPATH)%.c
	$(CC) $(CFLAGS) -lm -c $<

#
#clean for when you want to remove executables and object files
#
# if obj doesnt exit then make it, thats what the mkdir lines do

clean:
	mkdir -p -- "obj"
	rm obj/*
	rm bmm

