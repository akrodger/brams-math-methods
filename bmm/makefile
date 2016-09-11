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
	$(CC) $(CLFAGS) $(OBJS) -o $(APP)
	mv *.o $(OBJPATH)

%.o : $(SRCPATH)%.c
	$(CC) $(CFLAGS) -c $<

#
#clean for when you want to remove executables and object files
#

clean:
	rm obj/*
	rm bmm

