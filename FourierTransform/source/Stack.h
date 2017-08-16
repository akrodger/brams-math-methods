#include <stdint.h>
#include <stdlib.h>

#ifndef STACK_H
#define STACK_H

#ifndef NULL
#define NULL 0
#endif

//Some standard error codes
enum {
    SIZE_ERROR = -1,
    STANDARD_ERROR,
    SUCCESS
};



/**
 *  This is formalization of an integer Stack in C. declare as many as youwant,
 *  go wild! Proper usage of this struct is to declare as a NULL pointer of type 
 *  Stack*.Then pass that pointer into StackInit() for memory allocation.
 *  To get rid of a stack and free up the memory,  StackFree() on that stack.
 *
 * So:
 *  - MAX_SIZE: a long int saying the maximum number of members in the stack.
 *              we use long int because 
 *  - stackItems: Contains all the ints that are in the stack in ascending
 *                order.
 *  - currentItemIndex: Contains the index in stackItems of the top of the
 *                      stack.
 */
typedef struct{
    long int MAX_SIZE;
	long int *stackItems;
	long int currentItemIndex;
}Stack;

/**
 * This function handles initializing a dynamic stack of integers.
 * We consider a NULL pointer to be an unitialized Stack. This function handles
 * all the memory allocation for a given stack. The size argument is how long
 * you want the Stack to be. A stack cannot be longer than this.
 *
 * Usage:
 * 
 * Stack *s = NULL; //declare your stack and initialize it as NULL
 * StackInit(&s, <long int>); //initialize it. <long int> is the max size of s
 *
 * @param **stack the address of a pointer which was initialized to NULL.
 *                a Stack structure is loaded into the pointer at this address..
 *
 * @param size the maximum length of the Stack you want to initialize.
 */
void StackInit(Stack **stack, long int size);

/**
 * Frees the memoy used up by a Stack structure. Uses free() frunction 
 * to do this.
 *
 * @param *stack the stack to free.
 */
void StackFree(Stack *stack);

/**
 * This function pushes a value onto the stack. 
 *
 * @param *stack the Stack to push something onto.
 *
 * @param value the Integer to place onto the stack
 *
 * return STANDARD_ERROR if stack not initialize or stack is full
 *                      otherwise return SUCCESS.
 */
int StackPush(Stack *stack, long int value);

/**
 * This function pops a value off of the stack. 
 *
 * @param *stack the Stack to pop an int off of.
 *
 * @param *value this pointer is used to return the popped value into an integer 
 *
 * return STANDARD_ERROR if stack not initialized or stack is full
 *                      otherwise return SUCCESS.
 */
int StackPop(Stack *stack, long int *value);

/**  
 *  Copies the data of one stack into another..reads from *src 
 *  and pushes into *dest
 *  
 *  This uses a little bit cheating, treating *src as an array, but this is what
 *  must be done if you want to copy easily.
 *
 *  @param *src stack to copy from
 *  @param *dest stack to copy into
 *
 *  return STANDARD_ERROR if either not initialized or 
 *                        if not enough room in *dest
 *                        return SUCCESS otherwise
 */
int StackCopy(Stack *src, Stack *dest);

/**
 * This function checks for whether the stack is empty or not. 
 *
 * @param *stack the stack to check
 *
 * return 0 if stack is NULL/uninitialized or has 0 elements
 *           otherwise, return 1
 */
int StackIsEmpty(const Stack *stack);

/**
 * This function checks for whether the stack is full or not. 
 *
 * @param *stack the stack to check
 *
 * return 0 if stack is NULL/uninitialized or has under max number of elements
 *           otherwise, return 1
 */
int StackIsFull(const Stack *stack);

/**
 * Returns the current size of the stack in terms of how many active elements
 * are in it. 
 *
 * @param *stack the stack to check
 *
 * @return the number of elements in the stack. 
 *          if not initialized, return SIZE_ERROR. (see Stack.h for error codes)      
 */
long int StackGetSize(const Stack *stack);

#endif // STACK_H
