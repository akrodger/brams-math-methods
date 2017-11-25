#include "Stack.h"

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
void StackInit(Stack **stack, long int size)
{
    if(*stack != NULL) //this means the stack is already initialized.
    {
        return;
    }
    *stack = (Stack*) calloc(1, sizeof(Stack));
    (*stack)->currentItemIndex = -1;
    (*stack)->MAX_SIZE = size;
    (*stack)->stackItems = (long int*) calloc(size, sizeof(long int));

    return;
}

/**
 * Frees the memoy used up by a Stack structure. Uses free() frunction 
 * to do this.
 *
 * @param *stack the stack to free.
 */
void StackFree(Stack *stack)
{
    free(stack->stackItems);
    free(stack);
}

/**
 * This function pushes a value onto the stack. 
 *
 * @param *stack the Stack to push an int onto.
 *
 * @param value the Integer to place onto the stack
 *
 * return STANDARD_ERROR if stack not initialize or stack is full
 *                      otherwise return SUCCESS.
 */
int StackPush(Stack *stack, long int value)
{
    if ((stack == NULL) || (stack->currentItemIndex >= stack->MAX_SIZE))
    {
        //this if is the logical equivalent of the documentation above
        return STANDARD_ERROR; //it messed up
    } else {
        stack->currentItemIndex += 1; //increment current item
        stack->stackItems[stack->currentItemIndex] = value; //set value
        return SUCCESS;
    }
}

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
int StackPop(Stack *stack, long int *value)
{
    if ((stack == NULL) || (stack->currentItemIndex <= -1)) 
    {
        //this if is the logical equivalent of the documentation above
        return STANDARD_ERROR; //it messed up
    } else {
        *value = stack->stackItems[stack->currentItemIndex]; //set value
        stack->currentItemIndex -= 1; //decrement current item
        return SUCCESS;
    }
}

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
 * return STANDARD_ERROR if not initialized or if not enough room in *dest
 */
int StackCopy(Stack *src, Stack *dest) 
{
    int i = 0;
    if( (src == NULL) || (dest == NULL) || 
        (StackGetSize(src) > (dest->MAX_SIZE - StackGetSize(dest)) ) )
    {
        return STANDARD_ERROR;
    }
    for (i = 0; i < StackGetSize(src); i++)
    {
        StackPush(dest, src->stackItems[i]);
    }
    return SUCCESS;
}

/**
 * This function checks for whether the stack is empty or not. 
 *
 * @param *stack the stack to check
 *
 * return 0 if stack is NULL/uninitialized or has 0 elements
 *           otherwise, return 1
 */
int StackIsEmpty(const Stack *stack)
{
    if ((stack = NULL) || (stack->currentItemIndex >= -1)) {
        //this if is the logical equivalent of the documentation above
        return 0;
    } else {
        return 1;
    }
}

/**
 * This function checks for whether the stack is full or not. 
 *
 * @param *stack the stack to check
 *
 * return 0 if stack is NULL/uninitialized or has under max number of elements
 *           otherwise, return 1
 */
int StackIsFull(const Stack *stack)
{
    if ((stack = NULL) || (stack->currentItemIndex != (stack->MAX_SIZE - 1)) )
    {
        return 0;
    }else{
        return 1;
    }
}

/**
 * Returns the current size of the stack in terms of how many active elements
 * are in it. 
 *
 * @param *stack the stack to check
 *
 * @return the number of elements in the stack. 
 *          if not initialized, return SIZE_ERROR. (see Stack.h for error codes)      
 */
long int StackGetSize(const Stack *stack)
{
    if(stack == NULL){
        return SIZE_ERROR;
    }else{
        return (stack->currentItemIndex + 1);
    }

}
