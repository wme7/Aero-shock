'''Save this file as safeNumberInput.py and complete for the
SafeNumber Input Exercise.'''

def safeWholeNumber(prompt):
    '''Prompt the user to enter a whole number. As long as the input
    is not a legal whole number string, point out the error and prompt again.
    Return the int value of the first string that represents a legal
    whole number.  The shell illustrations show prompts and inputs, and then
    the value finally returned:

    >>> safeWholeNumber('Enter the number of people: ')
    Enter the number of people: xy
    Error!  Enter the number of people: 14a
    Error!  Enter the number of people: 14
    14
    >>> safeWholeNumber('Enter the number of dogs: ')
    Enter the number of dogs: 5
    5
    '''

    return 0 # change: this is just to make the stub runnable


def safeInt(prompt):
    '''Prompt the user to enter a integer, and return the first legal int value.
    This is the same idea as safeWholeNumber, except the integer string may
    be negative.
    '''

    return 0 # change: this is just to make the stub runnable

def safeDecimal(prompt):
    '''Prompt the user to enter a decimal, and return the first legal float
    value. This is the same idea as safeWholeNumber, except the
    decimal string may start with a '-' and contain one '.'
    '''

    return 0.0 # change: this is just to make the stub runnable

def main():
    print('Mess this up to test!')

    x = safeWholeNumber('Enter whole number: ')
    y = safeWholeNumber('Enter whole number: ')
    print('The sum is', x+y)

    x = safeInt('Enter integer: ')
    y = safeInt('Enter integer: ')
    print('The sum is', x+y) # should be right after you fill in safeInt

    x = safeDecimal('Enter decimal number: ')
    y = safeDecimal('Enter decimal number: ')
    print('The sum is', x+y) # should be right after you fill in safeDecimal

main()    


    
    
