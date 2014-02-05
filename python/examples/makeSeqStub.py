'''Stub for first strange sequence exercise.'''


def nextInSeq(x):
    '''Given an integer x,
    return x/2 if x is even, and 3x+1 if n is odd.
    >>> nextInSeq(5)
    16
    >>> nextInSeq(16)
    8
    '''
    #code here for makeSeq1!

def makeSeq(x, n):
    '''Print n elements of the sequence starting with integer x,
    where each later element in the sequence is nextInSeq applied to the
    previous element.
    >>> makeSeq(5, 3)
    5
    16
    8
    '''
    # Code here for makeSeq1!

def seqUntilOne(x):
    '''Start printing elements of the sequence starting with integer x,
    where each later element in the sequence is nextInSeq applied to the
    previous element.  Stop when 1 is printed.
    Return the total number of terms.
    >>> tot = seqUntilOne(5)
    5
    16
    8
    4
    2
    1
    >>> tot
    6
    '''
    # Code here for the second exercise makeSeq2!

main():
    x = input('Enter a positive integer: ')
    n = input('Enter a number of terms')    #   comment out for makeSeq2
    makeSeq(x, n)                           #   comment out for makeSeq2
##    tot = seqUntilOne(x)                  # uncomment for makeSeq2
##    print('Total terms:', tot)            # uncomment for makeSeq2

main()
