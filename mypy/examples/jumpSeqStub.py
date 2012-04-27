'''Save as jumpSeq.py and complete it for the Strange Sequence Exercise.'''

def jump(n):
    '''Assume n is an integer.  If n is even, return n//2.
    Otherwise (n is odd) return 3*n + 1.
    '''
    # write or copy from jumpFunc.py!

def printJumps(n):
    ''' Iterate the jump function, starting with the positive integer value n,
    stopping when the latest value is 1.  Print each value:

    >>> printJumps(1)
    1
    >>> printJumps(5)
    5
    16
    8
    4
    2
    1
    '''
    # code before or after coding listJumps
    

def listJumps(n):
    '''Starting from positive integer n, generate the same iterates of the
    jump function as in printJumps(n), but put them in a list and return them:

    >>> listJumps(1)
    [1]
    >>> listJumps(5)
    [5, 16, 8, 4, 2, 1]
    '''
    # code
    
    
def main():
    n = int(input('Enter a positive integer: '))
    printJumps(n)
    print(listJumps(n)) # prints None if yo have not coded listJump

main()

