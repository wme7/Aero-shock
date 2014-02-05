'''Save this file as isNumberString.py and complete for the
Is Number String Exercise.'''

def isIntStr(s):
    '''Return True if string s represents an integer in standard notation,
    with digits, possibly following a minus sign.  Otherwise return False.

    >>> isIntStr('-234')
    True
    >>> isIntStr('2-3')
    False
    >>> isIntStr('34a')
    False    
    ''' # code

def isDecimalStr(s):
    '''Return True if string s represents a decimal number in standard
    notation, with digits, possibly with a decimal point, possibly started by
    a minus sign.  Otherwise return False.

    >>> isDecimalStr('-234')
    True
    >>> isDecimalStr('23.456')
    True
    >>> isDecimalStr('3.4.')
    False    
    ''' # code


def main():
    for i in range(5):
        s = input('Enter a string, possibly a number: ')
        print('Is an int:', isIntStr(s))
        print('Is a decimal number:', isDecimalStr(s))

main()    


    
    
