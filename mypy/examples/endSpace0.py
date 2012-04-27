''' Test end keyword parameter for print.'''

def print1():
    print('all', 'on', 'same', 'line')
    print('different line')


def print2():
    '''Illustration of end keyword: same results as print1.'''
    print('all', 'on' , end=' ')
    print('same', end=' ')
    print('line')
    print('different line')

def main():
    print('print1 results:')
    print1()
    print('print2 has the same results:')
    print2()

main()
