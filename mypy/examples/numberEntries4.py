''' use a function to number the entries in any list'''

def numberList(items):
    '''Print each item in a list items, numbered in order.'''
    number = 1
    for item in items:
        print(number, item)
        number = number + 1

def main():
    numberList(['red', 'orange', 'yellow', 'green'])
    print()
    numberList(['apples', 'pears', 'bananas'])

main()
