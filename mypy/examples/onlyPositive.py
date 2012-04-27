'''Use an IF statement inside a FOR loop to select only positive numbers.'''

def printAllPositive(numberList):
    '''Print only the positive numbers in numberList.'''
    for num in numberList:
        if num > 0:
            print(num)

printAllPositive([3, -5, 2, -1, 0, 7])
