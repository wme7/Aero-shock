def sumList(nums):
    '''Return the sum of the numbers in the list nums.'''
    sum = 0
    for num in nums:
        sum = sum + num
    return sum

print(sumList([5, 2, 4, 7]))
