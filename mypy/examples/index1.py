'''first use of a variable for index values'''

s = 'word'
print('The full string is: ', s)
n = len(s)
for i in range(n):
    print()
    print('i =', i)
    print('The letter at index i:', s[i])
    print('The part before index i (if any):', s[:i])
    print('The part before index i+2:', s[:i+2])
    
