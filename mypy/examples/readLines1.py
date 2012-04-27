'''Interactive loop with verbose prompt each time.'''

lines = list()
testAnswer = input('Press y if you want to enter more lines: ')
while testAnswer == 'y':
    line = input('Next line: ')
    lines.append(line)
    testAnswer = input('Press y if you want to enter more lines: ')

print('Your lines were:')
for line in lines:
    print(line)
