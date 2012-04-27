''' A neater interactive loop with an agreed upon sentinal value
to end the looping.'''

lines = list()
print('Enter lines of text.')
print('Enter an empty line to quit.')

line = input('Next line: ') # initalize before the loop
while line != '':           # while NOT the termination condition
    lines.append(line)
    line = input('Next line: ')  # !! reset value at end of loop!

print('Your lines were:')
for line in lines:
    print(line)
