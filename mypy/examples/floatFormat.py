'''Test float formatting.'''

x = 23.457413902458498
print(x, format(x, '.5f'), format(x, '.2f'))
                                  
x = 2.876543
print('longer: {x:.5f}, shorter: {x:.3f}'.format(**locals()))

print('Python approximations to 20 digits for .1, .2, .1 + .2, and .3:')
print(format(.1, '.20f'))
print(format(.2, '.20f'))
print(format(.1 + .2, '.20f'))
print(format(.3, '.20f'))
print()
print('Note that the results for .1 + .2 and for .3 are different!')
