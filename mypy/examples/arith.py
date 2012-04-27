'''Fancier format string example with 
parameter identification numbers
-- useful when some parameters are used several times.'''

x = 20
y = 30
formatStr = '{0} + {1} = {2}; {0} * {1} = {3}.'
equations = formatStr.format(x, y, x+y, x*y)
print(equations)
