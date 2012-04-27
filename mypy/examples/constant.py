'''Illustrate a global constant being used inside functions.'''

PI = 3.14159265358979   # global constant -- only place the value of PI is set

def circleArea(radius):
    return PI*radius*radius    # use value of global constant PI

def circleCircumference(radius):
    return 2*PI*radius         # use value of global constant PI

print('circle area with radius 5:', circleArea(5))
print('circumference with radius 5:', circleCircumference(5))
