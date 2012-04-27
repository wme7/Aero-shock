'''Hello to you!  Illustrates locals() for formating in print.
'''

person = input('Enter your name: ')
greeting = 'Hello, {person}!'.format(**locals())
print(greeting) 
