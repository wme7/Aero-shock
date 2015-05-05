'''prints poorly numbered entries from the list'''

items = ['red', 'orange', 'yellow', 'green']
number = 1
for item in items:
    print(number, item)
    number = 2 # will change to 2 after printing 1
