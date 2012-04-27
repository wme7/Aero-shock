"""A tiny English to Spanish dictionary is created,
using the Python dictionary type dict.
Then the dictionary is used, briefly.
"""

def createDictionary():
    '''Returns a tiny Spanish dictionary'''
    spanish = dict()
    spanish['hello'] = 'hola'
    spanish['yes'] = 'si'
    spanish['one'] = 'uno'
    spanish['two'] = 'dos'
    spanish['three'] = 'tres'
    spanish['red'] = 'rojo'
    spanish['black'] = 'negro'
    spanish['green'] = 'verde'
    spanish['blue'] = 'azul'
    return spanish

def main():
    dictionary = createDictionary()
    print(dictionary['two'])
    print(dictionary['red'])

main()
