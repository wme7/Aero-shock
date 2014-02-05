'''Triggers error inside nested function call'''

def happyBirthday(person):
    print("Happy Birthday to you!")
    print("Happy Birthday to you!")
    print("Happy Birthday, dear " + person + ".")
    print("Happy Birthday to you!")

def main():
    happyBirthday('Emily')
    happyBirthday('Andre')
    # next line causes error, illustrates traceback
    happyBirthday(2) 

main()
