""" Simple if statement example
mentions extra charge for luggage over 50 lb.
"""

def main():    
    weight = float(input("How many pounds does your suitcase weigh? "))
    if weight > 50:
        print("There is a $25 charge for luggage that heavy.")
    print("Thank you for your business.")

main()
