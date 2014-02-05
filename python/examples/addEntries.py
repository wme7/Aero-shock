"""Example with two Entry objects and type conversion.
Do addition.
"""

from graphics import *

def main():
    win = GraphWin("Addition", 300, 300)
    win.yUp()

    instructions = Text(Point(win.getWidth()/2, 30),
                     "Enter two numbers.\nThen click the mouse.")
    instructions.draw(win)

    entry1 = Entry(Point(win.getWidth()/2, 250),25)
    entry1.setText('0')
    entry1.draw(win)

    Text(Point(win.getWidth()/2, 280),'First Number:').draw(win)                   
    
    entry2 = Entry(Point(win.getWidth()/2, 180),25)
    entry2.setText('0')
    entry2.draw(win)

    Text(Point(win.getWidth()/2, 210),'Second Number:').draw(win)
    
    win.getMouse()  # To know the user is finished with the text.

    numStr1 = entry1.getText()
    num1 = int(numStr1)

    numStr2 = entry2.getText()
    num2 = int(numStr2)
    sum = num1 + num2

    result = "The sum of\n{num1}\nplus\n{num2}\nis {sum}.".format(**locals())
    Text(Point(win.getWidth()/2, 110), result).draw(win)                     
    
    win.promptClose(instructions)

main()
