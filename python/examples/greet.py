"""Simple example with Entry objects.
Enter your name, click the mouse, and see greetings.
"""

from graphics import *

def main():
    win = GraphWin("Greeting", 300, 300)
    win.yUp()

    instructions = Text(Point(win.getWidth()/2, 40),
                     "Enter your name.\nThen click the mouse.")
    instructions.draw(win)

    entry1 = Entry(Point(win.getWidth()/2, 200),10)
    entry1.draw(win)

    Text(Point(win.getWidth()/2, 230),'Name:').draw(win) # label for the Entry
    
    win.getMouse()  # To know the user is finished with the text.

    name = entry1.getText() 

    greeting1 = 'Hello, ' + name + '!'
    Text(Point(win.getWidth()/3, 150), greeting1).draw(win)
                     
    greeting2 = 'Bonjour, ' + name + '!'
    Text(Point(2*win.getWidth()/3, 100), greeting2).draw(win)
    
    win.promptClose(instructions)

main()
