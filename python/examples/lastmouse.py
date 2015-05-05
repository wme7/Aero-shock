'''A simple test of dealing with the last (previous) mouse click.
The animation stops when the mouse is clicked.
'''

from graphics import *
import time

win = GraphWin()

circle = Circle(Point(100,100), 15)
circle.draw(win)

Text(Point(100, 40), "Click to quit.").draw(win)

while win.checkMouse() == None:
    circle.move(50, 0)
    time.sleep(.5)
    circle.move(-50, 0)
    time.sleep(.5)

win.close()
