'''Program: makeRectBad.py
Attempt a function makeRect (incorrectly),
which takes a takes a corner point and dimensions to construct a Rectangle.
'''

from graphics import *

def makeRect(corner, width, height):  # Incorrect!
    '''Return a new Rectangle given one corner Point and the dimensions.'''
    corner2 = corner
    corner2.move(width, height)
    return Rectangle(corner, corner2)

def main():
    win = GraphWin('Draw a Rectangle (NOT!)', 300, 300)
    win.yUp() 

    rect = makeRect(Point(20, 50), 250, 200)
    rect.draw(win)
    
    win.promptClose(win.getWidth()/2, 20)

main()

