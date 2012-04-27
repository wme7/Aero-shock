'''Program: makeRectangle.py
Demonstrate the use of the clone method in a function makeRect, which
takes a ctakes a corner point and dimensions to construct a Rectangle.
'''

from graphics import *

def makeRect(corner, width, height):
    '''Return a new Rectangle given one corner Point and the dimensions.'''
    
    corner2 = corner.clone()  # !! Note the .clone()  !!
    corner2.move(width, height)
    return Rectangle(corner, corner2)

def main():
    win = GraphWin('Draws a Rectangle', 300, 300)
    win.yUp() 

    rect = makeRect(Point(20, 50), 250, 200)
    rect.draw(win)
    
    win.promptClose(win.getWidth()/2, 20)

main()
