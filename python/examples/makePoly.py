"""Function allowsInteractively chose polygon vertices, ending when click above line.
Illustrates use of while loops, Text, Line, Polygon, getMouse.
"""

from graphics import *

def isBetween(x, end1, end2):
    '''Return True if x is between the ends or equal to either.
    The ends do not need to be in increasing order.'''
    
    return end1 <= x <= end2 or end2 <= x <= end1

def isInside(point, rect):
    '''Return True if the point is inside the Rectangle rect.'''
    
    pt1 = rect.getP1()
    pt2 = rect.getP2()
    return isBetween(point.getX(), pt1.getX(), pt2.getX()) and \
           isBetween(point.getY(), pt1.getY(), pt2.getY())

def polyHere(rect, win):
    '''  Draw a polygon interactively in Rectangle rect, in GraphWin win. 
    Collect mouse clicks inside rect into the vertices of a Polygon,
    and always draw the Polygon created so far.
    When a click goes outside rect, stop and return the final polygon. 
    The Polygon ends up drawn.  The method draws and undraws rect.
    ''' 
    rect.setOutline("red")
    rect.draw(win)
    vertices = list()
    pt = win.getMouse()
    while isInside(pt, rect):
        vertices.append(pt) 
        poly = Polygon(vertices)  
        poly.draw(win)
        pt = win.getMouse()
        poly.undraw() 
    poly.draw(win)
    rect.undraw()
    return poly

def main():
    win = GraphWin('Drawing Polygons', 400, 400)
    win.yUp() 
    
    instructions = Text(Point(win.getWidth()/2, 30),
                        "Click vertices inside the red rectangle."+
                        "\nClick outside the rectangle to stop.")
    instructions.draw(win)

    rect1 = Rectangle(Point(5, 55), Point(200, 120))
    poly1 = polyHere(rect1, win)
    poly1.setFill('green')
    rect2 = Rectangle(Point(210, 50), Point(350, 350))
    poly2 = polyHere(rect2, win)
    poly2.setOutline('orange')

    win.promptClose(instructions)

main()
