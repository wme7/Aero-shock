'''Make a choice of colors via mouse clicks in Rectangles --
A demonstration of Boolean operators and Boolean functions.'''

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

def makeColoredRect(corner, width, height, color, win):
    ''' Return a Rectangle drawn in win with the upper left corner
    and color specified.'''

    corner2 = corner.clone()  
    corner2.move(width, -height)
    rect = Rectangle(corner, corner2)
    rect.setFill(color)
    rect.draw(win)
    return rect

def main():
    win = GraphWin('pick Colors', 400, 400)
    win.yUp() # right side up coordinates

    redButton = makeColoredRect(Point(310, 350), 80, 30, 'red', win)
    yellowButton = makeColoredRect(Point(310, 310), 80, 30, 'yellow', win)
    blueButton = makeColoredRect(Point(310, 270), 80, 30, 'blue', win)

    house = makeColoredRect(Point(60, 200), 180, 150, 'gray', win)
    door = makeColoredRect(Point(90, 150), 40, 100, 'white', win)
    roof = Polygon(Point(50, 200), Point(250, 200), Point(150, 300))
    roof.setFill('black')
    roof.draw(win)
    
    msg = Text(Point(win.getWidth()/2, 375),'Click to choose a house color.')
    msg.draw(win)
    pt = win.getMouse()

    if isInside(pt, redButton):
        color = 'red'
    elif isInside(pt, yellowButton):
        color = 'yellow'
    elif isInside(pt, blueButton):
        color = 'blue'
    else :
        color = 'white'
    house.setFill(color)
    
    msg.setText('Click to choose a door color.')
    pt = win.getMouse()

    if isInside(pt, redButton):
        color = 'red'
    elif isInside(pt, yellowButton):
        color = 'yellow'
    elif isInside(pt, blueButton):
        color = 'blue'
    else :
        color = 'white'
    door.setFill(color)

    win.promptClose(msg)

main()
