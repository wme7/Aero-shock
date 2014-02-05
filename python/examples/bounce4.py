'''Show a ball bouncing off the sides of the window.
User always has mouse control over the direction and speed of the ball.
Version where all conditions that change the direction of the shape
are handled inside one animation loop.'''

from graphics import *
import time, random

def moveInBox(shape, stopHeight, xLow, xHigh, yLow, yHigh, win):
    ''' Animate a shape moving toward any mouse click below stopHeight and
    bouncing when its center reaches the low or high x or y coordinates.
    The animation stops when the mouse is clicked at stopHeight or above.'''
    
    scale = 0.01 
    delay = .001
    dx = 0                      #NEW dx and dy are no longer parameters
    dy = 0                      #NEW
    while True:                 #NEW exit loop at return statement
        center = shape.getCenter()
        x = center.getX()
        y = center.getY()
        isInside = True
        if x < xLow or x > xHigh:
            dx = -dx
            isInside = False
        if y < yLow or y > yHigh:
            dy = -dy            
            isInside = False
        if isInside:
            pt = win.checkMouse()
            if pt != None:              #NEW dealing with mouse click now here
                if pt.getY() < stopHeight: # switch direction
                    (dx, dy) = getShift(center, pt)
                    (dx, dy) = (dx*scale, dy*scale)
                else:                  #NEW exit from depths of the loop
                    return             #NEW
        shape.move(dx, dy)                
        time.sleep(delay)

def makeDisk(center, radius, win):
    '''Return a red disk that is drawn in win with given center and radius.'''
    disk = Circle(center, radius)
    disk.setOutline("red")
    disk.setFill("red")
    disk.draw(win)
    return disk

def getShift(point1, point2): 
    '''Returns a tuple (dx, dy) which is the shift from point1 to point2.'''
    dx = point2.getX() - point1.getX()
    dy = point2.getY() - point1.getY()
    return (dx, dy)
  
def bounceBall():
    '''Make a ball bounce around the screen, and react to mouse clicks.'''
    
    win = GraphWin('Ball Bounce 3', 290, 290)
    win.yUp()

    lineHeight = win.getHeight() - 40
    textHeight = win.getHeight() - 20
    Line(Point(0, lineHeight), Point(win.getWidth(), lineHeight)).draw(win)
    

    radius = 10
    xLow = radius # center is separated from the wall by the radius at a bounce
    xHigh = win.getWidth() - radius
    yLow = radius
    yHigh = lineHeight - radius

    center = Point(win.getWidth()/2, lineHeight/2)
    ball = makeDisk(center, radius, win)
    
    prompt = 'Click above the line to stop\nor below to move toward the click.'
    Text(Point(win.getWidth()/2, textHeight), prompt).draw(win)

    moveInBox(ball, lineHeight, xLow, xHigh, yLow, yHigh, win)    

    win.close()
    
bounceBall()
