'''
Show a ball bouncing off the sides of the window.
User always has mouse control over the direction and speed of the ball.
'''

from graphics import *
import time, random

def bounceInBox(shape, dx, dy, xLow, xHigh, yLow, yHigh, win):
    ''' Animate a shape moving in jumps (dx, dy), bouncing when
    its center reaches the low and high x and y coordinates.
    The animation stops when the mouse is clicked, and the
    last mouse click is returned.'''
    
    delay = .001
    pt = None                        #NEW
    while  pt == None:               #NEW
        shape.move(dx, dy)
        center = shape.getCenter()
        x = center.getX()
        y = center.getY()
        isInside = True              #NEW
        if x < xLow or x > xHigh:
            dx = -dx
            isInside = False         #NEW
        if y < yLow or y > yHigh:
            dy = -dy
            isInside = False         #NEW
        time.sleep(delay)
        if isInside: # NEW  don't mess with dx, dy when outside  
            pt = win.checkMouse()    #NEW
    return pt                        #NEW

def moveInBox(shape, stopHeight, xLow, xHigh, yLow, yHigh, win): #NEW
    '''Shape bounces in win so its center stays within the low and high
    x and y coordinates, and changes direction based on mouse clicks,
    terminating when there is a click above stopHeight.'''

    scale = 0.01 
    pt = shape.getCenter() # starts motionless
    while pt.getY() < stopHeight:
       (dx, dy) = getShift(shape.getCenter(), pt)
       pt = bounceInBox(shape, dx*scale, dy*scale,
                        xLow, xHigh, yLow, yHigh, win)
      
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

    #NEW to mark and label the area where a click stops the program
    lineHeight = win.getHeight() - 40
    textHeight = win.getHeight() - 20
    Line(Point(0, lineHeight), Point(win.getWidth(), lineHeight)).draw(win)
    prompt = 'Click above the line to stop\nor below to move toward the click.'
    Text(Point(win.getWidth()/2, textHeight), prompt).draw(win)

    radius = 10
    xLow = radius # center is separated from the wall by the radius at a bounce
    xHigh = win.getWidth() - radius
    yLow = radius
    yHigh = lineHeight - radius  #NEW lower top to bouncing limits

    center = Point(win.getWidth()/2, lineHeight/2)
    ball = makeDisk(center, radius, win)

    moveInBox(ball, lineHeight, xLow, xHigh, yLow, yHigh, win) #NEW

    win.close()
    
bounceBall()
