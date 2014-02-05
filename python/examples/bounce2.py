'''
Show a ball bouncing off the sides of the window.
User graphically shows the initial direction and speed of the ball.
Control the animation with a while loop.'''

from graphics import *
import time, random
                                                        
def bounceInBox(shape, dx, dy, xLow, xHigh, yLow, yHigh, win):
    ''' Animate a shape moving in jumps (dx, dy), bouncing when
    its center reaches the low and high x and y coordinates.
    The animation stops when the mouse is clicked.'''
    
    delay = .001
    while win.checkMouse() == None:  
        shape.move(dx, dy)
        center = shape.getCenter()
        x = center.getX()
        y = center.getY()
        if x < xLow or x > xHigh:
            dx = -dx
        if y < yLow or y > yHigh:
            dy = -dy            
        time.sleep(delay)

def makeDisk(center, radius, win):
    '''Return a red disk that is drawn in win with given center and radius.'''
    disk = Circle(center, radius)
    disk.setOutline("red")
    disk.setFill("red")
    disk.draw(win)
    return disk

def getShift(point1, point2): # NEW utility function
    '''Returns a tuple (dx, dy) which is the shift from point1 to point2.'''
    dx = point2.getX() - point1.getX()
    dy = point2.getY() - point1.getY()
    return (dx, dy)

def getUserShift(point, prompt, win): #NEW direction selection
    '''Return the change in position from the point to a mouse click in win.
    First display the prompt string under point.'''
    
    text = Text(Point(point.getX(), 60), prompt)
    text.draw(win)
    userPt = win.getMouse()
    text.undraw()
    return getShift(point, userPt)
  
def bounceBall():
    '''Make a ball bounce around the screen.  The user sets the initial speed.'''
    
    win = GraphWin('Ball Bounce', 290, 290)
    win.yUp()

    radius = 10
    xLow = radius # center is separated from the wall by the radius of a bounce
    xHigh = win.getWidth() - radius
    yLow = radius
    yHigh = win.getHeight() - radius

    center = Point(win.getWidth()/2, win.getHeight()/2) #NEW central starting point
    ball = makeDisk(center, radius, win)

    #NEW interactive direction and speed setting
    prompt = '''                            
Click to indicate the direction and
speed of the ball:  The further you
click from the ball, the faster it starts.'''
    (dx, dy) = getUserShift(center, prompt, win)
    scale = 0.01 # to reduce the size of animation steps    
    bounceInBox(ball, dx*scale, dy*scale, xLow, xHigh, yLow, yHigh, win)    
    win.close()
    
bounceBall()
