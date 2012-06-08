'''
Show a ball bouncing off the sides of the window.
'''

from graphics import *
import time, random

def makeDisk(center, radius, win, color='red'):
    '''return a red disk that is drawn in win with given center and radius.'''
    disk = Circle(center, radius)
    disk.setOutline(color)
    disk.setFill(color)
    disk.draw(win)
    return disk

def main():
    winWidth = 60
    winHeight = 60
    win = GraphWin('Ball Bounce', winWidth, winHeight)
    win.yUp()

    radius = 10
    xLow = radius # center is separated from the wall by the radius at a bounce
    xHigh = winWidth - radius
    yLow = radius
    yHigh = winHeight - radius

    center = Point(xLow+1, yHigh//3)
    makeDisk(center, radius, win)
    makeDisk(center, 2, win, 'black')
    Line(Point(0, yHigh//3), center).draw(win)
    
    win.getMouse()
    win.close()
    
main()
