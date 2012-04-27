"""Draw Random Circles.  Start and stop with mouse.
"""

from graphics import *
import random, time

def main():
    win = GraphWin("Random Circles", 300, 300)
    text = Text(Point(win.getWidth()/2, 30), "Click to start; click to end")
    text.draw(win)
    win.getMouse()
    text.undraw()
    
    while win.checkMouse() == None:   #NEW
        r = random.randrange(256)
        b = random.randrange(256)
        g = random.randrange(256)
        color = color_rgb(r, g, b)

        radius = random.randrange(3, 40)        
        x = random.randrange(5, 295)
        y = random.randrange(5, 295)
        
        circle = Circle(Point(x,y), radius)
        circle.setFill(color)
        circle.draw(win)
        time.sleep(.05)
        
    win.close()

main()
