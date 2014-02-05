'''Program: triangle.py or triangle.pyw (best name for Windows)
Interactive graphics program to draw a triangle,
with prompts in a Text object and feedback via mouse clicks.
'''

from graphics import *

def main():
    win = GraphWin('Draw a Triangle', 350, 350)
    win.yUp() # right side up coordinates
    win.setBackground('yellow')
    message = Text(Point(win.getWidth()/2, 30), 'Click on three points')
    message.setTextColor('red')
    message.setStyle('italic')
    message.setSize(20)
    message.draw(win)

    # Get and draw three vertices of triangle
    p1 = win.getMouse()
    p1.draw(win)
    p2 = win.getMouse()
    p2.draw(win)
    p3 = win.getMouse()
    p3.draw(win)
    vertices = [p1, p2, p3]

    # Use Polygon object to draw the triangle
    triangle = Polygon(vertices)
    triangle.setFill('gray')
    triangle.setOutline('cyan')
    triangle.setWidth(4)  # width of boundary line
    triangle.draw(win)

    message.setText('Click anywhere to quit') # change text message
    win.getMouse()
    win.close() 

main()
