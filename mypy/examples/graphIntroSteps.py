
input('''Each time you press return, you see the lines displayed,
and the same lines executed.
See the Tutorial for explanation.

''')

from graphics import *
win = GraphWin()
input('''from graphics import *
win = GraphWin()
''')

pt = Point(100, 50)
input('''pt = Point(100, 50)
''')

pt.draw(win)
input('''pt.draw(win)
''')

cir = Circle(pt, 25) 
cir.draw(win) 
input('''cir = Circle(pt, 25) 
cir.draw(win) 
''')

cir.setOutline('red') 
cir.setFill('blue') 
input('''cir.setOutline('red') 
cir.setFill('blue') 
''')

line = Line(pt, Point(150, 100)) 
line.draw(win) 
input('''line = Line(pt, Point(150, 100)) 
line.draw(win) 
''')

rect = Rectangle(Point(20, 10), pt) 
rect.draw(win) 
input('''rect = Rectangle(Point(20, 10), pt) 
rect.draw(win) 
''')

line.move(10, 40)
input('''line.move(10, 40)
''')

win.close()
input('''win.close()

# Done.
# Press return to end.''')
      
