"""Example with two Entry objects and type conversion.
Version 2 with a function.
Do addition.
"""

from graphics import *

def makeLabeledEntry(entryCenterPt, entryWidth, initialStr, labelText, win):
    '''Return an Entry object with specifed center, width in characters, and
    initial string value.  Also create a static label over it with
    specified text.  Draw everything in the GraphWin win.
    '''
    
    entry = Entry(entryCenterPt, entryWidth)
    entry.setText(initialStr)
    entry.draw(win)

    labelCenter = entryCenterPt.clone()
    labelCenter.move(0, 30)
    Text(labelCenter,labelText).draw(win)
    return entry

        
def main():
    win = GraphWin("Addition", 300, 300)
    win.yUp()

    instructions = Text(Point(win.getWidth()/2, 30),
                     "Enter two numbers.\nThen click the mouse.")
    instructions.draw(win)

    entry1 = makeLabeledEntry(Point(win.getWidth()/2, 250), 25,
                              '0', 'First Number:', win)
    entry2 = makeLabeledEntry(Point(win.getWidth()/2, 180), 25,
                              '0', 'Second Number:', win)
     
    win.getMouse()  # To know the user is finished with the text.

    num1 = int(entry1.getText())
    num2 = int(entry2.getText())
    sum = num1 + num2

    result = "The sum of\n{num1}\nplus\n{num2}\nis {sum}.".format(**locals())
    Text(Point(win.getWidth()/2, 110), result).draw(win)                     
    
    win.promptClose(instructions)

main()
