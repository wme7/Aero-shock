#!/usr/bin/python3

import cgi   # NEW

def main(): # NEW except for the call to processInput
    form = cgi.FieldStorage()      # standard cgi script lines to here!
    
    # use format of next two lines with YOUR names and default data
    numStr1 = form.getfirst("x", "0") # get the form value associated with form
                                   # name 'x'.  Use default "0" if there is none. 
    numStr2 = form.getfirst("y", "0") # similarly for name 'y'
    contents = processInput(numStr1, numStr2)   # process input into a page
    print(contents)
    
def processInput(numStr1, numStr2):  
    '''Process input parameters and return the final page as a string.'''
    num1 = int(numStr1) # transform input to output data
    num2 = int(numStr2)
    total = num1+num2
    return fileToStr('additionTemplate.html').format(**locals())

# standard code for future cgi scripts from here on
def fileToStr(fileName): 
    """Return a string containing the contents of the named file."""
    fin = open(fileName); 
    contents = fin.read();  
    fin.close() 
    return contents

try:   # NEW
    print("Content-type: text/html\n\n")   # say generating html
    main() 
except:
    cgi.print_exception()                 # catch and print errors


