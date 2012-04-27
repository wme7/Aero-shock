#!/usr/bin/env python3

import cgi   

def main(): 
    form = cgi.FieldStorage()      # standard cgi script lines to here!
    
    ## get data from web form with
    # pythonVariable = form.getfirst(formName, default) for each variable
    ## ...
    contents = processInput(inputVariables) #replace inputVariables with yours
    print(contents)
    
def processInput(inputVariables):  # replace inputVariables with yours
    '''Process input parameters and return the final page as a string.'''
    ## process input, generate output variables ...
    ## ...
    return fileToStr('YOUR_CHOICETemplate.html').format(**locals())
    ## Change YOUR_CHOICE to your chosen name; remember to create the Template!

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


