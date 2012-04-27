#!/usr/bin/env python3

import cgi, os.path

DEFAULT_NAME = "no name"
DEFAULT_ADDR = "nowhere"

def main():
    form = cgi.FieldStorage()   
    
    name = form.getfirst("name", DEFAULT_NAME)
    addr = form.getfirst("addr", DEFAULT_ADDR)
	
    contents = processInput(name, addr)   
    print(contents)
    
def processInput(name, addr):  
    '''Process input parameters and return the final page as a string.'''
    fileName= "namelist.txt"
    if name != DEFAULT_NAME or addr != DEFAULT_ADDR:
        line = "Name: {name}  Address:  {addr}\n".format(**locals())
        append(fileName, line)
    lines = fileLinesToHTMLLines(fileName)
    return fileToStr("nameListTemplate.html").format(**locals())

def append(fileName, s):
    """Append string s to file with name fileName.
    This fails if there are multiple people trying simultaneously.
    """
    fout = open(fileName,'a') # 'a' means append to the end of the file
    fout.write(s)
    fout.close()

def safePlainText(s):
    '''Return string s with reserved html markup characters replaced
    and newlines replaced by <br>.'''
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('\n', '<br>')

def fileLinesToHTMLLines(fileName):
    """Allow lines of the file with name fileName to be embedded in html.
    This fails if there are multiple people trying.
    Alters code with possible tags.
    Returns the empty string if fileName does not exist
    """
    safeLines = list()
    if os.path.exists(fileName): # test if the file exists yet
        lines = fileToStr(fileName).splitlines()
        for line in lines:
            safeLines.append(safePlainText(line))
    return "<br>\n".join(safeLines)
    
# standard functions and code from here on
def fileToStr(fileName): 
    """Return a string containing the contents of the named file."""
    fin = open(fileName); 
    contents = fin.read();  
    fin.close() 
    return contents

try:
    print("Content-type: text/html\n\n")   # say generating html
    main() 
except:
    cgi.print_exception()                 # catch and print errors
