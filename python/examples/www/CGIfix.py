#!/usr/bin/env python3
r'''Usage for Mac/Linux/Unix systems:
    CGIfix.py file1 file2 .... # listed files become executable Python scripts
    CGIfix.py              # make all file ending with .cgi executable Python
    CGIfix.py --help       # print this message

For Mac/Unix/Linux systems: Makes Python scripts executable directly 
from the operating system:  In particular make sure files have Unix
'\n' newlines, initial script interpreter line, and are executable.

With no parameters it is specialized for cgi scripts (assumed to be Python).
With any Python program as an explicit parameter it can be used to make the
program executable from the command line without direct reference to the
Python interpreter.

If /usr/bin/python3 does not exist, you can either change SCRIPTLINE
in this source code or make a symbolic link in your system once, with terminal
commands like the ones below for a Mac with Python 3.2 in
/usr/local/bin/python3.2:

cd /usr/bin
sudo ln -s /usr/local/bin/python3.2   python3
'''
import re
import os, os.path
import sys
from stat import *

fixStringRegexp = re.compile('\r')

SCRIPTLINE = '#!/usr/bin/env python3'

def main(fileList=None):
    if fileList:
        if isinstance(fileList, str): 
            fileList = fileList.strip().split()        
    else:
        fileList = [f for f in os.listdir('.') if f.endswith('.cgi')]

    badSymbols = []
    noFind = []
    noRead = []
    noWrite = []
    badFirstLine = []
    ok = []
    for filename in fileList:
        canRead = False
        if filename.find('"') >= 0 or filename.find('%') >= 0 \
           or filename.find('&') >= 0:
            badSymbols.append(filename)
        elif not os.path.exists(filename):
            noFind.append(filename)
        else:
            try:
                fileString = fileToStr(filename)
                canRead = True
            except:
                noRead.append(filename)
        if canRead:
            if '\n' in fileString: # replace all '\r' by '' (delete)
                fixedString = fixStringRegexp.sub('', fileString)
            else: # allow for old apple format with just '\r'
                fixedString = fixStringRegexp.sub('\n', fileString)
            if not fixedString.startswith(SCRIPTLINE):
                badFirstLine.append(filename)
                fixedString = SCRIPTLINE + '\n' + fixedString
                # fancier version might remove an old script line
            canWrite = True
            if fixedString != fileString:
                try:        
                    strToFile(fixedString, filename)
                except:
                    noWrite.append(filename)
                    canWrite = False
            if canWrite:
                mode = getPermission(filename)
                if mode == 'error':
                    noRead.append(filename)
                else:
                    # mode is the permission as bitmask. 
                    readMode = mode & 0o444
                    os.chmod(filename, 0o711 | readMode) # ? can fail?
                    ok.append(filename)
    message = ''
    if badSymbols:
        message = 'Bad names for files: {}\n'.format(badSymbols) 
    if noFind:
        message += "Cannot find: {}\n".format(noFind) 
    if noRead:
        message += "Cannot read: {}\n".format(noRead) 
    if badFirstLine:
        message += "Added first line {}to: {}\n".format(SCRIPTLINE,
                                                        badFirstLine)
    if noWrite:
        message += "Cannot write: {}\n".format(noWrite) 
    if not ok:
        message += "NO files processed.\n"
    elif message:
        message += "The files that were processed are below:\n"
    else:
        message =  "Files Fixed Successfully!\n"
    
    message += ' '.join(ok)        

    print(message)

def getPermission(fileName):
    try:
        st = os.stat(fileName)
    except:
        return "error"
    return st[ST_MODE]

def fileToStr(fileName): 
    """Return a string containing the contents of the named file."""
    fin = open(fileName);
    contents = fin.read();
    fin.close()
    return contents

def strToFile(text, filename):
    """Write a file with the given name and the given text."""
    output = open(filename,"w")
    output.write(text)
    output.close()

if __name__ == '__main__':
    args = sys.argv[1:]
    if args == ['--help']:
        print(__doc__)
    else:
        main(args)    
