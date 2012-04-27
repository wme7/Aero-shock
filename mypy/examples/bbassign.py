'''
Command line script to organize the files created by the blackboard
assignment download zip *after* it is unzipped into the desired
*parent* directory. The assignment name is passed as a commandline
parameter, or the main method can be called directly with it as
parameter. The script moves assignmentName_user_filename to
assignmentName/user/filename  and converts assignmentName_username.txt
to assignmentName/user/HW-comments.txt (replace slash with backslash
in Windows).'''

import os, os.path

def main(PROJECT_PREFIX):
    PROJECT_PREFIX = PROJECT_PREFIX.strip().replace(' ', '_')
    PROJECT_PREFIX_ = PROJECT_PREFIX + '_'
    file_list = os.listdir(".") # current directory file list
    def filesToMove(name):  # locally defined function, which can read variables
        return  # in the outer function
    file_list = [name for name in file_list # list comprehension for files
                 if name.startswith(PROJECT_PREFIX_)] # satisfying condition
    prefixLen = len(PROJECT_PREFIX_) # include underscore

    if(not os.access(PROJECT_PREFIX, os.F_OK)): # if directory dones not exist
        os.mkdir(PROJECT_PREFIX)                #    create the directory

    for filename in file_list:
        i = prefixLen
        while filename[i].isalnum(): # move past alphanumeric user ID
            i += 1
        username = filename[prefixLen: i]
        newFilename = filename[i+1 :] # +1: skip '_' or '.'; rest is filename
        if newFilename == 'txt': # special file created if comments added
            newFilename = 'HW-comments.txt'
        newDir = os.path.join(PROJECT_PREFIX, username) # join 2 path components
        if(not os.access(newDir, os.F_OK)): # if newDir does not exist
            os.mkdir(newDir)                #    create the directory
        newLocation= os.path.join(newDir, newFilename) # assemble path to file
        os.rename(filename, newLocation) # moves and renames the file

	
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2: # sys.argv lists command line parameters
        print("""
Usage:  bbassign.py assignmentName
     Move all files in the current directory starting with a blank-converted
     assignmentName into directories assignmentName/user.
     It is called from the directory into which assignmentName.zip
     has been unzipped (from a Blackboard assignment).
     The files from the zip archive have the assignment name modified
     if it contained blanks:  the blanks are by replaced by underscores.
     This program makes that substitution automatically.  If there are
     several blanks in a row in the assignment name, I am not sure of the
     conversion algorithm:
     In that case check the prefix on all the zip file contents and make the
     program parameter be the actual assignment prefix used.""")
    else:
        main(sys.argv[1]) # sys.argv[0] is the program being run
    
