'''Create an html file with  user input (a name) embedded,
and call the default web browser to display the file.'''

# NEW more appropriate name, now that it is a format string
pageTemplate = '''
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
  <meta content="text/html; charset=ISO-8859-1"
 http-equiv="content-type">
  <title>Hello</title>
</head>
<body>
Hello, {person}!
</body>
</html>''' # NEW note '{person}' two lines up

def main():    # NEW
    person = input("Enter a name: ")  
    contents = pageTemplate.format(**locals())   
    browseLocal(contents, 'tempname.html') 

def strToFile(text, filename):
    """Write a file with the given name and the given text."""
    output = open(filename,"w")
    output.write(text)
    output.close()

def browseLocal(webpageText, filename):
    """Start your webbrowser on a local file containing the text."""
    strToFile(webpageText, filename)
    import webbrowser
    webbrowser.open(filename)
main()
