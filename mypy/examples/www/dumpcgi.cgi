#!/usr/bin/env python3

import cgi

def safePlainText(s):
    '''Return string s with reserved html markup characters replaced
    and newlines replaced by <br>.'''
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('\n', '<br>')

def dumpFormRaw(form):
    '''Return html fragment with key, value data and line breaks.'''
    lines = []
    if form:
        lines.append('Raw CGI Data:<br>')  # <br> forces a newline in html
        for name in form:
            val = form.getlist(name) # returns values for key as a list
            if len(val) == 1: # generally there is just one value per key
                val = val[0]  #   if so I choose to display it without the list
            line = '{name}: {val}'.format(**locals())# show single value or list
            lines.append(safePlainText(line)) # replace special characters
    else:
        lines.append('No CGI data received')
    return ' <br>\n'.join(lines)  # need <br> to force new line in HTML

def dumpFormPage(form):
    '''Return a formatted HTML page with all the data dumped from the form.'''

    return '''<html><title>Raw CGI Data</title><body>
                 {0}
              </body></html>'''.format(dumpFormRaw(form))
    
def main():
    form = cgi.FieldStorage()
    # lines above are totally standard
    print(dumpFormPage(form)) # easy check of all your HTML form data

# this is the immediately executed code below - copy it exactly    
try:
   print("Content-type: text/html\n\n") 
   main() 
except:
   cgi.print_exception()
