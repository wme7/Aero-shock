#!/usr/bin/env python3

'''Script to process pizza orders.  The output page depends on whether
the client name is a special OWNER string.  If so, past orders are read and
displayed from an order file.  Otherwise a normal order is taken and recorded.
This version has the data on sizes, toppings and costs hard coded in the form
and a separate copy in the script.  See pizza.cgi for a revision.
'''

import cgi, os.path   

OWNER = 'owner777'  # owner must know this string to use as the name!
ORDER_FILE = 'pizzaOrders.txt'

def main():
    form = cgi.FieldStorage()      
    
    client = form.getfirst('client', '')
    if client == OWNER:
        print(processOwner())
    else:
        size = form.getfirst('size', '') 
        toppings = form.getlist('topping') # note getting *list*
        pastState = form.getfirst('pastState', '')
        print(processCustomer(client, size, toppings, pastState))
    
def processCustomer(client, size, toppings, state):  
    '''Process customer input parameters and return the final page as a string.
    The same basic form is used for the initial page as well as later pages,
    since the customer is encouraged to order multiple pizzas.
    Size data here is a copy of the data on the form.'''

    # pizza size data
    costData = dict() # values for sizes: (baseCost, extra cost per topping)
    costData['small'] = (5.50, .50)
    costData['medium'] = (7.85, .75)
    costData['large'] = (10.50, 1.00)
        
    if not state == 'order':  # initial display of order form
        msg = ''
        invitation = 'Please fill out and submit your order.'
        state = 'order'
    elif not client or not size:  # must have a name and size entered
        msg = '<p><b>You must enter your name and a pizza size with an order!</b></p>'
        invitation = 'Please fill out and submit your order.'
        state = 'order'
    else:  # with a name and size, assume a real order
        if toppings:
            toppingString = ', '.join(toppings)
        else:
            toppingString = 'None'

        (baseCost, toppingCost) = costData[size]
        cost = baseCost + toppingCost * len(toppings)
        msg = '''<p>{client}:  Your ordered a {size} pizza; <br>
        toppings: {toppingString}.<br>
        The cost is ${cost:.2f}.  Please pick it up in 30 minutes.<br>
        Thank you for your order!</p>
        '''.format(**locals())
        line = '{client}; ${cost:.2f}; {size}; {toppings}\n'.format(**locals())
        append(ORDER_FILE, line)
        invitation = 'If you want another pizza, please order again.'
        state = 'order'
    return fileToStr('pizzaOrderTemplate1.html').format(**locals())

def processOwner():
    '''Display all orders from the order file.'''
    text = fileLinesToHTMLLines(ORDER_FILE) or 'No orders presently' #fancy or
    return fileToStr('pizzaReportTemplate.html').format(**locals())
    
# standard functions and code from here on
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


def fileToStr(fileName): 
    """Return a string containing the contents of the named file."""
    fin = open(fileName); 
    contents = fin.read();  
    fin.close() 
    return contents

try:   
    print("Content-type: text/html\n\n")   
    main() 
except:
    cgi.print_exception()               

