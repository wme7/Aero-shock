#! /usr/bin/python3

'''Run a local cgi server from the current directory that treats *.cgi files
as executable python cgi scripts.'''

import http.server, sys, os

class CGIExtHTTPRequestHandler(http.server.CGIHTTPRequestHandler):
    '''This request handler mimics the Loyola server, which looks for CGI files
    to end in '.cgi' and be in any directory as opposed to the CGIHTTPServer
    expectation that the cgi script are of the form /cgi-bin/*.py.'''
    
    def is_python(self, path):
        """Test whether argument path is a Python script: allow .cgi"""
        return path.lower().endswith('.cgi')

    def is_cgi(self):
        '''As on xenon, go by extension only.'''
        base = self.path
        query = ''
        i = base.find('?')
        if i != -1:
            query = base[i:]
            base = base[:i]
        if not base.lower().endswith('.cgi'):
            return False
        [parentDirs, script] = base.rsplit('/', 1)
        self.cgi_info = (parentDirs, script+query)
        return True
   

def run_server():
    dirName = os.getcwd()
    blanks = dirName.count(' ')
    if 0 < blanks:  # server cannot handle blanks in path names
        print("""The path to this directory contains {blanks} space(s):
{dirName}
Either rename directories to remove the blanks or
move this directory to a place with no blanks in the path.
Aborting the local server run!""".format(**locals()))
        input("Press return after reading this message.")
        return
        
    server_addr = ('localhost', 8080)
    cgiServer = http.server.HTTPServer(server_addr, CGIExtHTTPRequestHandler)
    sys.stderr.write('Localhost CGI server started\n.')
    cgiServer.serve_forever()

run_server()
