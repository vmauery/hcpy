'''
Copyright (c) 2009, Don Peterson
Copyright (c) 2011, Vernon Mauery
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.
* The names of the contributors may not be used to endorse or
promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''


from sys import stdout, stderr
import time

nl = "\n"

class Display(object):
    '''Derive from this object if you'd like to use other methods of
    display.  This base object just uses stdout and stderr and thus
    should work with any console.
    '''
    def __init__(self, out_stream=stdout, err_stream=stderr):
        self.out = out_stream.write
        self.error = err_stream.write
        self.enabled = True
        self.streams = []

    def msg(self, string, suppress_nl=False):
        '''Normal message string to the user.
        '''
        if self.enabled == True:
            if suppress_nl:
                self.out(string)
            else:
                self.out(string + nl)
            self.log(string, suppress_nl)

    def err(self, string, suppress_nl=False):
        '''Error message string to the user.
        '''
        if suppress_nl:
            self.error(string)
        else:
            self.error(string + nl)
        self.log(string, suppress_nl)

    def on(self):
        'Enable normal messages.'
        self.enabled = True

    def off(self):
        '''Disable normal messages.  Error messages can't be turned off.'''
        self.enabled = False

    def logon(self, stream):
        '''Add a stream to send output to.'''
        self.streams.append(stream)
        stream.write("<< On " + time.asctime(time.localtime()) + ">>" + nl)

    def log(self, string, suppress_nl=False):
        for stream in self.streams:
            try:
                if suppress_nl:
                    stream.write(string)
                else:
                    stream.write(string + nl)
            except:  # Ignore bad streams
                pass

    def logoff(self):
        '''Turn off all streams.'''
        self.log("<< Off " + time.asctime(time.localtime()) + ">>")
        self.streams = []

    def __del__(self):
        '''logoff in case we're exitting because of an exception.'''
        self.logoff()
