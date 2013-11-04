# Logger class for managing log levels for errors and info
#
# Written by Kristian Weegink, Nov 2013

class logger:
    def __init__(self,loglevel=0):
        self.loglevel = loglevel
    
    def info(self,text):
        if (self.loglevel > 0):
            print '[INFO] ' + text
    
    def error(self,text):
        if (self.loglevel > 0):
            print '[ERROR] ' + text
