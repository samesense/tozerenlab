#-----------------------------------------
#  Author:     Perry Evans
#              evansjp@mail.med.upenn.edu
#  2008
#-----------------------------------------
"""
Functions in this module aid in writing
small scripts
"""
import sys

def checkstart(sys_args, req_args, examples):
    """Place this at the start of cmd line scripts to check for
       the correct # of script argumtns.  If you do not recieve
       the correct arguments, the script exits and suggests 
       which arguments should be used.

    @param sys_args: sys.argv
    @param req_args: list of the arguments you want for the 
                     script
    @param examples: list of examples arguments for the script
    """

    if not len(sys_args) == len(req_args) + 1:
        print 'ENTER'
        for arg in req_args:
            print '\t', arg
        print 'EXAMPLE'
        print 'python', sys_args[0], str(examples)[1:-1].replace(',', '')
        sys.exit(0)


        
