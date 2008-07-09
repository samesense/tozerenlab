import sys

def checkStart(sys_args, req_args, examples):
    if not len(sys_args) == len(req_args) + 1:
        print 'ENTER'
        for arg in req_args:
            print '\t', arg
        print 'EXAMPLE'
        print 'python', sys_args[0], str(examples)[1:-1].replace(',', '')
        sys.exit(0)


        
