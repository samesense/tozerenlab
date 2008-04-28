from threading import *
from Queue import *

finishedDict={}

def worker():
    while not(q.empty()):
        item = q.get()
        if item in finishedDict:
            print 'already processd'
        else:
            finishedDict[item]=item*2
        q.task_done()

q = Queue(-1)

for i in range(0,20,3):
    q.put(i)

for i in range(0,30,2):
    q.put(i)

t = Thread(target=worker)
t.setDaemon(True)
t.start()

t = Thread(target=worker)
t.setDaemon(True)
t.start()


q.join()


