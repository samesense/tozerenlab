from numpy import searchsorted,array
from Queue import *
class PriorityQueue(Queue):
    def _put(self,item):

        if len(item) != 2:
            raise IndexError
        
        thisVal = item[1]
        thisPos = array(self.priority).searchsorted(thisVal)
       
        self.priority.insert(thisPos,thisVal)
        self.queue.insert(thisPos,item[1])

    def _get(self):
        self.priority.pop(0)
        return self.queue.pop(0)

    def _init(self,maxsize):
        self.maxsize = maxsize
        self.queue = list()
        self.priority = list()

