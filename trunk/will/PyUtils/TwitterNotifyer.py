#import twitter
import threading
import time

import logging

logging.basicConfig(level=logging.DEBUG,format='[%(threadName)-10s] %(message)s %(asctime)s')

def testNotif(waittime,numReps):
    for i in range(numReps):
        time.sleep(waittime)




class NotifContext:
    def __enter__(self):
        #cont = NotifContext()

        self.DoAnnounce(self.DefStartAnnounce)
        self.PaceThread = threading.Thread(name='PaceThread',target=self.PaceAnnounce)
        self.PaceThread.start()
        return self


    def __init__(self,method='twitter',user='MozartComputer',password='easypass'):

        logging.debug('Initializing')
#        self.tweetBird = twitter.Api(username = user, password = password)
        self.user=user
        self.password = password

        self.method = method

        self.updateTimeLimit = 60

        self.UpdatewaitTime = 60
        self.checkTime = 15

        self.currentWaitTime = 0

        self.DefStartAnnounce = 'Starting Computation'
        self.DefContAnnounce = 'Still Computing'
        self.DefFinished = 'Done Computing'

        self.ThisContAnnounce = [self.DefContAnnounce]

        self.stillRunning = threading.Event()

        self.MaxTry = 5
    
    
        
    def __exit__(self,typeIn,valueIn,tracebackIn):
        """
        __exit__
            Called when the main function has ended.  This sets the self.stillRunning Event,
        .join()'s the pacing thread and then sends outputs based on whether the function
        errored or finished properly.
        """
        logging.debug('Made to Exit')

        #Set event so PaceThread will exit next time it checks
        self.stillRunning.set()

        #while waiting, send out final message
        if typeIn == None:
            self.DoAnnounce('Finished Computation')
        else:
            self.DoAnnounce('Something Bad Happened:' + str(typeIn))

        #wait until the thread exits otherwise you get crazyness when python exits
        self.PaceThread.join()

        

        


    def PaceAnnounce(self):
        """
        PaceAnnounce
            Acts as a thread which keeps the time and sends periodic calls of DoAnnounce.
        It continually checks self.stillRunning to see if the function has ended.


        """
        logging.debug('In PaceAnnouce')
        while not(self.stillRunning.isSet()):
            logging.debug('Starting Sleep')
            time.sleep(self.checkTime)
            self.currentWaitTime += self.checkTime
            if (not(self.stillRunning.isSet()) & self.currentWaitTime>self.UpdatewaitTime):
                self.currentWaitTime=0
                self.DoAnnounce(self.ThisContAnnounce.pop())
                self.ThisContAnnounce = [self.DefContAnnounce]


    def UpdateMessage(self,message,now=False):
        """
        UpdateMessage
            This allows the program to send a status updates based on the progress of the function.

        UpdateMessage(message, now=False)
            If now is False then the message will be sent out at the next scheduled update.  If now
        is True then the timer is reset and the message will be sent immediately.


        """

        if now:
            self.currentWaitTime=0
            self.DoAnnounce(message)
        else:
            self.ThisContAnnounce.pop()
            self.ThisContAnnounce = [message]
        
           

    def DoAnnounce(self,message):
        """
        DoAnnouce
            Sends the provided message through the desired channel ... either Twitter or 'Logging'

        DoAnnouce(self,message)
            Sends the provided message.  Will retry self.MaxTry times incase twitter is unresponsive

        Returns:
            True|False
                Returns the success or failure of the message delivery.

        """
        if self.method == 'test':
            logging.info('Sending Tweet: ' + message)
            return True
        elif self.methof == 'twitter':
            for i in xrange(self.MaxTry):
                try:
                    self.tweetBird.PostUpdate(message)
                    return True
                except:
                    time.sleep(15)
            return False
