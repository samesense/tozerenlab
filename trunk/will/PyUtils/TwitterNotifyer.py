import twitter
import threading
import time

import logging

logging.basicConfig(level=logging.DEBUG,format='[%(threadName)-10s] %(message)s %(asctime)s')

class NotifContext:
    def __enter__(self):
        #cont = NotifContext()

        self.DoAnnounce(self.DefStartAnnounce)
        self.PaceThread = threading.Thread(name='PaceThread',target=self.PaceAnnounce)
        self.PaceThread.start()
        return self


    def __init__(self,user='MozartComputer',password='easypass'):

        logging.debug('Initializing')
        self.tweetBird = twitter.Api(username = user, password = password)
        self.user=user
        self.password = password

        self.updateTimeLimit = 60

        self.UpdatewaitTime = 600
        self.checkTime = 15

        self.currentWaitTime = 0

        self.DefStartAnnounce = 'Starting Computation'
        self.DefContAnnounce = 'Still Computing'
        self.DefFinished = 'Done Computing'

        self.stillRunning = threading.Event()
    
    
        
    def __exit__(self,typeIn,valueIn,tracebackIn):
        logging.debug('Made to Exit')

        self.stillRunning.set()
        self.PaceThread.join()

        if typeIn == None:
            self.DoAnnounce('Finished Computation')
        else:
            self.DoAnnounce('Something Bad Happened:' + str(typeIn))

        


    def PaceAnnounce(self):
        logging.debug('In PaceAnnouce')
        while not(self.stillRunning.isSet()):
            logging.debug('Starting Sleep')
            time.sleep(self.checkTime)
            self.currentWaitTime += self.checkTime
            if self.currentWaitTime>self.UpdatewaitTime:
                self.currentWaitTime=0
                self.DoAnnounce(self.DefContAnnounce)
           

    def DoAnnounce(self,message):
        logging.debug('Sending Tweet')
        try:
            self.tweetBird.PostUpdate(message)
        except:
            logging.debug('Could Not Send Tweet')
