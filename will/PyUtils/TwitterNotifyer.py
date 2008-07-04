"""
TwitterNotifyer.py
    A set of methods which can notify a twitter account with the computation
    progress.

import TwitterNotifyer

with TwitterNotify.NotifContext(USER = 'T_User', PASSWORD = 'T_Pass') as notif:
    Do stuff ....
    #you can send updates
    notif.UpdateMessage('Got this Far')


"""

import twitter
import threading
import time
from collections import deque
import logging
import getpass

logging.basicConfig(level = logging.DEBUG,
                    format = '[%(threadName)-10s] %(message)s %(asctime)s')

class NotifContext:
    def __enter__(self):
        """
        __enter__
            Called on entry into the "with" statement.  Starts the pace_thread.
        """

        self.DoAnnounce(self.default_start_ann)
        
        self.pace_thread.start()
        return self


    def __init__(self, METHOD = 'twitter', USER = None, PASSWORD = None):

        if (METHOD == 'twitter') & (USER == None):
            self.user = getpass.getpass(promt = 'Twitter Account Name: ')

        if (METHOD == 'twitter') & (PASSWORD == None):
            self.password = getpass.getpass(prompt = 'Passw for ' + self.user)

        if METHOD == 'twitter':
            self.tweet_bird = twitter.Api(username = self.user,
                                          password = self.password)
        else:
            self.tweet_bird = None

        self.method = METHOD

        if self.method == 'test':
            self.update_time_limit = 60
        else:
            self.update_time_limit = 60*60*2
            
        self.default_start_ann = 'Starting Computation'
        self.default_cont_ann = 'Still Computing'
        self.default_finished = 'Done Computing'

        self.this_cont_ann = deque()
        self.this_cont_ann.append(self.default_cont_ann)

        self.finished_running = threading.Event()

        self.max_try = 5
    
        self.pace_thread = threading.Thread(name = 'pace_thread',
                                            target = self.PaceAnnounce)
        
    def __exit__(self, TYPE_IN, VALUE_IN, TRACEBACK_IN):
        """
        __exit__
            Called when the main function has ended.  This sets the
        self.finished_running Event, .join()'s the pacing thread and then sends
        outputs based on whether the function errored or finished properly.
        """

        #logging.debug('Made to Exit')

        #Set event so pace_thread will exit next time it checks
        self.finished_running.set()

        #while waiting, send out final message
        if TYPE_IN == None:
            self.DoAnnounce('Finished Computation')
        else:
            self.DoAnnounce('Something Bad Happened:' + str(TYPE_IN))

        #wait until the thread exits otherwise you get
        #crazyness when python exits
        self.pace_thread.join()

        

    def Start(self):
        """
        Start
            Calls the .__enter__ method, allows easy usage with pre-"with"
            versions of python
        """
        self.__enter__()

    def Finished(self):
        """
        Finished
            Calls the .__exit__ method, allows easy usage with pre-"with"
            versions of python
        """
        self.__exit__(None, None, None)

    def PaceAnnounce(self):
        """
        PaceAnnounce
            Acts as a thread which keeps the time and sends periodic calls of
            DoAnnounce.  It continually checks self.finished_running to see if
            the function has ended.

        """
        while not(self.finished_running.isSet()):
            #wait until the Event has been triggered or until update limit is up
            self.finished_running.wait(self.update_time_limit)
            
            if not(self.finished_running.isSet()):
                #if the time-limit is up and the Event hasn't been .set() then
                #display an update message
                self.DoAnnounce(self.this_cont_ann.pop())
                if len(self.this_cont_ann) == 0:
                    self.this_cont_ann.append(self.default_cont_ann)


    def UpdateMessage(self, MESSAGE, NOW = False):
        """
        UpdateMessage
            This allows the program to send a status updates based on the
            progress of the function.

        UpdateMessage(message, now=False)
            If now is False then the message will be sent out at the next
            scheduled update.  If now is True then the timer is reset and the
            message will be sent immediately.
        """

        if NOW:
            return self.DoAnnounce(MESSAGE)
        else:
            self.this_cont_ann.append(MESSAGE)
            return False

    def UpdateNow(self, MESSAGE):
        """
        UpdateNow
            A short-cut for .UpdateMessage(MESSAGE, NOW = True)
        """
        return self.DoAnnounce(MESSAGE)
        

    def DoAnnounce(self, MESSAGE):
        """
        DoAnnouce
            Sends the provided message through the desired channel ... either
            Twitter or 'Logging'

        DoAnnouce(self, MESSAGE)
            Sends the provided message.  Will retry self.max_try times incase
            twitter is unresponsive

        Returns:
            True|False
                Returns the success or failure of the message delivery.

        """
        if self.method == 'test':
            print 'Sending Tweet: ' + MESSAGE
            return True
        elif self.method == 'twitter':
            for i in xrange(self.max_try):
                try:
                    self.tweet_bird.PostUpdate(MESSAGE)
                    return True
                except:
                    time.sleep(15)
            return False

def testNotif():
    """
    An example snipet of code.
    """
    cont = NotifContext(METHOD = 'test')
    
    #with NotifContext(method='twitter') as notif:
    for i in range(NUM_REPS):
        #do some code ...
        print i
        if i == UPDATE_ITER:
            #send an update NOW
            NOTIF.UpdateMessage('Update message', now = True)

        elif i == UPDATE_ITER+1:
            #send an update at the next possible occasion
            NOTIF.UpdateMessage('Wait Update message')

        elif i == ERROR_ITER:
            #give an error
            raise IndexError
        else:
            time.sleep(WAIT_TIME)
