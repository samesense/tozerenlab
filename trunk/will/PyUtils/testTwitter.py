from __future__ import with_statement
import TwitterNotifyer as Notif
import nose.tools
import threading


def testNotif():
    """
    Test basic importing
    """
    context = Notif.NotifContext(METHOD = 'test')
    nose.tools.assert_false(context == None)

def testTwitterAPI():
    """
    Test twitter API
    """
    context = Notif.NotifContext(METHOD = 'twitter',
                                 USER = 'testuser',
                                 PASSWORD = 'testpass')

    nose.tools.assert_not_equal(context, None)


def testDoAnnnounce():
    """
    Test that DoAnnounce returns a value
    """
    context = Notif.NotifContext(METHOD = 'test')
    nose.tools.assert_true(context.DoAnnounce('Just a Test'))

def testWithContext():
    """
    Test "with" and the creation and destruction of the pace_thread
    """
    with Notif.NotifContext(METHOD = 'test') as context:
        nose.tools.assert_true(context.pace_thread.isAlive(),
                               'Did not create PaceThread')
        nose.tools.assert_false(context.finished_running.isSet(),
                                'Event not set properly')

    nose.tools.assert_false(context.pace_thread.isAlive(),
                            'PaceThread not stopped')
    nose.tools.assert_true(context.finished_running.isSet(),
                           'Event not set properly')

def testNoWithContext():
    """
    Test the .Start() and .Finished()
    """
    context = Notif.NotifContext(METHOD = 'test')

    context.Start()
    nose.tools.assert_true(context.pace_thread.isAlive(),
                       'Did not create PaceThread')

    context.Finished()
    nose.tools.assert_false(context.pace_thread.isAlive(),
                            'PaceThread not stopped')

def testUpdateDelayed():
    """
    Test the UpdateMessage(MESSAGE)
    """    
    context = Notif.NotifContext(METHOD = 'test')

    nose.tools.assert_false(context.UpdateMessage('Just A Test'),
                           'UpdateMessage returned wrong value')
    
    nose.tools.assert_equal(context.this_cont_ann.pop(),'Just A Test',
                            'Test message not appended')

def testUpdateImmediate():
    """
    Test the UpdateMessage(MESSAGE, NOW = True)
    """
    context = Notif.NotifContext(METHOD = 'test')
    nose.tools.assert_true(context.UpdateMessage('Just A Test', NOW = True),
                           'UpdateMessage returned wrong value')


def testUpdateNow():
    """
    Test the UpdateNow(MESSAGE)
    """
    context = Notif.NotifContext(METHOD = 'test')
    nose.tools.assert_true(context.UpdateNow('Just A Test'),
                           'UpdateNow returned wrong value')
    






