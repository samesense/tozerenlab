"""
PyGoogle
    A module for controlling the google-docs spreadsheet.
"""
import PyMozilla
from urllib2 import URLError
import re

class GoogleSpreadInteract():
    """
    GoogleSpreadInteract
        A class for interfacing with the google-docs spreadsheet elements
    """
    def __init__(self, FORMLINK = None, TRUSTME = 0):
        """
        GoogleSpreadInteract()
            Instanciates the interaction module with no URL.  Use DefineURL to
            add a URL.
            

        GoogleSpreadInteract(FORMLINK = URL, TRUSTME = 0)
            If provided then the class is instanciated with a URL and the
                submit URL and feilds are determined automatically.
            If TRUSTME is set to 1 then no error checks will be performed
                before uploading.

        """
        self.form_link = FORMLINK
        self.form_feilds = None
        self.form_dict = {}
        self.moz_emu = PyMozilla.MozillaEmulator(cacher = None, trycount = 0)
        self.submit_link = None

        self.trustme = TRUSTME
        

        if self.form_link != None:
            try:
                formText = self.moz_emu.download(self.form_link)
            except URLError:
                print 'Could not access form URL: ', self.form_link
                raise URLError
            
            self.ParseForm(formText)

    def ParseForm(self, FORMTEXT):
        """
        ParseForm(formText)
            Parses out the text from the FORM URL to determine the available
            feilds and determines the submit URL.
        """

        form_extractor = re.compile('form action="(.*?)"')
        
        self.submit_link = form_extractor.findall(FORMTEXT)[0]

        if self.submit_link == None:
            print 'Could not find the Submit URL from the form'
            raise IndexError

        title_extractor = re.compile('"ss-q-title">(.*)</label>')
        self.form_feilds = title_extractor.findall(FORMTEXT)

        if self.form_feilds == None:
            print 'Could not find Any feilds'
            raise IndexError

        for i in range(0, len(self.form_feilds)):
            self.form_dict[self.form_feilds[i]] = 'single%3A' + str(i)

    def DefineURL(self, FORMLINK):
        """
        DefineURL(formURL)
            Redefines the FORM URL and calls ParseForm
        """

        self.form_link = FORMLINK
        self.form_feilds = None
        self.form_dict = {}
        formText = self.moz_emu.download(self.form_link)
        self.ParseForm(formText)
    

    def SubmitTupleList(self, TUPLELIST):
        """
        SubmitTupleList(TUPLELIST)
            Takes a list of len(2) tuples where the 0'th element is a
            feildname and the 1'th element is the data to be submitted.

        """
        post_data = ''

        for thistup in TUPLELIST:
            if thistup[0] in self.form_dict:
                post_data += self.form_dict[thistup[0]] + '='
                post_data += str(thistup[1]) + '&'
            else:
                print 'An unknown feild was provided: ', thistup[0]


        self.moz_emu.download(self.submit_link, postdata = post_data)

    def SubmitMultiRecords(self, RECORDLIST):
        """
        SubmitMultiRecords(RECORDLIST)
            Takes a list of lists of tuples.  Each element in the list is
            indiviually submitted to SubmitTupleList.
        """

        if self.trustme != 1:
            val = self.__PerformMultiRecordCheck(RECORDLIST)
            if val != 1:
                print 'A bad record was found: ', str(val)
                raise IndexError
            

        for thisrecord in RECORDLIST:
            self.SubmitTupleList(thisrecord)

        
        
    def __PerformTupleCheck(self, TUPLELIST):
        """
        __PerformTupleCheck(TUPLELIST)
            Checks the tupleList against the feilds and returns 1 if all
            provided feilds are recognized and 0 otherwise
        """

        for thistup in TUPLELIST:
            if not(thistup[0] in self.form_dict):
                return 0
        return 1

    def __PerformMultiRecordCheck(self, RECORDLIST):
        """
        __PerformMultiRecordCheck(recordList)
            Checks the recordList against the feilds and returns 0 if all
            provided feilds are recognized.  Otherwise it returns the
            first bad record.
        """
        for thisrecord in RECORDLIST:
            if self.__PerformTupleCheck(thisrecord) == 0:
                return thisrecord
        return 1

    def PrintFormFields(self):
        """
        PrintFormFields()
            Lists the fields found in the form.
        """        
        for field in self.form_dict.keys():
            print field

    def GetFormFields(self):
        """
        GetFormFields()
            Accessor for the form field dictionary.
        """        
        return self.form_dict
            





















