from PyMozilla import *
import re

class GoogleSpreadInteract():
    def __init__(self,formLink=None,TrustMe=0):
        """
        GoogleSpreadInteract()
            Instanciates the interaction module with no URL.  Use DefineURL to add a URL.
            

        GoogleSpreadInteract(formLink=URL,TrustMe=0)
            If provided then the class is instanciated with a URL and the submit URL and feilds are determined automatically.
            If TrustMe is set to 1 then no error checks will be performed before uploading.

        """
        self.formLink = formLink
        self.formFeilds = None
        self.formDict = {}
        self.MozEmu = MozillaEmulator(cacher=None,trycount=0)
        self.submitLink = None

        self.TrustMe = TrustMe
        

        if self.formLink != None:
            try:
                formText = self.MozEmu.download(self.formLink)
            except:
                print 'Could not access form URL: ', self.formLink
                raise IndexError
            
            self.ParseForm(formText)

    def ParseForm(self,formText):
        """
        ParseForm(formText)
            Parses out the text from the FORM URL to determine the available feilds and determines the submit URL.
        """
        
        self.submitLink = re.compile('form action="(.*?)"').findall(formText)[0]

        if self.submitLink == None:
            print 'Could not find the Submit URL from the form'
            raise IndexError

        self.formFeilds = re.compile('"ss-q-title">(.*)</label>').findall(formText)

        if self.formFeilds == None:
            print 'Could not find Any feilds'
            raise IndexError

        for i in range(0,len(self.formFeilds)):
            self.formDict[self.formFeilds[i]] = 'single%3A' + str(i)

    def DefineURL(self,formLink):
        """
        DefineURL(formURL)
            Redefines the FORM URL and calls ParseForm
        """

        self.formLink = formLink
        self.formFeilds = None
        self.formDict = {}
        formText = self.MozEmu.download(self.formLink)
        self.ParseForm(formText)
    

    def SubmitTupleList(self,tupleList):
        """
        SubmitTupleList(tupleList)
            Takes a list of len(2) tuples where the 0'th element is a feildname and the 1'th element is the data to be submitted.

        """
        postData = ''

        for thisTup in tupleList:
            if thisTup[0] in self.formDict:
                postData += self.formDict[thisTup[0]] + '=' + str(thisTup[1]) + '&'
            else:
                print 'An unknown feild was provided: ', thisTup[0]


        outputData = self.MozEmu.download(self.submitLink,postdata=postData)

    def SubmitMultiRecords(self,recordList):
        """
        SubmitMultiRecords(recordList)
            Takes a list of lists of tuples.  Each element in the list is indiviually submitted to SubmitTupleList.
        """

        if self.TrustMe != 1:
            val = self.__PerformMultiRecordCheck(recordList)
            if val != 1:
                print 'A bad record was found: ', str(val)
                raise IndexError
            

        for thisRecord in recordList:
            self.SubmitTupleList(thisRecord)

        
        
    def __PerformTupleCheck(self,tupleList):
        """
        __PerformTupleCheck(tupleList)
            Checks the tupleList against the feilds and returns 1 if all provided feilds are recognized and 0 otherwise
        """

        for thisTup in tupleList:
            if not(thisTup[0] in self.formDict):
                return 0
        return 1

    def __PerformMultiRecordCheck(self,recordList):
        """
        __PerformMultiRecordCheck(recordList)
            Checks the recordList against the feilds and returns 0 if all provided feilds are recognized.  Otherwise it returns the
            first bad record.
        """
        for thisRecord in recordList:
            if self.__PerformTupleCheck(thisRecord) == 0:
                return thisRecord
        return 1
            





















