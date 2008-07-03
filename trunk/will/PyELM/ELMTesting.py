import PyELM
import pickle
import os
import nose.tools


class SeqAns:
    def __init__(self):
        self.seq = []
        self.SimpleAns = []
        self.MultiAns = []
        self.ELM = []

        

def testImport():
    test = PyELM.PyELM()
    nose.tools.assert_not_equal(test, None)
        
