import nose.tools
import PyVirus
from testELM import SeqAns, FileSetup


def testViralSeq():
    genome = PyVirus.ViralSeq('AAAAAAAACTAAA')
    nose.tools.assert_not_equal(genome, None)
    
def testBkgSeq():
    genome = PyVirus.BkgSeq('AAAAAAAACTAAA')
    nose.tools.assert_not_equal(genome, None)

def testPatSeq():
    genome = PyVirus.PatSeq('AAAAAAAACTAAA')
    nose.tools.assert_not_equal(genome, None)

def testRefSeq():
    genome = PyVirus.RefSeq('AAAAAAAACTAAA')
    nose.tools.assert_not_equal(genome, None)

