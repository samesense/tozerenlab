import nose.tools
import PyGenome


def testGenomeLoad():
    genome = PyGenome.PyGenome('AAAAAAAACTAAA', 'AAAAAAAACGAAA')
    nose.tools.assert_not_equal(genome, None)

def testGenomeAlign():
    genome = PyGenome.PyGenome('AAAAAAAACTAAA', 'AAAAAAAACAAA')
    genome.AlignSeqs()
    nose.tools.assert_not_equal(genome.alignment, None)
    nose.tools.assert_equal(genome.alignment[1][-1],'-')
