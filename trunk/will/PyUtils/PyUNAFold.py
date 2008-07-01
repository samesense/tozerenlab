"""
PyUNAFold
    An interface for the UNAFold program for finding siRNA templates.
"""

import subprocess
import re
import numpy

def HybridMin(SEQ1, SEQ2, TYPE = 'RNA', TEMP = 37, SODIUM = 1, MAGNESIUM = 0,
               UNAFOLDPATH = 'C:\\UNAFold\\bin\\', RETURNTYPE = 'tuple'):
    """
    HybridMin(Seq1,Seq2)
        Interfaces with the hybrid_min function of UNAFold and returns the
        free-energy binding of Seq1 with Seq2.  Returns a tuple containing
        the three output results.

        TYPE=['RNA']|'DNA'
        TEMP=[37]
        SODIUM=[1]
        MAGNESIUM=[0]
        UNAFOLDPATH=['c:\\UNAFold\\bin\\']
        RETURNTYPE=['tuple']|'array'
            Returns the data either as a tuple of numpy array

    """

    command = UNAFOLDPATH + 'hybrid-min '
    command += ' --tmin=' + str(TEMP)
    command += ' --tmax=' + str(TEMP)
    command += ' --sodium=' + str(SODIUM)
    command += ' --magnesium=' + str(MAGNESIUM)
    command += ' --NA=' + TYPE
    command += ' -q '

    command += SEQ1 + ' ' + SEQ2

#    print command
    sys_call = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
    sys_call.wait()
    output = sys_call.communicate()[0]

    output_list = re.split('\t', output)

    if RETURNTYPE == 'tuple':
        return (float(output_list[0]), float(output_list[1]),
                float(output_list[2]))
    elif RETURNTYPE == 'array':
        return numpy.array([float(output_list[0]), float(output_list[1]),
                            float(output_list[2])])
    

def Hybrid(SEQ1, SEQ2, TMIN, TMAX, TSTEP = 1, TYPE = 'RNA', SODIUM = 1,
           MAGNESIUM = 0, UNAFOLDPATH = 'C:\\UNAFold\\bin\\',
           RETURNTYPE = 'tuple'):
    """
    Hybrid(SEQ1,SEQ2,TMIN,TMAX,TSTEP=1)
        Interfaces with the hybrid_min function of UNAFold and returns the
        free-energy binding of Seq1 with Seq2 across multiple temperature
        ranges.

        Returns a list of tuples containing:
            (Temp,Out1,Out2,Out3)

        TYPE=['RNA']|'DNA'
        SODIUM=[1]
        MAGNESIUM=[0]
        UNAFOLDPATH=['c:\\UNAFold\\bin\\']
        RETURNTYPE=['tuple']|'array'
            Returns the data either as a tuple of numpy array
    """

    command = UNAFOLDPATH + 'hybrid-min '
    command += ' --tmin=' + str(TMIN)
    command += ' --tmax=' + str(TMAX)
    command += ' --tinc=' + str(TSTEP)
    command += ' --sodium=' + str(SODIUM)
    command += ' --magnesium=' + str(MAGNESIUM)
    command += ' --NA=' + TYPE
    command += ' -q '

    command += SEQ1 + ' ' + SEQ2

#    print command
    sys_call = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
    sys_call.wait()
    output = sys_call.communicate()[0]

    output_list = re.split('\n', output)

    if RETURNTYPE == 'tuple':
        final_list = []
        for row_num in range(len(output_list)-1):
            temp = re.split('\t', output_list[row_num])
            print temp
            final_list.append((TMIN+row_num*TSTEP, float(temp[0]),
                               float(temp[1]), float(temp[2])))

        return final_list
    elif RETURNTYPE == 'array':
        final_array = numpy.zeros((len(output_list)-1, 4))
        for row_num in range(len(output_list)-1):
            temp = re.split('\t', output_list[row_num])
            final_array[row_num] = [TMIN+row_num*TSTEP, float(temp[0]),
                                   float(temp[1]), float(temp[2])]

        return final_array







        
