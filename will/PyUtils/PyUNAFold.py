import subprocess
import os
import re

def hybrid_min(Seq1,Seq2,Type='RNA',Temp=37,Sodium=1,Magnesium=0,UNAFoldPath='C:\\UNAFold\\bin\\'):
    """
    hybrid_min(Seq1,Seq2)
        Interfaces with the hybrid_min function of UNAFold and returns the free-energy binding of Seq1 with Seq2.
        Returns a tuple containing the three output results.

        Type=['RNA']|'DNA'
        Temp=[37]
        Sodium=[1]
        Magnesium=[0]
        UNAFoldPath=['c:\\UNAFold\\bin\\']

    """

    command = UNAFoldPath + 'hybrid-min '
    command += ' --tmin=' + str(Temp)
    command += ' --tmax=' + str(Temp)
    command += ' --sodium=' + str(Sodium)
    command += ' --magnesium=' + str(Magnesium)
    command += ' --NA=' + Type
    command += ' -q '

    command += Seq1 + ' ' + Seq2

#    print command
    SysCall = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
    SysCall.wait()
    output = SysCall.communicate()[0]

    outputList = re.split('\t',output)

    return (float(outputList[0]),float(outputList[1]),float(outputList[2]))
    

def hybrid(Seq1,Seq2,Tmin,Tmax,Tstep=1,Type='RNA',Sodium=1,Magnesium=0,UNAFoldPath='C:\\UNAFold\\bin\\'):
    """
    hybrid(Seq1,Seq2,Tmin,Tmax,Tstep=1)
        Interfaces with the hybrid_min function of UNAFold and returns the free-energy binding of Seq1 with Seq2 across multiple temperature ranges.
        Returns a list of tuples containing:
            (Temp,Out1,Out2,Out3)

        Type=['RNA']|'DNA'
        Sodium=[1]
        Magnesium=[0]
        UNAFoldPath=['c:\\UNAFold\\bin\\']

    """

    command = UNAFoldPath + 'hybrid-min '
    command += ' --tmin=' + str(Tmin)
    command += ' --tmax=' + str(Tmax)
    command += ' --tinc=' + str(Tstep)
    command += ' --sodium=' + str(Sodium)
    command += ' --magnesium=' + str(Magnesium)
    command += ' --NA=' + Type
    command += ' -q '

    command += Seq1 + ' ' + Seq2

#    print command
    SysCall = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
    SysCall.wait()
    output = SysCall.communicate()[0]

    outputList = re.split('\n',output)

    finalList = []
    for rowNum in range(len(outputList)-1):
        temp = re.split('\t',outputList[rowNum])
        print temp
        finalList.append((Tmin+rowNum*Tstep,float(temp[0]),float(temp[1]),float(temp[2])))

    return finalList
    
