import os
import re
import subprocess
import string

BASE_DIR = 'C:\\Documents and Settings\\Will\\My Documents\\'
CURRENT_DIRECTORIES = [BASE_DIR + 'PyUtils\\', BASE_DIR + 'PyELM\\']
RC_FILE = BASE_DIR + 'PyUtils\\WillRC.rc'
DEST_DIREC = BASE_DIR + 'PyLintResults\\'

SKIPPING_FILES = ['PyMozilla.py', 'CurrentWorking.py', 'testing_script.py']


files = []
for thisdir in CURRENT_DIRECTORIES:
    for thisfile in os.listdir(thisdir):
        if re.match('.*?\.py$',thisfile) != None:
            if thisfile not in SKIPPING_FILES:
                files.append(thisdir+thisfile)


os.chdir(DEST_DIREC)

command = 'pylint.bat '
command += '--rcfile="' + RC_FILE + '"'
command += ' "' + string.join(files, '" "') + '"'

print command

test = subprocess.Popen(command, shell = True)
test.wait()
