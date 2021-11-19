# process-sample.py

# extract reads and make mitochondrial calls


import callmito
import os
import sys
import argparse
import time

# SETUP

parser = argparse.ArgumentParser(description='process-sample.py')

parser.add_argument('--ref', type=str,help='genome fasta with dictionary and .fai',required=True)
parser.add_argument('--finaldir', type=str,help='final dir for output',required=True)
parser.add_argument('--name', type=str,help='name of sample to process',required=True)
parser.add_argument('--cram', type=str,help='aligned cram file',required=True)
parser.add_argument('--coords',type=str,help='coordinates to extract, numts + chrM',required=True)

args = parser.parse_args()

#####################################################################

myData = {} # dictionary for keeping and passing information

myData['finalDir'] = args.finaldir
myData['ref'] = args.ref
myData['sampleName'] = args.name

myData['cramFileName'] = args.cram
myData['coordsFileName'] = args.coords


if myData['finalDir'][-1] != '/':
   myData['finalDir'] += '/'


if os.path.isdir(myData['finalDir']) is False:
    print('Error! output dir %s not does not exist' % myData['finalDir'])
    sys.exit()
    
# setup the output dir
myData['finalDirSample'] = myData['finalDir']  + myData['sampleName']

if os.path.isdir(myData['finalDirSample']) is False:
    print('making ',myData['finalDirSample'])
    cmd = 'mkdir ' + myData['finalDirSample']
    print(cmd)
    callmito.runCMD(cmd)
myData['finalDirSample'] += '/'    

myData['logFileName'] = myData['finalDirSample'] + myData['sampleName'] + '.mito.log'
myData['logFile'] = open(myData['logFileName'],'w')

# add initial infoto log
callmito.init_log(myData)
callmito.check_prog_paths(myData)

# get reads to extract
callmito.extract_reads(myData)



myData['logFile'].close()