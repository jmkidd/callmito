# callmito.py
# functions for going from aligned CRAM to mitochondrial calls

import sys
import subprocess
import os
import argparse
import time
import socket
import shutil
import gzip


###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed')
        print(cmd)
        sys.exit(1)
###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD_output(cmd):
    val = subprocess.Popen(cmd, text=True, shell=True, stdout = subprocess.PIPE)
    resLines = []
    for i in val.stdout:
       i = i.rstrip()
       resLines.append(i)
    return resLines
#############################################################################        
# setup paths to default programs to use and checks for required programs
def check_prog_paths(myData):        
    myData['logFile'].write('\nChecking for required programs...\n')
    
    for p in ['bwa','gatk','samtools']:
        if shutil.which(p) is None:
            s = p + ' not found in path! please fix (module load?)'
            print(s, flush=True)
            myData['logFile'].write(s + '\n')
            myData['logFile'].close()        
            sys.exit()
        else:
            myData['logFile'].write('%s\t%s\n' % (p,shutil.which(p)))
            
    myData['logFile'].flush()              
############################################################################# 
def init_log(myData):
    k = list(myData.keys())
    k.sort()
    myData['startTime'] = time.localtime()
    myData['tStart'] = time.time()
    t = time.strftime("%a, %d %b %Y %H:%M:%S", myData['startTime'])        
    myData['logFile'].write(t + '\n')
    
    hn = socket.gethostname()
    myData['logFile'].write('Host name: %s\n' % hn)
    print('Host name: %s\n' % hn,flush=True)
    
    myData['logFile'].write('\nInput options:\n')
    for i in k:
        if i in ['logFile']:
            continue        
        myData['logFile'].write('%s\t%s\n' % (i,myData[i]))                
    myData['logFile'].flush()  
############################################################################# 
def parse_sam_line(myLine):
    res = {}
    res['seqName'] = myLine[0]
    res['flag'] = int(myLine[1])
    res['chrom'] = myLine[2]
    res['chromPos'] = int(myLine[3])
    res['mapQ'] = int(myLine[4])
    res['cigar'] = myLine[5]
    res['seq'] = myLine[9]
    
    #
    res['seqLen'] = len(myLine[9])
    
    
    res['cigarExpand'] = expand_cigar(res['cigar'])
    res['qual'] = myLine[10]
    res['mateChrom'] = myLine[6]
    res['matePos'] = myLine[7]
    
    res['fragLen'] = int(myLine[8])
    
    res['cigarCounts']={}
    res['cigarCounts']['M'] = 0
    res['cigarCounts']['D'] = 0
    res['cigarCounts']['I'] = 0
    res['cigarCounts']['S'] = 0
    res['cigarCounts']['H'] = 0

    
    if res['flag'] & 0x10 != 0:
        res['reverseStrand'] = True
    else:
        res['reverseStrand'] = False

    if res['flag'] & 0x4 != 0:
        res['unMapped'] = True
    else:
        res['unMapped'] = False

    if res['flag'] & 0x400 != 0:
        res['isDuplicate'] = True
    else:
        res['isDuplicate'] = False

    if res['flag'] & 0x100 != 0:
        res['notPrimaryAlignment'] = True
    else:
        res['notPrimaryAlignment'] = False

    if res['flag'] & 0x800 != 0:
        res['isSupplementaryAlignment'] = True
    else:
        res['isSupplementaryAlignment'] = False


    if res['flag'] & 0x1 != 0:
        res['isPaired'] = True
    else:
        res['isPaired'] = False


    if res['flag'] & 0x8 != 0:
        res['mateUnmapped'] = True
    else:
        res['mateUnmapped'] = False

    if res['flag'] & 0x40 != 0:
        res['isFirst'] = True
    else:
        res['isFirst'] = False

        
    for i in res['cigarExpand']:
        res['cigarCounts'][i[1]] += i[0]
        
    # check for proper seqlen to update, 2015-05-05
    if myLine[9] == '*':  #not actually sequence present in SAM line
        res['seqLen'] = res['cigarCounts']['M'] + res['cigarCounts']['I']  + res['cigarCounts']['S'] + res['cigarCounts']['H']
     
     
    res['otherTags'] = myLine[11:]
             
    return res
#####################################################################
#returns lists of [int,flag]
def expand_cigar(cigar):
    res = []
    if cigar == '*':
        return res
    digits = ['0','1','2','3','4','5','6','7','8','9']
    accumulate = ''
    i = 0
    while True:
        if i == len(cigar):
            break
        if cigar[i] in digits:
            accumulate += cigar[i]
            i += 1
        else:
            d = int(accumulate)
            res.append([d,cigar[i]])
            i += 1
            accumulate = ''
    return res
#####################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c   
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''    
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
##############################################################################
def get_seq_from_sam(samRec):
    name = samRec['seqName']
    if samRec['isFirst'] is True:
        readNum = 1
    else:
        readNum = 2
    
    # get if the seq...
    if samRec['reverseStrand'] is True:
        seq = samRec['seq']
        seq = revcomp(seq)
        qual =  samRec['qual']
        qual = qual[::-1] #reverse
    else:
        seq = samRec['seq']
        qual =  samRec['qual']
    return([name,readNum,seq,qual])         
############################################################################# 
# define critera to extract reads for remapping
def to_extract(samRec):
    if samRec['unMapped'] is True:
        return False

    if samRec['isDuplicate'] is True:
        return False

    if samRec['notPrimaryAlignment'] is True:
        return False

    if samRec['isPaired'] is False:
        return False

    return True
###############################################################################        


###############################################################################        
def extract_reads(myData):
    searchDelta = 100 # space to search across
    
    myData['logFile'].write('\nstarting extraction of fastq\n')
    
    cmd = 'samtools view -T %s -M -L %s %s ' % (myData['ref'], myData['coordsFileName'], myData['cramFileName'])
    print(cmd)
    myData['logFile'].write(cmd + '\n')
    
    myData['readData'] = {}  # dictionary to store all of read 1 and read 2 info
    myData['readsToExtract'] = {} # dictionary of reads to extract
    
    samLines = runCMD_output(cmd)
    for samLine in samLines:
        samLine = samLine.rstrip()
        samLine = samLine.split()
        samRec = parse_sam_line(samLine)
        
        seqInfo = get_seq_from_sam(samRec)
        if seqInfo[0] not in myData['readData']:
            myData['readData'][samRec['seqName']] = ['Empty','Empty']
        myData['readData'][samRec['seqName']][seqInfo[1]-1] = [seqInfo[2],seqInfo[3],samRec]  # seq, qual, then samRec to keep
        
        if to_extract(samRec) is True:
            myData['readsToExtract'][samRec['seqName']] = 1   
    
    
    s = 'Have total of %i reads pass extraction criteria' % len(myData['readsToExtract'])    
    print(s,flush=True)
    myData['logFile'].write(s + '\n')

    # get how many need extraction
    nMissing = 0
    rnsToGetMate = []    
    for rn in myData['readsToExtract']:
        if myData['readData'][rn][0] == 'Empty' or myData['readData'][rn][1] == 'Empty':
            nMissing += 1
            rnsToGetMate.append(rn)
    s = 'After initial read through, there are %i with missing mates' % nMissing
    print(s,flush=True)
    myData['logFile'].write(s + '\n')
    
    s = 'Starting pass 1 of cleanup of other read ends'            
    print(s,flush=True)
    myData['logFile'].write(s + '\n')
    
    for rn in rnsToGetMate:
        if myData['readData'][rn][0] != 'Empty' and myData['readData'][rn][1] != 'Empty': # already found
            continue
        
        if myData['readData'][rn][0] == 'Empty':
            rec = myData['readData'][rn][1][2]
        else:
            rec = myData['readData'][rn][0][2]
        
        # check out where the mate is
        mateChrom = rec['mateChrom']
        if mateChrom == '=':
            mateChrom = rec['chrom']
        matePos = int(rec['matePos'])
        
        searchStart = matePos - searchDelta
        searchEnd = matePos + searchDelta
        searchInt = '%s:%i-%i' % (mateChrom,searchStart,searchEnd)
        
        # search region
        cmd = 'samtools view -T %s %s %s ' % (myData['ref'], myData['cramFileName'], searchInt)
        
        samLines = runCMD_output(cmd)
        for samLine in samLines:
            samLine = samLine.rstrip()
            samLine = samLine.split()
            samRec = parse_sam_line(samLine)
        
            seqInfo = get_seq_from_sam(samRec)
            if seqInfo[0] not in myData['readData']:
                myData['readData'][samRec['seqName']] = ['Empty','Empty']
            myData['readData'][samRec['seqName']][seqInfo[1]-1] = [seqInfo[2],seqInfo[3],samRec]  # seq, qual, then samRec to keep
        
    # get how many need extraction after second pass
    nMissing = 0
    rnsToGetMate = []    
    for rn in myData['readsToExtract']:
        if myData['readData'][rn][0] == 'Empty' or myData['readData'][rn][1] == 'Empty':
            nMissing += 1
            rnsToGetMate.append(rn)
    s = 'After second pass, there are %i with missing mates' % nMissing
    print(s,flush=True)
    myData['logFile'].write(s + '\n')

    if nMissing != 0:
        s = 'ERROR! there are still misisng reads after the second pass!'
        print(s,flush=True)
        myData['logFile'].write(s + '\n')
        sys.exit()

    # ready to write the fastq to out
     
    myData['fastq1OutName'] = myData['finalDirSample'] + 'r1.fq.gz'
    myData['fastq2OutName'] = myData['finalDirSample'] + 'r2.fq.gz'
     
    out1 = gzip.open(myData['fastq1OutName'],'wt')
    out2 = gzip.open(myData['fastq2OutName'],'wt')
    for rn in myData['readsToExtract']:
        s = myData['readData'][rn][0][0]
        q = myData['readData'][rn][0][1]
        out1.write('@%s\n%s\n+\n%s\n' % (rn,s,q))

        s = myData['readData'][rn][1][0]
        q = myData['readData'][rn][1][1]
        out2.write('@%s\n%s\n+\n%s\n' % (rn,s,q))
    out1.close()
    out2.close()
     
    s = 'reads written to output fastq files!'
    print(s,flush=True)
    myData['logFile'].write(s + '\n')          
    myData['logFile'].flush() 
###############################################################################        




