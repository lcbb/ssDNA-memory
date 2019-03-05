#!/usr/local/bin/python3.4
#
# C. Tyson R. Shepherd, PhD
#
# Digital files stored in DNA format.
# Input from web.
#
# Output: Scaffold and staple sequences to order
#
import sys
from struct import calcsize
from struct import pack
import dnaDecrypt

if len(sys.argv) != 2:
    print("Usage: python3 ssDNAdecode.py <input_file>");
    sys.exit(1)
#print("Content-Type: text/html\n\n")
print("Output");
#
# Oligos used to amplify, sequence or extract all messages.
#
masterStart='CTTGGGTGGAGAGGCTATTC';
masterStop='GATCTCCTGTCATCTCACCT';
msgStop='GTACTAGTCGACGCGTGGCC';
#
inPath = sys.argv[1]
inFile = open(inPath)
inSeq = ''.join(line.strip() for line in inFile.readlines())
print("Input");
print(str(inSeq));
fileSeq=inSeq[:inSeq.find(msgStop)];
fType=dnaDecrypt.seqTypeBarcode(fileSeq[0:20]);
salt=fileSeq[20:36]
pageNum=dnaDecrypt.DNA2int(fileSeq[36:40])
fSize=dnaDecrypt.DNA2int(fileSeq[40:48])
fileSeq=fileSeq[48:]
fileBits=dnaDecrypt.dnaDecoder(fileSeq)
fileBits=[fileBits[8*i:8*(i+1)] for i in range(int(len(fileBits)/8))]
fileBits=[int(i, 2) for i in fileBits]
fileBytes=b''
for i in fileBits:
    fileBytes=fileBytes+pack('B', i);
fOut=dnaDecrypt.fileDecryptor(fileBytes, salt)
if fType=='txt':
	print(fOut.decode("utf-8").rstrip());
if fType=='bmp' or fType=='jpg' or fType=='png':
	fOutSt=salt+'.'+fType
	fBMPout=open(fOutSt, 'wb')
	fBMPout.write(fOut)
	fBMPout.close()
