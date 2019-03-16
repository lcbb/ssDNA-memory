# Digital files stored in DNA format
#    Copyright (C) 2019  Tyson R. Shepherd, PhD
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.#
# 
# Input: Name of file containing DNA sequence
#
# Output: Decrytped and decoded DNA to digital file
#
#
import sys
from struct import calcsize
from struct import pack
import dnaDecrypt

if len(sys.argv) != 2:
    print("Usage: python3 ssDNAdecode.py <input_file>");
    sys.exit(1)
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
print("Input:");
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
print("\n\nOutput:");
if fType=='txt':
	print(fOut.decode("utf-8").rstrip());
if fType=='bmp' or fType=='jpg' or fType=='png':
	fOutSt=salt+'.'+fType
	fBMPout=open(fOutSt, 'wb')
	fBMPout.write(fOut)
	fBMPout.close()
