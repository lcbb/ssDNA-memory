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
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import bitstring
import random
import math
import os, sys
import re
import dnaMemFuncs
#
# Oligos used to amplify, sequence or extract all messages.
#
masterStart='CTTGGGTGGAGAGGCTATTC';
masterStop='GATCTCCTGTCATCTCACCT';
msgStop='GTACTAGTCGACGCGTGGCC';
#
# Process form to python
#
inPath = sys.argv[1]
inFile = open(inPath)
#inText = ''.join(line.rstrip() for line in inFile.readlines())
#print("Input");
#print(str(inText));

fileType = "text"
msgStart = str(dnaMemFuncs.fileTypeBarcode(fileType));
blockLen = int(sys.argv[2])
ScafSeqs=[]
gBlocks=[]
eflag=1
ivSeq=dnaMemFuncs.randDNA(16);
cCount=0
#
# File conversion to scaffold seqs
#
if inFile:
	while 1:
		oriFileCont=inFile.read(65536)
		if not oriFileCont: break
		cCount=cCount+1
		if cCount>255:
			print("File is too large, it will take too long to process")
			break
		fPage=dnaMemFuncs.NTintFill(cCount, 4)
		csSeq=dnaMemFuncs.NTintFill(sys.getsizeof(oriFileCont), 8)
		txtEnc=dnaMemFuncs.fileEncryptor(oriFileCont, ivSeq);
		l=0
		flag=1
		dnaNewDat=0
		while flag:
			if dnaNewDat==0:
				dnaNewDat=dnaMemFuncs.dnaEncoder(txtEnc);
				l=l+1
				dnaDat=masterStart+masterStop+msgStart+ivSeq+fPage+csSeq+dnaNewDat+msgStop
			else:
				dnaDat=dnaNewDat
			dnaNewDat=dnaMemFuncs.blastTester(dnaDat);
			if dnaNewDat == dnaDat:
				flag=0
			if l>50:
				flag=0
				print("Too many repeats - Somthing went wrong, try again.")
				break
		dnaDat=dnaNewDat[40:]
		print('Full bitstream sequence: '+dnaDat)
		lDat=len(str(dnaDat));
#
		if blockLen>=lDat:
			ScafSeqs.append(str(dnaDat)+dnaMemFuncs.randDNA(blockLen-(lDat)))
		else:
			numObjs=int(math.ceil((float(lDat)-float(blockLen))/float(blockLen-50)))+1
			print("\nNumber of Objects: "+str(numObjs)+"\n\n")
			mplr=objLen-(numObjs*blockLen-lDat)/(numObjs-1)
			for i in range(0,numObjs):
				cntr=int(blockLen/2+i*mplr)
				ScafSeqs.append(dnaDat[int(cntr-(blockLen/2)):int(cntr+(blockLen/2))])
#
		for i in ScafSeqs:
			gBlocks.append(masterStart+i+masterStop);
#
		for i in range(len(gBlocks)):
			print('Block sequences '+str(i+1)+":\n"+gBlocks[i]+"\n\n");
#
else:
	print('Error: no file');
