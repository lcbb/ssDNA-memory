# Function container for digital file storage in DNA
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
#
import sys
import os
import string
from struct import unpack
from struct import pack
from struct import calcsize
from bitstring import BitArray
from Crypto.Cipher import AES
from io import StringIO
#
def DNA2int(seq):
	lA=len(seq)
	numOut=0
	convert = {'T': 0, 'G': 1, 'A': 2, 'C': 3}
	for i in range(lA):
		numOut=numOut+(4**(lA-i-1)*int(convert[seq[i]]))
	return(numOut)
#
def dnaDecoder(fileSeq):
	fileBit=''
	for i in fileSeq:
		if (i=='C' or i=='A' or i=='M'):
			fileBit=fileBit+'0';
		if (i=='T' or i=='G' or i == 'K'):
			fileBit=fileBit+'1';
#	print(str(calcsize('H'))+'<br>')
	return(fileBit);
def seqTypeBarcode(fTypeSeq):
	ttb={'GTTTAAGGTCACATCGCATG': 'txt', 'TCTTGACAAACGTGTGCTTG': 'mp3', 'GTAGTTCAGACGCCGTTAAG': 'jpg', 'GGTTCCGTTTTACATTCCAG': 'tif', 'GTATGGCACGCCTAATCTGG': 'wav', 'GGTCTTTCATGCGTATAGTC': 'png', 'CTCTTAAAACTGGTATCACC': 'bmp', 'CTATTACGAGCGCTTGGATC': '7z', 'CTTGATGCTTTACAAGATCG': 'aes'}
	if str(fTypeSeq) in ttb:
		return (ttb[str(fTypeSeq)]);
	else:
		return ('');
def fileDecryptor(fileSeq, ivSeq):
	iv=b'';
	for i in ivSeq:
		iv=iv+pack('B', ord(i));
	key='ry%Tr*>2Y><NFv5aqAEhU@Q046Cy$n92'
	mode=AES.MODE_CBC
	decryptor = AES.new(key, mode, iv);
	return(decryptor.decrypt(fileSeq));
