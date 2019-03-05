#
# C. Tyson R. Shepherd, PhD
#
# Function container for digital file storage in DNA
#
#
#
import sys
import os
import string
from struct import unpack
from struct import pack
from struct import calcsize
from bitstring import BitArray
from lzma import decompress
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
