#
# C. Tyson R. Shepherd, PhD 2019
#
# Function container for digital file storage in DNA
#
#
#
from Crypto.Random.random import randint
import sys
import os
import string
from struct import pack
from bitstring import BitArray
from Crypto.Cipher import AES
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from io import StringIO
#import matlab.engine
#
def randDNA(numB):
	randSeq=''
	for x in range(0,numB):
		i=randint(0,1)
		k=len(randSeq);
		if int(i)==0:
			if (int(randint(0,99))<45) and (randSeq[k-5:k-1]!='CCCC'):
				randSeq+='C'
			else:
				if (randSeq[k-5:k-1]=='AAAA'):
					randSeq+='C'
				else:
					randSeq+='A'
		elif int(i)==1:
			if (int(randint(0,99))<45) and (randSeq[k-5:k-1]!='GGGG'):
				randSeq+='G'
			else:
				if (randSeq[k-5:k-1]=='TTTT'):
					randSeq+='G'
				else:
					randSeq+='T'
	return (randSeq);
#
def newSubSeq(seq):
	outSeq=''
	for x in seq:
		i=randint(0,2)
		if i==0:
			outSeq=outSeq+x
		else:
			if x == 'C':
				outSeq=outSeq+'A'
			elif x=='A':
				outSeq=outSeq+'C'
			elif x=='G':
				outSeq=outSeq+'T'
			elif x=='T':
				outSeq=outSeq+'G'
	if len(outSeq) != len(seq):
		print("Something went wrong with the newSubSeq");
	return (outSeq)

def blastTester(seq):
#	blastn -query seq -task blastn-short -subject seq
# Create two sequence files
	seq1 = SeqRecord(Seq(str(seq)), id="seq1")
	tfOut=open("seq.fasta", "w")
	SeqIO.write(seq1, tfOut, "fasta")
	tfOut.close()
	newSeq=seq
#	print("blasting <br>")

# Run BLAST and parse the output as XML
	output = NcbiblastnCommandline(query="seq.fasta", task="blastn-short", word_size=7, evalue=10.0, gapextend=2, gapopen=5, penalty=-3, reward=2, subject="seq.fasta", outfmt=5)()[0]
	blast_result_record = NCBIXML.read(StringIO(output))
	out=1

# Print some information on the result
	for alignment in blast_result_record.alignments:
#		print(str(len(alignment.hsps)))
		if len(alignment.hsps) > 10:
			os.remove("seq.fasta")
			out=0
			return(out)
		elif len(alignment.hsps) == 1:
			os.remove("seq.fasta")
			return(seq)
		for hsp in alignment.hsps:
			if hsp.expect>0.1:
				qS=hsp.query_start
				qE=hsp.query_end
				sS=hsp.sbjct_start
				sE=hsp.sbjct_end
				if qS>68 and qE<(len(seq)-20):
					newSeq=newSeq[:qS]+newSubSeq(seq[qS:qE])+newSeq[qE:];
				elif sS>68 and sE<(len(seq)-20):
					if sS<sE:
						newSeq=newSeq[:sS]+newSubSeq(seq[sS:sE])+newSeq[sE:];
					elif sS>sE:
						newSeq=newSeq[:sE]+newSubSeq(seq[sE:sS])+newSeq[sS:];
				else:
					print('Problem with blastSeq, rerun the program')
					os.remove("seq.fasta");
					out=0
					return(out)
	os.remove("seq.fasta")
	return(newSeq);
def Int2DNA(numB):
	convert = {'0': 'T', '1': 'G', '2': 'A', '3': 'C'}
	newNum=''
	current=numB
	A=''
	while int(current)!=0:
		newNum=str(int(current%4))+newNum
		current=int(current/4)
	for i in newNum:
		A=A+convert[i]
	return(A);
#
def NTintFill(numB, fillNum):
	A=Int2DNA(numB)
	B=''
	for i in range(fillNum-len(A)):
		B=B+'T'
	for i in A:
		B=B+i
	return(B);
#
def Str2DNA(strA):
	strB=strA[1:5]
	strB=strB.upper()
	strA=strA.upper()
	A=''
	allowedStr=string.ascii_uppercase+string.digits
	for i in strB:
		if i not in allowedStr:
			print('String tags can only contain: '+allowedStr);
			sys.exit();
	convert = {'A': 'AGG', 'B': 'AGT', 'C': 'ATA', 'D': 'ATC', 'E': 'ATG', 'F': 'ATT', 'G': 'CAA', 'H': 'CAC', 'I': 'CAG', 'J': 'CAT',
	'K': 'CCA', 'L': 'CCC', 'M': 'CCG', 'N': 'CCT', 'O': 'CGA', 'P': 'CGC', 'Q': 'CGG', 'R': 'CGT', 'S': 'CTA', 'T': 'CTC', 'U': 'CTG', 'V': 'CTT',
	'W': 'GAA', 'X': 'GAC', 'Y': 'GAG', 'Z': 'GAT', '0': 'AAA', '1': 'AAC', '2': 'AAG', '3': 'AAT', '4': 'ACA', '5': 'ACC', '6': 'ACG', '7': 'ACT',
	'8': 'AGA', '9': 'AGC', '@': 'TAC', '$': 'TAG', '>': 'TAA', '#': 'TAT'}
	for i in strA:
		A=A+convert[i]
	return(A);
#
def dnaEncoder(txtEnc):
	dnaDat=''
	k=0
	for i in txtEnc:
		k=k+1
		if int(i)==0:
			if (int(randint(0,99))<45) and (dnaDat[k-6:k-1]!='CCCCC'):
				dnaDat+='C'
			else:
				if (dnaDat[k-6:k-1]=='AAAAA'):
					dnaDat+='C'
				else:
					dnaDat+='A'
		elif int(i)==1:
			if (int(randint(0,99))<45) and (dnaDat[k-5:k-1]!='GGGG'):
				dnaDat+='G'
			else:
				if (dnaDat[k-5:k-1]=='TTTT'):
					dnaDat+='G'
				else:
					dnaDat+='T'
	return(dnaDat);
def fileTypeBarcode(fType):
	ttb={'text': 'GTTTAAGGTCACATCGCATG', 'mp3': 'TCTTGACAAACGTGTGCTTG', 'jpg': 'GTAGTTCAGACGCCGTTAAG', 'tif': 'GGTTCCGTTTTACATTCCAG', 'wav': 'GTATGGCACGCCTAATCTGG', 'png': 'GGTCTTTCATGCGTATAGTC', 'bmp': 'CTCTTAAAACTGGTATCACC', '7z': 'CTATTACGAGCGCTTGGATC'}
	if str(fType) in ttb:
		return (ttb[str(fType)]);
	else:
		return ('CTTGATGCTTTACAAGATCG');
def fileEncryptor(fileIn, ivSeq):
	iv=b''
#	key='ThisIsADNABlock'
	key='ry%Tr*>2Y><NFv5aqAEhU@Q046Cy$n92'
	mode=AES.MODE_CBC
#	fsize=len(fileIn)
#	fileOut=pack('<H', fsize);
	for i in ivSeq:
		iv=iv+pack('B', ord(i))
	encryptor=AES.new(key, mode, iv);
	fileIn += ' ' * (16-len(fileIn)%16)
	fileOut=encryptor.encrypt(fileIn);
	fOut=BitArray(fileOut).bin[:];
	return(fOut);
def filePackager(fileIn, ivSeq):
	oriFileEnc=compress(fileIn,1,1);
#	key='ThisIsADNABlock1'
	key='ry%Tr*>2Y><NFv5aqAEhU@Q046Cy$n92'
	mode=AES.MODE_CBC
	fsize=len(oriFileEnc)
	fileOut=pack('<H', fsize);
	iv=b''
	for i in ivSeq:
		iv=iv+pack('B', ord(i))
	encryptor=AES.new(key, mode, iv);
	oriFileEnc += b' ' * (16-len(oriFileEnc)%16)
	fileOut=fileOut+encryptor.encrypt(oriFileEnc);
	fOut = BitArray(fileOut).bin[:];
##	fOut = ivBin+' '+bitstring.BitArray(fileOut).bin[:];
	return (fOut);
