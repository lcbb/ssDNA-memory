# ssDNA-memory
Single-stranded, kilobase-scale DNA memory encoding scheme for archival in phage

You will need Python3 with the pycrypto and bitstring libraries installed to run the decryption algorithm.

Decode the file on the command line by typing:

python3 ssDNAMdecode.py filename
  
Filename is the name of the text file containing the sequence of DNA cloned to the phage. Three example text files are provided encoding lines from The Crucible, Waiting for Godot, and Hamlet.


To encode a new file on the command line type:

python3 ssDNAMcreate.py filename blockSize

filename is the name of the file to encode in DNA
blockSize is the size of the DNA to generate, either breaking into blocks of that size, or filling with slack space

You will need pycrypto and biopython python packages installed and NCBI blast command line tools to run the encoding algorithm.
