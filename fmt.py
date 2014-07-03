__author__ = 'ralf'

from Bio import SeqIO
from re import sub
import os


def abi2fasta(abi_file, fasta_file=None):
        # convert ab1 file to fasta file
        if not fasta_file:
                fasta_file = sub(r'\.ab1', '.fa', abi_file)
                if fasta_file == abi_file:
                        fasta_file = abi_file + '.fa'

        abi_handle = open(abi_file, 'rb')
        fasta_handle = open(fasta_file, 'w')

        seq = SeqIO.parse(abi_handle, 'abi')
        SeqIO.write(seq, fasta_handle, 'fasta')

        abi_handle.close()
        fasta_handle.close()


def abi2fastq(abi_file, fastq_file=None):
        #convert ab1 file to fastq file
        if not fastq_file:
                fastq_file = sub(r'\.ab1', '.fq', abi_file)
                if fastq_file == abi_file:
                        fastq_file = abi_file + '.fq'

        abi_handle = open(abi_file, 'rb')
        fastq_handle = open(fastq_file, 'w')

        seq = SeqIO.parse(abi_handle, 'abi')
        SeqIO.write(seq, fastq_handle, 'fastq-sanger')

        abi_handle.close()
        fastq_handle.close()


def flat2fasta(flat_file, fasta_file=None):
        #convert flat file to fasta
        if not fasta_file:
                fasta_file = sub(r'.txt', '.fa', flat_file)
                if fasta_file == flat_file:
                        fasta_file = flat_file + '.fa'

        flat_handle = open(flat_file, 'r')
        fasta_handle = open(fasta_file, 'w')

        seq = SeqIO.parse(flat_handle, 'genbank')
        SeqIO.write(seq, fasta_handle, 'fasta')

        flat_handle.close()
        fasta_handle.close()


def abi_converter(abi_file, _fmt, out_file=None):
        #convert abi file to fasta or fastq
        if not os.path.isfile(out_file):

                if _fmt == 'fa':
                        abi2fasta(abi_file, out_file)
                elif _fmt == 'fq':
                        abi2fastq(abi_file, out_file)
                else:
                        print("Format should be 'fa' or 'fq'")
