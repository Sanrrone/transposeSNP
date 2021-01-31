#!/usr/bin/python3
from __future__ import with_statement
# ==============================================================================
# 						multiGenomicContext
#
# Author: Sandro Valenzuela (sandrolvalenzuead@gmail.com) 
#
# Please type "python multiGenomicContext.py -h" for usage help
#
# ==============================================================================


import sys, os, re, subprocess, csv, glob
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from pysam import VariantFile


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


# basic input
parser = OptionParser(
    usage="Usage: python tsnp.py -v myvcffile.vcf -r myroginalref.fasta -t mytargetref.fasta -o tnsp.vcf")
parser.add_option("-v", "--vcf", dest="vcffile", help="Your VCF file, it can also be a BCF file")
parser.add_option("-r", "--ref", dest="reffile", help="Your reference used to obtain your VCF (fasta)")
parser.add_option("-t", "--target", dest="target", help="target fasta file to search the coords from the original VCF")
parser.add_option("-o", "--output", dest="output", help="Output file name", default='tsnp.vcf')
# advance input
parser.add_option("-d", "--dp", dest="dp", help="Min DP to consider valid to search", default=0)
parser.add_option("-w", "--windowsize", dest="wsize",
                  help="window size to search in target fasta. Consider the snp inside and in the middle of the windows size",
                  default=50)
parser.add_option("-b", "--blastidx", dest="blastidx",
                  help="The Blast index for your target fasta, it will be used instead of (-t) for 'blastn -db'",
                  default=None)
# parser.add_option("-l","--length",dest="minLen",help="minimum length for resulting sequences after remove the query match sequence", default=0)
# parser.add_option("-v","--inverse",dest="inverse",help="inverse the match to keep the query sequience or equivalent in subject sequences", default=False, action='store_true')
# parser.add_option("-t","--translate",dest="translate",help="translate the output sequence if they are DNA", default=False, action='store_true')


(options, args) = parser.parse_args()
# basic inputs
vcffile = options.vcffile
reffile = options.reffile
targetfile = options.target
output = options.output
# advance inputs
dp = options.dp
wsize = options.wsize
blastidx = options.blastidx
# alignL = int(options.alignL)
# minLen = int(options.minLen)
# inverse = options.inverse
# translate = options.translate

#######################################################################################################################
# check variables
#######################################################################################################################
if not vcffile:
    print("ERROR: No input provided (-v), use -h for help")
    sys.exit()
else:
    if not os.path.isfile(vcffile):
        print(str("ERROR: " + vcffile + " doesn't exist, check the file directory"))
        sys.exit()

if not reffile:
    print("ERROR: No reference provided (-r), use -h for help")
    sys.exit()
else:
    if not os.path.isfile(reffile):
        print(str("ERROR: " + reffile + " doesn't exist, check the file directory"))
        sys.exit()

if not targetfile:
    if not blastidx:
        print("ERROR: No target provided (-t or -b for blast index), use -h for help")
        sys.exit()
else:
    if not os.path.isfile(targetfile):
        print(str("ERROR: " + targetfile + " doesn't exist, check the file directory"))
        sys.exit()

if not output or output == "":
    output = "tsnp.vcf"

if wsize < 10:
    wsize = 10

if wsize % 2 == 1:
    wsize = wsize + 1

# searching for blastn
if which('blastn') is None:
    print("ERROR: No blastn found, install it before continue")
    sys.exit()

#######################################################################################################################

print("* reading vcf file")
vcf_in = VariantFile(vcffile)  # auto-detect input format
print("* reading original reference")
reffasta = SeqIO.to_dict(SeqIO.parse(reffile, "fasta"))
print("* reading target reference")

if not blastidx:
    command = str(
        "blastn -query tmp.fasta -subject " + targetfile + " -out tmp.out -evalue 0 -outfmt 10 -qcov_hsp_perc 100 -perc_identity 100 -max_target_seqs 1")
else:
    command = str(
        "blastn -query tmp.fasta -db " + blastidx + " -out tmp.out -evalue 0.001 -outfmt 10 -num_threads " + str(
            os.cpu_count()) + " -qcov_hsp_perc 100 -perc_identity 100 -max_target_seqs 1")

outfile = open(output, 'w')
outfile.write("%s" % vcf_in.header)

for rec in vcf_in.fetch():
    #outfile.write("%s" % rec)
    pos = rec.pos
    contig = rec.chrom
    print("* working on contig", contig + ":", pos)
    # print(reffasta.keys())
    if rec.info['DP'] >= dp:
        seq = reffasta[contig]
        clength = len(seq.seq)

        # check if the windows size is between a valid reference position
        if pos + wsize <= clength:
            seqpos = 0 #means the snp is on the begining of pos.
            refseq = seq.seq[pos-1:pos+wsize]
        elif pos - wsize >= 0:
            seqpos = int(pos - wsize)
            refseq = seq.seq[pos - 1 - wsize: pos] # means the snp is after wsize nucleotides
        else:
            print("* position", pos, "is in a smaller sequence than the winows size (" + wsize + ")")
            continue
        #####################################################################


        tmp = open('tmp.fasta', 'w')
        tmp.write(">%s\n%s\n" % (contig + "___" + str(pos), refseq))
        tmp.close()

        # calling blast
        subprocess.call(command, shell=True)

        ### parsing results
        if os.path.getsize("tmp.out") > 0:
            tmp = open("tmp.out", "r")
            for uniquerow in csv.reader(tmp, delimiter=','):
                # uniquerow[0] is our query sequence
                # uniquerow[1] is the name of the sequence that match with our query (header of the fasta to be specific)
                # uniquerow[2] is identity
                # uniquerow[3] is alignment coverage (length)
                # uniquerow[8] subject sequence start
                # uniquerow[9] subject sequence end

                #print(rec.filter, rec.info, rec.format, rec.samples)
                if int(uniquerow[8]) < int(uniquerow[9]):
                    subjectpos = uniquerow[8]
                else:
                    break


                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t.\t" % (uniquerow[1], subjectpos, '.' if rec.id is None else rec.id, rec.ref, ','.join(rec.alts), str(int(rec.qual))))
                info = []
                #print(dir(rec.filter))
                for field, value in rec.info.items():
                    if isinstance(value, tuple):
                        value = ','.join(str(i) for i in value)

                    if isinstance(value, float):
                        value = round(value,7)

                    info.append(field+"="+str(value))

                outfile.write("%s\t" % (';'.join(info)))

                recformat = ':'.join(rec.format.keys())
                outfile.write("%s" % recformat)

                for sample, recsample in rec.samples.items():
                    outfile.write("\t%s:%s" % ('/'.join(str(i) for i in rec.samples[sample]['GT']), ','.join(str(i) for i in rec.samples[sample]['PL'])))

                outfile.write("\n")

        else:
            print("snp", contig, pos, "not found in target file")


    else:
        print("* " + contig, pos, "with low DP, next")

outfile.close()
os.remove("tmp.out")
os.remove("tmp.fasta")
