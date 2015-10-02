# callMom
Python script for calling SNPs, MNPs, and indels in mitochondrial DNA (mtDNA). MtDNA is maternally inherited, hence "call Mom". Script takes as input a reference sequence, and a fasta file with individual mtDNA sequences. The script identifies variants with respect to the reference and outputs in VCF format.

How to use:

./callMom.VERSION.py REFERENCE INDIVIDUALS

Where:

VERSION is the downloaded version of callMom

The REFERENCE is a single fasta file, such a the rCRS (Cambridge Reference Sequence)

The INDIVIDUALS file has a series of fasta files in the format:

\>sample id
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT

The sequence should be the same length as the reference, and all on a single line.
