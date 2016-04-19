#!/usr/local/bin/python3
import sys
import re
import os
import copy
import pickle

if len(sys.argv) < 2:
    sys.exit('./callMom.py <reference> <fasta>')


### CLASS DEFINITIONS
# Site variables: id, pos, end, ref, alt, len, out, vtype, maxlen
# create site class
class Site():
    def __init__(self, a, b, c, d):
        self.id = a
        self.pos = int(b)
        self.end = int(b)
        self.ref = c
        self.alt = d
        self.len = 1
        self.out = 0
        self.vtype = ''
        self.maxlen = 1

    def printsite(self):
        self.updatesite()
        out = self.id + '\t' + str(self.pos) + '\t' + str(self.end) + '\t' + self.ref + '\t' + self.alt + '\t' + str(
            self.len) + '\t' + str(self.maxlen) + '\t' + self.vtype + '\t' + str(self.out)
        return out

    def updatesite(self):
        self.len = len(self.alt)
        self.end = self.pos + self.len -1


# create object for storing all variant alleles at a position in refseq
class Var():
    def __init__(self, a):
        self.pos = int(a)
        self.len = []
        self.ref = []
        self.alt = []
        self.altcode = []
        self.ac = []
        self.vtype = []
        self.maxlen = 0
        self.hasvar = {}

    # extend ref/alt by 1 bp upstream
    def addbp(self, r):
        for i in range(len(self.ref)):
            self.ref[i] = r[self.pos - 1] + self.ref[i]
        for i in range(len(self.alt)):
            self.alt[i] = r[self.pos - 1] + self.alt[i]
            self.len[i] = len(self.alt[i])
            if self.len[i] > self.maxlen:
                self.maxlen = self.len[i]
        self.vtype = 'M'

    # function for outputting a string of the object contents
    # varfrqout1.write("pos\tref\talt\tlen\tac\tvtype\thas\n")
    def printvar(self):
        out = str(self.pos)
        tref = []
        for a1 in self.ref:
            tref.append(a1)
        out = out + '\t' + ','.join(tref)
        talt = []
        for a2 in self.alt:
            talt.append(a2)
        out = out + '\t' + ','.join(talt)
        tlen = []
        for a3 in self.len:
            tlen.append(str(a3))
        out = out + '\t' + ','.join(tlen)
        tac = []
        for a4 in self.ac:
            tac.append(str(a4))
        out = out + '\t' + ','.join(tac)
        tvtype = []
        for a5 in self.vtype:
            tvtype.append(str(a5))
        out = out + '\t' + ','.join(tvtype)
        thasvar = []
        done = []
        for a6 in range(len(self.ref)):
            for a7 in range(len(self.alt)):
                k = self.ref[a6] + ':' + self.alt[a7]
                if k in self.hasvar.keys() and k not in done:
                    thasvar.append(self.ref[a6] + ':' + self.alt[a7] + ':' + ','.join(
                        self.hasvar[self.ref[a6] + ':' + self.alt[a7]]))
                    done.append(k)
        out = out + '\t' + '\t'.join(thasvar)
        return out

    def vcf(v, ids):
        # combine genotypes into a string
        gtar = []
        # loop through ids
        for i in ids:
            hasalt = 0
            # loop through alts
            for a in range(len(v.alt)):
                k = v.ref[0] + ':' + v.alt[a]
                if i in v.hasvar[k]:
                    # individual has this alt allele
                    # add genotype code
                    gtcode = a + 1
                    gtar.append(str(gtcode))
                    hasalt = 1
            if hasalt == 0:
                gtar.append('0')
        gtstr = '\t'.join(gtar)
        # put together string of ac
        acar = []
        for a in v.ac:
            acar.append(str(a))
        acstr = ','.join(acar)
        # compile final output string
        out = '\t'.join(['MT', str(v.pos), '.', v.ref[0], ','.join(v.alt), '100', 'fa',
                         'VT=' + ','.join(v.vtype) + ';' + 'AC=' + acstr, 'GT', gtstr])
        return out

    # check for errors
    def err(v1):
        if len(set(v1.alt)) != len(v1.alt) or len(set(v1.ref)) != len(v1.ref):
            return 1
        if len(set(v1.alt)) == len(v1.alt) and len(set(v1.ref)) == len(v1.ref):
            return 0

    # check if two variants 1 bp away are the same indel
    def samevar(v1, v2):
        for vt1 in range(len(v1.vtype)):
            for vt2 in range(len(v2.vtype)):
                if v1.vtype[vt1] == 'I' and v2.vtype[vt2] == 'I':
                    return 1
                if v1.vtype[vt1] == 'M' and v2.vtype[vt2] == 'M':
                    return 1

    # lengthen snp alt when overlapping indel
    def updatealt(self):
        if len(self.alt) > 1:
            # multiple alternate alleles
            if 'I' not in set(self.vtype) and 'M' not in set(self.vtype):
                # multiple SNP alleles
                self.ref = list(set(self.ref))
            if 'I' in set(self.vtype) or 'M' in set(self.vtype):
                # contains complex variants
                # get longest ref
                maxlen = 0
                maxref = ''
                for r in self.ref:
                    if len(r) > maxlen:
                        maxref = r
                        maxlen = len(r)
                self.maxlen = maxlen
                # store the original variant
                ogv = copy.deepcopy(self)
                # set the new ref
                self.ref = []
                self.ref.append(maxref)
                # loop through alts
                for i in range(len(self.alt)):
                    # add bases to the alt for snps
                    if 'S' in self.vtype[i]:
                        j = len(self.alt[i])
                        while len(self.alt[i]) < len(self.ref[0]):
                            self.alt[i] = self.alt[i] + list(self.ref[0])[j]
                            j = j + 1
                    # add bases to the alt for mnps
                    if 'M' in self.vtype[i]:
                        j = len(self.alt[i])
                        while len(self.alt[i]) < len(self.ref[0]):
                            self.alt[i] = self.alt[i] + list(self.ref[0])[j]
                            j = j + 1
                    # add bases to the alt for indels
                    if 'I' in self.vtype[i]:
                        # no need to edit when ref len is same
                        if len(self.ref[0]) != len(ogv.ref[i]):
                            # add bases to alt
                            j = len(ogv.ref[i])
                            while j < len(self.ref[0]):
                                self.alt[i] = self.alt[i] + list(self.ref[0])[j]
                                j = j + 1
                                # update all the hasvar keys
                temp = {}
                for i in range(len(ogv.alt)):
                    k0 = ogv.ref[i] + ':' + ogv.alt[i]
                    k1 = self.ref[0] + ':' + self.alt[i]
                    if k1 not in temp.keys():
                        temp[k1] = self.hasvar[k0]
                self.hasvar = temp

    # sort alt alleles in order of frequency
    # assign numbers to the alleles
    # self.pos=int(a)
    # self.len=[]
    # self.ref=[]
    # self.alt=[]
    # self.altcode=[]
    # self.ac=[]
    # self.vtype=[]
    # self.maxlen=0
    # self.hasvar={}
    def sortalt(self):
        # create a duplicate variant
        ogv = copy.deepcopy(self)
        # only needed if n>1
        if len(self.ac) > 1:
            # store the unsorted list
            usoac = copy.deepcopy(self.ac)
            # sort the list
            soac = copy.deepcopy(self.ac)
            soac = sorted(soac, reverse=True)
            # update the sorted list in self
            self.ac.sort(reverse=True)
            # reset the other items to be sorted
            self.alt = []
            self.len = []
            self.vtype = []
            # loop through new ac order
            for i in range(len(soac)):
                # pop the first item
                ac = soac.pop(0)
                # get first index of item with ac
                ui = usoac.index(ac)
                # clear the value
                usoac[ui] = ''
                # update other lists
                self.alt.append(ogv.alt[ui])
                self.len.append(ogv.len[ui])
                self.vtype.append(ogv.vtype[ui])


### OPERATION CODE

# STEP 1
# read in reference sequence to list
refseq = ['>']
f1 = open(sys.argv[1], 'r')
for l in f1:
    l = l.rstrip()
    if '>' not in l:
        for b in l:
            refseq.append(b.upper())
f1.close()

bases = ['A', 'C', 'G', 'T']
# STEP 2
# identify all variant bases
# store in database indexed by id then position
initvar = open('cmo.0.initvar.txt', 'w')
ids = []
vars = {}
curid = ''
f2 = open(sys.argv[2], 'r')
acgtn = ['A', 'C', 'G', 'T', 'N']
indel = ['-']
fastas = {}
for l in f2:
    l = l.rstrip()
    # read the sequence
    if '>' not in l:
        seq = '>' + l
        fastas[curid]=seq
        # loop through fasta file
        # check if variant in individual
        p = 0
        while p < len(seq):
            # matches reference, continue
            if refseq[p] == seq[p].upper():
                p = p + 1
                continue
            # does not match reference, variant
            if refseq[p] != seq[p].upper():
                # record variant site
                initvar.write(curid + '\t' + str(p) + '\n')
                # is it a SNP or MNP?
                if seq[p].upper() in acgtn:
                    # initalize variant record for id at position
                    vars[curid][p] = ''
                    q = p
                    # check all consecutive positions that differ
                    done = 0
                    while refseq[q] != seq[q].upper() and done == 0:
                        if seq[q].upper() in acgtn:
                            # store the variant position and base
                            vars[curid][p] = vars[curid][p] + seq[q].upper()
                            # go to the next base
                            q = q + 1
                        if seq[q].upper() not in acgtn:
                            done = 1
                    p = q
                # found a deletion
                if seq[p] in indel:
                    vars[curid][p] = ''
                    q = p
                    # check all consecutive positions that differ
                    done = 0
                    while refseq[q] != seq[q].upper() and done == 0:
                        if seq[q].upper() in indel:
                            # store the variant position and base
                            vars[curid][p] = vars[curid][p] + seq[q].upper()
                            q = q + 1
                        if seq[q].upper() not in indel:
                            done = 1
                    p = q
                if seq[p] not in indel and seq[p].upper() not in acgtn:
                    # skip these sites
                    p = p + 1
    # read the id
    if '>' in l:
        l = l.replace('>', '')
        ids.append(l)
        # initialize database of variants
        vars[l] = {}
        curid = l
f2.close()
initvar.close()

# DATA 
# vars[id][position] contains list of variants for each individual 

# STEP 3
# populate list of  all individual variant alleles
sites = []
# loop through ids
for i in ids:
    # loop through variant sites in this individual
    for j in sorted(vars[i].keys()):
        sites.append(Site(i, j, ''.join(refseq[j:j + len(vars[i][j])]), vars[i][j]))

siteout0 = open('cmo.1.sites.0.txt', 'w')
for s in range(len(sites)):
    if sites[s].out == 0:
        siteout0.write(sites[s].printsite() + '\n')
siteout0.close()

## identify long variants
longsites=[]
longsitebp=[]
maxlen=0
for s in range(len(sites)):
    if len(sites[s].ref) > 1:
        longsites.append(sites[s])
        for i in range(sites[s].pos,sites[s].end+1):
                longsitebp.append(i)
        if len(sites[s].ref) > maxlen:
                maxlen=len(sites[s].ref)


## extend variants overlaping long variants
for s in range(len(sites)):
    if sites[s].pos in longsitebp:
        if s+1 < len(sites):
            if sites[s+1].pos in longsitebp:
                if sites[s].id == sites[s+1].id:
                    if maxlen > (sites[s+1].pos - sites[s].pos):
                        #print(sites[s].printsite()+' AND '+sites[s+1].printsite())
                        #sites.pop(s+1)
                        newref=''
                        newalt=''
                        for j in range(sites[s].pos,sites[s+1].end+1):
                            newref=newref+refseq[j]
                        for j in range(sites[s].pos,sites[s+1].end+1):
                            newalt=newalt+fastas[sites[s].id][j]
                        sites[s].ref=newref
                        sites[s].alt=newalt
                        sites[s].end=sites[s+1].end
                     #   print(sites[s].printsite()+' AND '+sites[s+1].printsite())
                        sites[s+1].out=1

siteout0.close()

for s in sites:
    if s.out==1:
        sites.remove(s)
#   out = self.id + '\t' + str(self.pos) + '\t' + str(self.end) + '\t' + self.ref + '\t' + self.alt + '\t' + str(
#    self.len) + '\t' + str(self.maxlen) + '\t' + self.vtype + '\t' + str(self.out)

# Site variables: id, pos, end, ref, alt, len, out, vtype, maxlen

# STEP 4
# shift indels over by 1 bp and extend reference
# also, identify variant type
# also, update mnp and indel length
for s in range(len(sites)):
    # check if indel
    if '-' not in sites[s].alt:
        if len(sites[s].ref) == 1 and len(sites[s].alt) == 1:
            sites[s].vtype = 'S'
        if len(sites[s].ref) > 1 and len(sites[s].alt) >= 1:
            sites[s].vtype = 'M'
            sites[s].updatesite()
    if '-' in sites[s].alt:
        # set variant type
        sites[s].vtype = 'I'
        # set position over 1 bp for indel
        sites[s].pos = sites[s].pos - 1
        # get length of the indel
        sites[s].len = len(sites[s].alt)
        # for deletion, alt is single base
        sites[s].alt = refseq[sites[s].pos]
        sites[s].ref = ''
        p = sites[s].pos
        # ref is length of alt + 1
        while len(sites[s].ref) < sites[s].len + 1:
            sites[s].ref = sites[s].ref + refseq[p]
            p = p + 1
        sites[s].updatesite()


# DATA
# vars[id][pos] list of variants per individual
# sites[s] list of individual variant sites

# STEP 5
# print site list after indel position and mnp length update
siteout1 = open('cmo.2.sites.1.txt', 'w')
for s in range(len(sites)):
    if sites[s].out == 0:
        siteout1.write(sites[s].printsite() + '\n')
siteout1.close()

# STEP 6
# calculate frequency of all ref/alt alleles at a site
varfrq = {}
site2var = open('cmo.3.site2var.txt', 'w')
for s in range(len(sites)):
    if sites[s].out == 0:
        # existing site, compare ref/alt
        if sites[s].pos in varfrq.keys():
            # does the ref allele match?
            if sites[s].ref in varfrq[sites[s].pos].ref:
                # does the alt allele match?
                if sites[s].alt in varfrq[sites[s].pos].alt:
                    # exact match
                    # only one ref/alt, increment ac
                    if len(varfrq[sites[s].pos].ref) == 1 and len(varfrq[sites[s].pos].alt) == 1:
                        # increment variant count
                        varfrq[sites[s].pos].ac[0] = 1 + varfrq[sites[s].pos].ac[0]
                        # append to list of individuals with the variant
                        varfrq[sites[s].pos].hasvar[sites[s].ref + ':' + sites[s].alt].append(sites[s].id)
                    # multiple ref/alt alleles
                    if len(varfrq[sites[s].pos].ref) > 1 or len(varfrq[sites[s].pos].alt) > 1:
                        for r in range(len(varfrq[sites[s].pos].ref)):
                            a = r
                            # check the list for ref and alt
                            if sites[s].ref in varfrq[sites[s].pos].ref and sites[s].alt in varfrq[sites[s].pos].alt:
                                # loop through list looking for the ref/alt
                                if varfrq[sites[s].pos].ref[r] == sites[s].ref and varfrq[sites[s].pos].alt[a] == sites[
                                    s].alt:
                                    # found a match in the list, check that same ref/alt combination
                                    # increment ac for that ref/alt combination
                                    varfrq[sites[s].pos].ac[a] = 1 + varfrq[sites[s].pos].ac[a]
                                    # increment list of individuals with this ref/alt conmbination
                                    varfrq[sites[s].pos].hasvar[sites[s].ref + ':' + sites[s].alt].append(sites[s].id)
                                    site2var.write(varfrq[sites[s].pos].printvar() + '\n')
                            if sites[s].ref not in varfrq[sites[s].pos].ref or sites[s].alt not in varfrq[
                                sites[s].pos].alt:
                                # new ref/alt combination to be added
                                varfrq[sites[s].pos].len.append(sites[s].len)
                                varfrq[sites[s].pos].ref.append(sites[s].ref)
                                varfrq[sites[s].pos].alt.append(sites[s].alt)
                                varfrq[sites[s].pos].ac.append(1)
                                varfrq[sites[s].pos].hasvar[sites[s].ref + ':' + sites[s].alt].append(sites[s].id)
                                varfrq[sites[s].pos].vtype.append(sites[s].vtype)
                                site2var.write(varfrq[sites[s].pos].printvar() + '\n')
                # same ref, diff alt, add to pos
                if sites[s].alt not in varfrq[sites[s].pos].alt:
                    varfrq[sites[s].pos].len.append(sites[s].len)
                    varfrq[sites[s].pos].ref.append(sites[s].ref)
                    varfrq[sites[s].pos].alt.append(sites[s].alt)
                    varfrq[sites[s].pos].ac.append(1)
                    varfrq[sites[s].pos].hasvar[sites[s].ref + ':' + sites[s].alt] = [sites[s].id]
                    varfrq[sites[s].pos].vtype.append(sites[s].vtype)
                    site2var.write(varfrq[sites[s].pos].printvar() + '\n')
            # novel ref allele at pos
            if sites[s].ref not in varfrq[sites[s].pos].ref:
                varfrq[sites[s].pos].len.append(sites[s].len)
                varfrq[sites[s].pos].ref.append(sites[s].ref)
                varfrq[sites[s].pos].alt.append(sites[s].alt)
                varfrq[sites[s].pos].ac.append(1)
                varfrq[sites[s].pos].hasvar[sites[s].ref + ':' + sites[s].alt] = [sites[s].id]
                varfrq[sites[s].pos].vtype.append(sites[s].vtype)
                site2var.write(varfrq[sites[s].pos].printvar() + '\n')
        # new site, creat variant
        if sites[s].pos not in varfrq.keys():
            varfrq[sites[s].pos] = Var(sites[s].pos)
            varfrq[sites[s].pos].len.append(sites[s].len)
            varfrq[sites[s].pos].ref.append(sites[s].ref)
            varfrq[sites[s].pos].alt.append(sites[s].alt)
            varfrq[sites[s].pos].ac.append(1)
            varfrq[sites[s].pos].hasvar[sites[s].ref + ':' + sites[s].alt] = [sites[s].id]
            varfrq[sites[s].pos].vtype.append(sites[s].vtype)
            site2var.write(varfrq[sites[s].pos].printvar() + '\n')
site2var.close()

# STEP 7
# Output initial set of variants
dropl = []
mergel1 = []
mergel2 = []
varfrqout0 = open('cmo.4.varfrq.0.txt', 'w')
varfrqout0.write("pos\tref\talt\tlen\tac\tvtype\thas\n")
for p in sorted(varfrq.keys()):
    varfrqout0.write(varfrq[p].printvar() + '\n')
varfrqout0.close()

# STEP 8
# updade variant reference and alternate alleles
varfrqout1 = open('cmo.5.varfrq.1.txt', 'w')
varfrqout1.write("pos\tref\talt\tlen\tac\tvtype\thas\n")
for p in sorted(varfrq.keys()):
    varfrq[p].updatealt()
    varfrqout1.write(varfrq[p].printvar() + '\n')
varfrqout1.close()

# STEP 9
# identify consecutive variants
adjacentsites = []
mergeout1 = open('cmo.6.merge.1.txt', 'w')
mergeout1.write('bp1\tbp2\n')
# identify consecutive variants 
for p in sorted(varfrq.keys()):
    if p + 1 in varfrq.keys():
        mergel1.append(str(p) + '\t' + str(p + 1))
        mergeout1.write(str(p) + '\t' + str(p + 1) + '\n')
        adjacentsites.append(str(p))
        adjacentsites.append(str(p + 1))
mergeout1.close()

disite = []
trisite = []

# STEP 10
# identify consecutive variants, round 2
mergeout2 = open('cmo.7.merge.2.txt', 'w')
mergeout2.write('bp1\tbp2\tbp3\n')
for i in range(len(mergel1)):
    if i + 1 < len(mergel1):
        a, b = mergel1[i].split('\t')
        c, d = mergel1[i + 1].split('\t')
        if int(b) == int(c):
            mergel2.append(a + '\t' + b + '\t' + d)
            mergeout2.write(a + '\t' + b + '\t' + d + '\n')
            trisite.append(a)
            trisite.append(b)
            trisite.append(d)
        if int(b) != int(c):
            mergel2.append(a + '\t' + b + '\tN')
            mergeout2.write(a + '\t' + b + '\tN\n')
            disite.append(a)
            disite.append(b)
mergeout2.close()

# STEP 10.1
# remove trisite from disite
for t in trisite:
    if t in disite:
        disite.remove(t)

# STEP 11
# Merge adjacent SNPs into MNPs

# DATA
# vars[id][pos] list of variants per individual
# sites[s] list of individual variant sites
# varfrq [pos] list of alleles and individuals having them

# STEP 12
# sort alternate alleles in decreasing prevalence
varfrqout3 = open('cmo.10.varfrq.3.txt', 'w')
varfrqout3.write("pos\tref\talt\tlen\tac\tvtype\thas\n")
for p in sorted(varfrq.keys()):
    varfrq[p].sortalt()
    varfrqout3.write(varfrq[p].printvar() + '\n')
varfrqout3.close()

# STEP 13
# VCF output
vcfout = open('ALL.chrMT.phase3_callmom.20130502.genotypes.vcf', 'w')
vcfout.write('##fileformat=VCFv4.2\n')
vcfout.write('##fileDate=20151002\n')
vcfout.write('##source=callMomV0.2\n')
vcfout.write('##reference=gi|251831106|ref|NC_012920.1| Homo sapiens mitochondrion, complete genome\n')
vcfout.write('##contig=<ID=MT,length=16569,assembly=b37>\n')
vcfout.write('##INFO=<ID=VT,Number=.,Type=String,Description="Alternate allele type. S=SNP, M=MNP, I=Indel">\n')
vcfout.write(
    '##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate allele counts, comma delimited when multiple">\n')
vcfout.write('##FILTER=<ID=fa,Description="Genotypes called from fasta file">\n')
vcfout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
vcfout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ids) + '\n')

for p in sorted(varfrq.keys()):
    #if p not in disite and p not in trisite:
    vcfout.write(Var.vcf(varfrq[p], ids) + '\n')
    #if str(p) in disite:
    #    vcfout.write(Var.vcf(varfrq[p], ids) + '\n')
    #if str(p) in trisite:
    #    vcfout.write(Var.vcf(varfrq[p], ids) + '\n')
vcfout.close()
