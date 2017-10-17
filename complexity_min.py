"""
Consolidate molecular barcodes unless a certain entropy threshold is met, then split molecular reads to try to
minimize complexity until entropy_min_reads or entropy_max is satisfied


"""
__author__ = 'Allison MacLeay'


import HTSeq
import sys
import os
import logging
import numpy as np
from collections import Counter
from scipy.spatial.distance import hamming
import itertools


# Configure logger
logging.basicConfig()
logger = logging.getLogger('complexity_min')

logger.level = logging.DEBUG

def read_bins(fastq_file):
    infile = HTSeq.FastqReader(fastq_file)
    read_num = 0
    cur_molecular_id = ''
    cur_sample_id = ''
    for read in infile:
        read_num += 1
        read_name, sample_id, molecular_id = read.name.split(' ')
        if molecular_id == cur_molecular_id:
            bin_reads.append(read)
        else:
            if cur_molecular_id != '':
                yield cur_molecular_id, cur_sample_id, bin_reads
            cur_molecular_id = molecular_id
            cur_sample_id = sample_id
            bin_reads = [read]
    yield cur_molecular_id, cur_sample_id, bin_reads  # The last bin


def get_base_counts(bases, quals, min_qual):
    num = {}
    qual = {}
    num['A'] = num['C'] = num['G'] = num['T'] = num['N'] = 0
    qual['A'] = qual['C'] = qual['G'] = qual['T'] = qual['N'] = 0
    for bb, qq in zip(bases, quals):
        if qq > min_qual:
            num[bb] += 1
        if qq > qual[bb]:
            qual[bb] = qq
    return num, qual


def consolidate_position(bases, quals, min_qual, min_freq):
    num, qual = get_base_counts(bases, quals, min_qual)
    most_common_base = max(num.iterkeys(), key=(lambda key: num[key]))
    freq = float(num[most_common_base]) / len(bases)
    entropy = calc_entropy([num[key] for key in num if key != 'N'])
    if freq >= min_freq:
        return True, most_common_base, qual[most_common_base], entropy
    else:
        return False, 'N', 0, .5


def calc_entropy(counts):
    """ Calculate entropy by position """
    counts = np.array(counts).astype(float)
    features = np.count_nonzero(counts)
    if features <= 1:
        return 0
    probs = counts/sum(counts)
    ent = 0.
    for i in [pr for pr in probs if pr != 0]:
        ent -= i * np.log(i)/np.log(features)
    return ent

class HaltingConditions(object):
    def __init__(self, reads, depth, limit):
        self.reads = reads
        self.depth = depth
        self.limit = limit
    def evaluate(self):
        raise NotImplementedError

class SequenceLimit(HaltingConditions):
    def __init__(self, reads, depth, limit=1):
        if limit is None:
            limit = 1
        super(SequenceLimit, self).__init__(reads, depth, limit)

    def evaluate(self):
        return len(set([read.seq for read in self.reads])) > self.limit

class DepthLimit(HaltingConditions):
    def __init__(self, reads, depth, limit=2):
        logger.debug('Level %d', depth)
        if limit is None:
            limit = 2
        super(DepthLimit, self).__init__(reads, depth, limit)

    def evaluate(self):
        return self.depth < self.limit and len(set([read.seq for read in self.reads])) > 1





class SequencePicker(object):
    def __init__(self, reads):
        self.reads = reads
        self.data = None
        self.name = None

    def evaluate(self):
        raise NotImplementedError


class MostCommon(SequencePicker):

    def __init__(self, reads):
        super(MostCommon, self).__init__(reads)
        self.name = 'most common 2 sequences'

    def pick(self):
        seqs = Counter([read.seq for read in self.reads])
        return [seq[0] for seq in seqs.most_common(2)]


class MostDifferent(SequencePicker):
    def __init__(self, reads):
        super(MostDifferent, self).__init__(reads)

    def pick(self):
        max_hamming = 0
        ma, mb = None, None
        for a, b in itertools.combinations(self.reads, 2):
            this_hamming = hamming([char for char in a.seq], [char for char in b.seq])
            if this_hamming > max_hamming or ma is None:
                max_hamming = this_hamming
                ma = a.seq
                mb = b.seq
        return [ma, mb]

class MinEntropy(SequencePicker):
    def __init__(self, reads):
        super(MinEntropy, self).__init__(reads)

    @staticmethod
    def _get_entropy_array(reads):
        read_bases = zip(*[list(read) for read in reads])
        read_quals = zip(*[[70] * len(read) for read in reads])
        # Iterate position by position
        bases, quals = zip(*[get_base_counts(bases, quals, 0) for bases, quals in zip(read_bases, read_quals)])
        entr_arr = []
        for num, qual in zip(bases, quals):
            counts = [num[key] for key in num if key != 'N']
            entr_arr.append(calc_entropy(counts))
        return entr_arr

    def pick(self):
        readset = list(set([read.seq for read in self.reads]))
        if len(readset) == 2:
            return readset[0], readset[1]
        entr_arr = self._get_entropy_array([read.seq for read in self.reads])
        combos = {}
        for a, b in itertools.combinations(self.reads, 2):
            if a.seq == b.seq:
                continue
            reads0, reads1 = split_reads(self.reads, [a.seq, b.seq])
            combos[(a.seq, b.seq)] = [sum(self._get_entropy_array([read.seq for read in reads0])), sum(self._get_entropy_array([read.seq for read in reads1]))]

        # abs(len(reads1) - (len(reads0) + len(reads1))/2)
        max_hamming = 0
        ma, mb = None, None
        for a, b in itertools.combinations(self.reads, 2):
            this_hamming = hamming([char for char in a.seq], [char for char in b.seq])
            if this_hamming > max_hamming or ma is None:
                max_hamming = this_hamming
                ma = a.seq
                mb = b.seq
        return [ma, mb]


def split_reads(reads, new_seqs):
    reads0 = []
    reads1 = []
    for read in reads:
        read_chars = [char for char in read.seq]
        ham0 = hamming([char for char in new_seqs[0]], read_chars)
        ham1 = hamming([char for char in new_seqs[1]], read_chars)
        if ham0 <= ham1:
            reads0.append(read)
        else:
            reads1.append(read)

    return reads0, reads1

def minimize_entropy(bins, halting_condition=SequenceLimit, sequence_picker=MostCommon, depth=0, limit=None):
    new_bins = []
    for cur_mol_id, cur_sample_id, reads in bins:
        if halting_condition(reads, depth, limit).evaluate():
            reads0, reads1 = split_reads(reads, sequence_picker(reads).pick())
            bin0 = (cur_mol_id, cur_sample_id, reads0)
            bin1 = (cur_mol_id, cur_sample_id, reads1)
            new_bins.extend(minimize_entropy([bin0, bin1], halting_condition, sequence_picker, depth=depth+1))
        else:
            new_bins.extend([(cur_mol_id, cur_sample_id, reads)])
    return new_bins

def consolidate(r1_fastq, r2_fastq, r1_consolidated_fastq, r2_consolidated_fastq, min_qual, min_freq,
                min_entropy_freq, max_entropy):
    outfolder = os.path.dirname(r1_consolidated_fastq)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    outfile = open(r1_consolidated_fastq, 'w')
    bins = list(read_bins(r1_fastq))
    num_bins_starting = len(bins)
    bins = minimize_entropy(bins, DepthLimit, MinEntropy)

    num_input_reads = 0
    num_consolidated_reads = 0
    num_successes = 0  # Bases with successful consolidation
    num_bases = 0
    for cur_molecular_id, cur_sample_id, reads in bins:
        num_input_reads += len(reads)
        num_consolidated_reads += 1
        # Get all the bases and quals in the read
        read_bases = zip(*[list(read.seq) for read in reads])
        read_quals = zip(*[list(read.qual) for read in reads])
        # Iterate position by position
        consolidation_sucess, cons_seq, cons_qual, cons_entropy = zip(*[consolidate_position(bases, quals, min_qual, min_freq) for bases, quals in zip(read_bases, read_quals)])
        # Count consolidation successes and failures
        num_successes += sum(consolidation_sucess)
        num_bases += len(consolidation_sucess)
        # Write consolidated FASTQ read
        outfile.write('@%s_%d %s\n' % (cur_molecular_id, len(reads), cur_sample_id)) # Header: Molecular id, number of reads, 2nd incoming header field (includes sample id)
        outfile.write(''.join(cons_seq) +'\n')
        outfile.write('+\n')
        outfile.write(''.join([chr(q+33) for q in cons_qual]) + '\n')


    num_bins_end = len(bins)
    logger.info("Created %d bins from %d molecular barcodes", num_bins_end, num_bins_starting)
    logger.info("Read %d input reads", num_input_reads)
    logger.info("Wrote %d consolidated reads", num_consolidated_reads)
    logger.info("Successfully consolidated %d bases out of %d (%.2f%%)", num_successes, num_bases, 100*float(num_successes)/num_bases)
    outfile.close()


def main():
    if len(sys.argv) < 5:
        print 'Usage: python consolidate.py r1_fastq r2_fastq r1_consolidated_fastq r2_consolidated_fastq min_qual min_freq min_entropy_freq max_entropy'
        sys.exit()

    r1_fastq_file = sys.argv[1]
    r1_consolidated_fastq_file = sys.argv[3]
    r2_fastq_file = sys.argv[2]
    r2_consolidated_fastq_file = sys.argv[4]
    min_qual = int(sys.argv[5])
    min_freq = float(sys.argv[6])
    min_entropy_freq = int(sys.argv[7])
    max_entropy = float(sys.argv[8])
    consolidate(r1_fastq_file, r2_fastq_file, r1_consolidated_fastq_file, r2_consolidated_fastq_file,
                min_qual, min_freq, min_entropy_freq, max_entropy)


if __name__ == '__main__':
    main()
