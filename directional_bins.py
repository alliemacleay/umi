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
logger = logging.getLogger('directional_bins')

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


def get_base_counts(bases, quals, min_qual, read_counts):
    num = {}
    qual = {}
    num['A'] = num['C'] = num['G'] = num['T'] = num['N'] = 0
    qual['A'] = qual['C'] = qual['G'] = qual['T'] = qual['N'] = 0
    for bb, qq, count in zip(bases, quals, read_counts):
        if qq > min_qual:
            num[bb] += (1 * count)
        if qq > qual[bb]:
            qual[bb] = qq

    return num, qual


def consolidate_position(bases, quals, read_counts, min_qual, min_freq):
    num, qual = get_base_counts(bases, quals, min_qual, read_counts)
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


class Node(object):
    def __init__(self, seq, count, quals, mol_id, sample_id, mismatches_allowed=1):
        self.seq = seq
        self.count = count
        self.quals = quals
        self.molecule_id = mol_id
        self.sample_id = sample_id
        self.children = []
        self.mismatches_allowed = mismatches_allowed

    def add_child(self, node):
        self.children.append(node)

    def is_connected(self, node):
        if self.distance(node) < (self.mismatches_allowed + 1) and self.counts_proportional(node):
            return True
        return False

    def counts_proportional(self, node):
        counta = max(node.count, self.count)
        countb = min(node.count, self.count)
        return counta >= 2 * countb - 1

    def is_parent_of(self, node):
        return self.count > node.count or (node.count == 1 and self.count == 1)

    def distance(self, node):
        differences = 0
        for a, b in zip([basea for basea in node.seq], [baseb for baseb in self.seq]):
            if a != b:
                differences += 1
        return differences

    def show_distance(self, node):
        """ self: atCg
            node: atTg
        """
        seqa = []
        seqb = []
        for b, a in zip([basea for basea in node.seq], [baseb for baseb in self.seq]):
            if a != b:
                seqa.append(a.upper())
                seqb.append(b.upper())
            else:
                seqa.append(a.lower())
                seqb.append(b.lower())
        print('{}\n{}'.format(''.join(seqa), ''.join(seqb)))



class Bin(object):
    def __init__(self, node):
        self.head_node = node
        self.nodes = [self.head_node]

    def __len__(self):
        return len(self.nodes)

    def get_seq(self, min_qual=0, min_freq=.001):
        max_counts = self.head_node.count
        for child in self.head_node.children:
            if child.count > max_counts:
                raise ValueError('Children have greater counts than parent')
            elif child.count == max_counts:
                if max_counts == 1 and len(self.nodes) > 1:
                    read_bases = zip(*[list(read.seq) for read in self.nodes])
                    read_quals = zip(*[list(read.quals) for read in self.nodes])
                    read_counts = [read.count for read in self.nodes]
                    # Iterate position by position
                    consolidation_sucess, cons_seq, cons_qual, cons_entropy = zip(
                        *[consolidate_position(bases, quals, read_counts, min_qual, min_freq) for bases, quals in
                          zip(read_bases, read_quals)])
                    return ''.join(cons_seq), cons_qual, cons_entropy
                raise ValueError('Can not choose consensus. {} == {}\n{}'.format(child.count, max_counts, child.show_distance(self.head_node)))
        return self.head_node.seq, self.head_node.quals, 0.

    def accepts(self, node):
        for existing_node in self.nodes:
            if existing_node.is_connected(node):
                return True
        return False

    def get_reads(self):
        return [node.seq for node in self.nodes]

    def add(self, node):
        for existing_node in self.nodes:
            if existing_node.is_connected(node):
                if existing_node.is_parent_of(node):
                    existing_node.add_child(node)
                    self.nodes.append(node)
                else:
                    raise ValueError('Check for adding nodes in descending order of counts')
                return

    def connect(self, bin, node):
        if self.head_node.count < bin.head_node.count:
            raise ValueError('Must merge lower bin with lower head node into self')

        found_node = None
        if not self.node_exists(node):
            self.add(node)
        for mynode in self.nodes:
            if mynode.seq == node.seq:
                found_node = mynode
                break
        for binnode in bin.nodes:
            if binnode.distance(found_node) == 1:
                found_node.add_child(binnode)
                break
        self.nodes.extend(bin.nodes)
        try:
            self.get_seq()
        except ValueError as e:
            pass

    def node_exists(self, node):
        return len([mynode.seq for mynode in self.nodes if mynode.seq == node.seq]) > 0



def directional_bins(bins, mismatches_allowed=1):
    new_bins = []
    merged = 0
    for mol_id, sample_id, reads in bins:
        mol_bins = []
        seqs = [read.seq for read in reads]
        quals = {read.seq: read.qual for read in reads}
        counts = Counter(seqs)
        for seq, count in counts.most_common():
            this_node = Node(seq, count, quals[seq], mol_id, sample_id, mismatches_allowed)
            added_to = None
            offset = 0
            for i in range(len(mol_bins)):
                bin = mol_bins[i - offset]
                if bin.accepts(this_node):
                    if added_to is not None:
                        merged += 1
                        if added_to.head_node.count > bin.head_node.count:
                            mol_bins.remove(bin)
                            added_to.connect(bin, this_node)
                            offset += 1
                        else:
                            mol_bins.remove(added_to)
                            bin.connect(added_to, this_node)
                            added_to = bin
                            offset += 1
                    else:
                        bin.add(this_node)
                        added_to = bin
            if added_to is None:
                mol_bins.append(Bin(this_node))
        new_bins.extend(mol_bins)
    logger.debug('Merged %d', merged)
    return new_bins


def directional_consolidation(r1_fastq, r2_fastq, r1_consolidated_fastq, r2_consolidated_fastq, min_qual, min_freq,
                              mismatches_allowed=2):
    outfolder = os.path.dirname(r1_consolidated_fastq)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    outfile = open(r1_consolidated_fastq, 'w')
    mol_bins = list(read_bins(r1_fastq))
    num_bins_starting = len(mol_bins)
    bins = directional_bins(mol_bins, mismatches_allowed=mismatches_allowed)
    for bin in bins:
        node = bin.head_node
        cur_molecular_id = node.molecule_id
        cur_sample_id = node.sample_id
        cons_seq, cons_qual, cons_entropy = bin.get_seq()
        outfile.write('@%s_%d %s\n' % (cur_molecular_id, len(bin), cur_sample_id)) # Header: Molecular id, number of reads, 2nd incoming header field (includes sample id)
        outfile.write(''.join(cons_seq) +'\n')
        outfile.write('+\n')
        outfile.write(''.join([chr(q+33) for q in cons_qual]) + '\n')
    outfile.close()
    num_bins_end = len(bins)
    logger.info("Created %d bins from %d molecular barcodes", num_bins_end, num_bins_starting)
    return num_bins_end

def main():
    if len(sys.argv) < 5:
        print 'Usage: python directional_bins.py r1_fastq r2_fastq r1_consolidated_fastq r2_consolidated_fastq min_qual min_freq'
        sys.exit()

    r1_fastq_file = sys.argv[1]
    r1_consolidated_fastq_file = sys.argv[3]
    r2_fastq_file = sys.argv[2]
    r2_consolidated_fastq_file = sys.argv[4]
    min_qual = int(sys.argv[5])
    min_freq = float(sys.argv[6])
    directional_consolidation(r1_fastq_file, r2_fastq_file, r1_consolidated_fastq_file, r2_consolidated_fastq_file,
                min_qual, min_freq)


#if __name__ == '__main__':
#    main()

f1='test/data/complex/r1_out.fastq'
f2='test/data/complex/r2_out.fastq'
o1='test/data/complex/r1_cons.fastq'
o2='test/data/complex/r2_cons.fastq'

res = {}
for i in range(60):
    res[i + 1] = directional_consolidation(f1, f2, o1, o2, 0, .6, i+1)
print res
import matplotlib.pyplot as plt
plt.plot(res.keys(), res.values())
plt.ylabel('bins')
plt.xlabel('mismatches allowed')
plt.show()
plt.savefig('/Users/Admin/Documents/mgh/projects/umi_data/bins_by_mismatches.png')


