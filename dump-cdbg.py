#! /usr/bin/env python3
from __future__ import print_function

import sys
import khmer
from sourmash_lib import MinHash
import screed
import argparse
from collections import OrderedDict, defaultdict
import os, os.path
from spacegraphcats import graph_parser
from cydata import intset, int_string_map, int_int_map


# graph settings
DEFAULT_KSIZE=31
DEFAULT_MEMORY = 1e8

# minhash settings
MH_SIZE_DIVISOR=50
MH_MIN_SIZE=5

class Pathfinder(object):
    "Track segment IDs, adjacency lists, and MinHashes"
    def __init__(self, ksize, mxtfile, node_offset=0):
        self.ksize = ksize

        self.node_counter = 1 + node_offset
        self.nodes = int_int_map()             # node IDs (int) to size
        self.nodes_to_kmers = int_int_map()    # node IDs (int) to kmers
        self.kmers_to_nodes = int_int_map()    # kmers to node IDs
        self.adjacencies = defaultdict(intset) # node to node
        self.labels = defaultdict(intset)      # nodes to set of labels
        self.mxtfp = open(mxtfile, 'wt')
        #self.assemblyfp = open(mxtfile + '.assembly', 'wt')

    def new_hdn(self, kmer):
        "Add a new high-degree node to the cDBG."
        if kmer in self.kmers_to_nodes:
            return self.kmers_to_nodes[kmer]

        this_id = self.node_counter
        self.node_counter += 1

        self.nodes[this_id] = 1
        self.nodes_to_kmers[this_id] = kmer
        self.kmers_to_nodes[kmer] = this_id

        return this_id

    def new_linear_node(self, visited, size):
        "Add a new linear path to the cDBG."
        node_id = self.node_counter
        self.node_counter += 1
        self.nodes[node_id] = size

        kmer = min(visited)               # identify linear nodes by min(hash)
        self.nodes_to_kmers[node_id] = kmer
        self.kmers_to_nodes[kmer] = node_id

        return node_id

    def add_adjacency(self, node_id, adj):
        "Add an edge between two nodes to the cDBG."
        node_id, adj = min(node_id, adj), max(node_id, adj)
        
        x = self.adjacencies[node_id]
        x.add(adj)

    def add_label(self, kmer, label):
        x = self.labels[kmer]
        x.add(label)

    def add_minhash(self, path_id, mh):
        # save minhash info to disk
        mins = " ".join(map(str, mh.get_mins()))
        self.mxtfp.write('{0},{1}\n'.format(path_id, mins))

    def add_path_assembly(self, path_id, assembly):
        self.assemblyfp.write('>{0}\n{1}\n'.format(path_id, assembly))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('output')
    args = p.parse_args()

    gxtfile = os.path.basename(args.output) + '.gxt'
    gxtfile = os.path.join(args.output, gxtfile)

    dumpfile = gxtfile + '.pickle'

    from pickle import load
    with open(dumpfile, 'rb') as fp:
        pathy = load(fp)

    print(len(pathy.nodes))
                    
    # save to GXT/MXT.
    print('saving gxtfile', gxtfile)

    all_labels = set()
    label_counts = {}
    with open(gxtfile, 'w') as fp:
        w = graph_parser.Writer(fp, ['labels'], [])

        for k, v in pathy.nodes.items():
            print('1', k, v)
            kmer = pathy.nodes_to_kmers.get(k)
            l = ""
            if kmer:
                labels = pathy.labels.get(kmer)
                if labels:
                    for x in labels:
                        label_counts[x] = label_counts.get(x, 0) + 1
                    all_labels.update(labels)
                    l = " ".join(map(str, labels))
            w.add_vertex(k, v, [l])

        for k, v in pathy.adjacencies.items():
            print('2', k, v)
            for edge in v:
                w.add_edge(k, edge, [])


if __name__ == '__main__':
    main()
