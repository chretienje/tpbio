#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file,'rt') as f:
        sentences = f.readlines()
        for sentence in sentences:
            sentence = sentence.replace("\n", "")
            for letter in sentence:
                if letter not in "TGCA":
                    break
            else:
                yield sentence

def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size+1):
        yield read[i : i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for sentence in read_fastq(fastq_file):
        for kmer in cut_kmer(sentence, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer in kmer_dict.keys():
        graph.add_weighted_edges_from([(kmer[:-1], kmer[1:], kmer_dict[kmer])])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    if delete_entry_node:
        for path in path_list:
            starting_node = path[0]
            graph.remove_node(starting_node)
    if delete_sink_node:
        for path in path_list:
            sink_node = path[-1]
            graph.remove_node(sink_node)
    for path in path_list:
        for node in path[1:-1]:
            graph.remove_node(node)
    return graph


def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):

    same_weight_path_index = set()
    weight_max = max(weight_avg_list)
    is_weight_unique = 0
    for i in range(len(weight_avg_list)):
        weight = weight_avg_list[i]
        if weight == weight_max:
            is_weight_unique += 1
            same_weight_path_index.add(i)
    if is_weight_unique == 1:
        indice = weight_avg_list.index(weight_max)
        graph = remove_paths(graph, path_list[:indice]+path_list[indice+1:], delete_entry_node, delete_sink_node)
        return graph
    same_length_path_index = set()
    length_max = max(path_length)
    is_length_unique = 0
    for i in range(len(path_length)):
        if i not in same_weight_path_index:
            continue
        length = path_length[i]
        if length == length_max:
            is_length_unique += 1
            same_length_path_index.add(i)
    if is_length_unique == 1:
        indice = path_length.index(length_max)
        graph = remove_paths(graph, path_list[:indice]+path_list[indice+1:], delete_entry_node, delete_sink_node)
        return graph

    path_list_remaining = [path[i] for i in range(len(path_list)) if i in same_length_path_index]
    indice = randint(0, len(path_list_remaining))
    graph = remove_paths(graph, path_list_remaining[:indice]+path_list_remaining[indice+1:], delete_entry_node, delete_sink_node)
    return graph



def path_average_weight(graph, path):
    s = 0
    for i in range(len(path)-1):
        node = path[i]
        boo = False
        for n, nbrs in graph.adj.items():
            if boo:
                break
            if n == node:
                for nbr, eattr in nbrs.items():
                    if path[i+1] == nbr:
                        s+=eattr["weight"]
                        boo = True
                        break
    return s/(len(path)-1)

def solve_bubble(graph, ancestor_node, descendant_node):
    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    if paths == []:
        return graph
    return select_best_path(graph, paths, [len(path) for path in paths], [path_average_weight(graph, path) for path in paths])

def simplify_bubbles(graph):
    ancestor_nodes = []
    descendant_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) > 1:
            ancestor_nodes.append(node)
        if len(list(graph.predecessors(node))) > 1:
            descendant_nodes.append(node)
    for ancestor_node in ancestor_nodes:
        for descendant_node in descendant_nodes:
            if ancestor_node in graph.nodes and descendant_node in graph.nodes:
                graph = solve_bubble(graph, ancestor_node, descendant_node)
    return graph

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    starting_nodes = list(graph.nodes)
    for n, nbrs in graph.adj.items():
        for nbr, eattr in nbrs.items():
            if nbr in starting_nodes:
                starting_nodes.remove(nbr)
    return starting_nodes

def get_sink_nodes(graph):
    sink_nodes = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            sink_nodes.append(node)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            paths = list(nx.all_simple_paths(graph, starting_node, ending_node))
            if len(paths) > 0:
                path = paths[0]
                path = path[0] + "".join([path[i][-1] for i in range(1,len(path))])
                contigs.append((path, len(path)))
    return contigs


def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as f:
        for i in range(len(contigs_list)):
            contig = contigs_list[i]
            f.write(">contig_" + str(i) + " " + "len=" + str(contig[1]) + "\n" + fill(contig[0]) + "\n")
        f.close

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    d=build_kmer_dict(args.fastq_file, args.kmer_size)
    #print(d)
    graph = build_graph(d)
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    contigs = get_contigs(graph, starting_nodes, sink_nodes)
    save_contigs(contigs, "debruijn/test.txt")
    graph_2 = nx.DiGraph()
    graph_2.add_weighted_edges_from([(1, 2, 10), (3, 2, 10), (2, 4, 10),
                                     (4, 5, 10), (2, 10,10), (10, 5,10),
                                     (2, 8, 10), (8, 9, 10), (9, 5, 10),
                                     (5, 6, 10), (5, 7, 10)])
    graph_2 = solve_bubble(graph_2, 2, 5)
    print(graph_2.nodes)
    print(graph_2.edges)

if __name__ == '__main__':
    main()
