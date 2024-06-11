# genome co-evolutionary networks stuff
from Bio import AlignIO
from Bio import Phylo
from ete3 import PhyloTree
import subprocess
import numpy as np
import pandas as pd
from io import StringIO

"""Takes in an alignment object and its output filename, runs IQTREE, and imports the resulting tree
as a newick string"""
def infer_tree(ali_obj, ali_filename):
    AlignIO.write(ali_obj, ali_filename, "fasta")
    subprocess.run(["iqtree2","-s",ali_filename, "-redo"])
    tree_filename=ali_filename+".treefile"
    outtree = treefile_to_str(tree_filename)
    # Maybe add a command that removes all the IQtree files, to keep the folders tidy
    return(outtree)
    
"""Makes a range of indices that split a given range into windows of the same size. If the window size
doesn't split the range evenly, it creates a longer final window if the modulo is less than 1/2 of the
window size, and creates a small addtional window otherwise."""
def rangemaker(begin=0, end=None, step=None):
    # Make range of start values and rename the endpoint
    starts = range(begin, end, step)
    tot = end
    # Exit early if the step divides the total cleanly
    if tot % step ==0:
        stops = np.array(starts)+step-1
        return(list(starts),stops)
    # Check whether the remainder of dividing the total by the step is less than 1/2 the step, if so, remove the last value of the starts list
    if tot % step < step/2:
        starts = starts[0:len(starts)-1]
    stops = []
    # For each start value, let's create a stop value
    for i in range(len(starts)):
        # The stop value is the step times the window number, minus one (so it doesn't overlap with the next start)
        stp_val = step*(i+1)-1
        # If the stop value being calculated is the last one, it should be the end point, regardless of the step size
        if i == len(starts)-1:
            stp_val = end
        stops.append(stp_val)
    return(list(starts),stops)

""" Takes an input alignment filename, and a window size, and outputs the tree for each window as
inferred by IQTREE (as a newick string)"""
def window_treemaker(ali_filename,w_size):
    full_ali = AlignIO.read(ali_filename, "fasta")
    ali_len = len(full_ali[0].seq)
    starts, stops = rangemaker(0, ali_len, w_size)
    ali_chunks=[]
    chunks_fnames = []
    for i in range(len(starts)):
        curr_start = starts[i]
        curr_stop = stops[i]
        curr_chunk = full_ali[:,curr_start:curr_stop]
        curr_file = "chunk_"+str(i)+".fas"
        ali_chunks.append(curr_chunk)
        chunks_fnames.append(curr_file)
    trees=list(map(infer_tree,ali_chunks,chunks_fnames))
    return(trees)

""" Calculate all the tip to tip distances in a tree"""
def calculate_tip2tip_distances(tree):
    terminals = tree.get_terminals()
    num_terminals = len(terminals)
    distances = np.zeros((num_terminals, num_terminals))
    for i in range(num_terminals):
        for j in range(i+1, num_terminals):
            distance = tree.distance(terminals[i], terminals[j])
            distances[i][j] = distance
            distances[j][i] = distance
    return distances

""""Scales all branch lengths to the longest tip-tip distance"""
def scale_tree_to_max_dist(tree):
    distances = calculate_tip2tip_distances(tree)
    scale = distances.max()
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = clade.branch_length/scale
    return(tree)

"""Reads a tree as a newick string and outputs a Biopython Phylo tree"""
def strtoBio(nwck_str):
    biopyth_t=Phylo.read(StringIO(nwck_str), "newick")
    return biopyth_t

"""Reads a tree as a newick string and outputs an ETE3 PhyloTree object"""
def strtoEte3(nwk_str):
    outtree = PhyloTree(nwk_str)
    return outtree

"""Reads a the filename of a treefile, and outputs the tree as a string in newick format"""
def treefile_to_str(filename):
    with open(filename, 'r') as file:
        nwk_str = file.read().replace("\n","")
    return nwk_str

"""Extract all branch lengths from a Biopython Phylo tree and outputs a list of blengths"""
def get_blengths(Biotree):
    out = []
    for clade in Biotree.find_clades():
        if clade.branch_length:
            out.append(clade.branch_length)
    return out

"""Takes in two trees as newick strings, and outpus the RF-distance between them"""
def RF_distance(t1_str, t2_str):
    t1, t2 = strtoEte3(t1_str), strtoEte3(t2_str)
    return t1.robinson_foulds(t2, unrooted_trees=True)[0]

"""Takes in two trees as newick strings (with the same topology) and outputs the absolute value of
the Pearson correlation coefficient of a regression of their corresponding branch lengths."""
def blength_corr(t1_str, t2_str):
    t1 = scale_tree_to_max_dist(strtoBio(t1_str))
    t2 = scale_tree_to_max_dist(strtoBio(t2_str))
    t1.ladderize()
    t2.ladderize()
    t1_blen = get_blengths(t1)
    t2_blen = get_blengths(t2)
    return abs(np.corrcoef(t1_blen, t2_blen)[0,1])

"""Calculates the distance between two trees - if the RF distance is 0, it gives the correlation
calculated with the 'blength_corr()' function"""
def dist_getter(t1_str, t2_str):
    rf = RF_distance(t1_str, t2_str)
    if rf > 0:
        output = rf
    else:
        output = blength_corr(t1_str, t2_str)
    return output

##################### main #####################
if False:
    def main(infile, k):
        trees_list = window_treemaker(infile, k)
        out = []
        for i in trees_list:
            for j in trees_list:
                out.append(dist_getter(i,j))
            out_arr = np.array(out).reshape(len(trees_list), len(trees_list))
        return out_arr