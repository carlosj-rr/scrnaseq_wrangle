# genome co-evolutionary networks stuff
from Bio import AlignIO
from Bio import Phylo
from ete3 import PhyloTree
import subprocess
import numpy as np
import pandas as pd

# Function that infers a tree with IQtree for any given alignment:
def infer_tree(ali_obj, ali_filename):
    AlignIO.write(ali_obj, ali_filename, "fasta")
    subprocess.run(["iqtree2","-s",ali_filename, "-redo"])
    tree_filename=ali_filename+".treefile"
    outtree = Phylo.read(tree_filename, "newick")
    # Maybe add a command that removes all the IQtree files, to keep the folders tidy
    return(outtree)
    
# Large function to get something very simple done
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

""" All in one function that can takes an input alignment filename, and a window size, and outputs
the distance table"""
def window_maker(ali_filename,w_size):
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
    rez = []
    for i in trees:
        t1 = i
        col = []
        for j in trees:
            t2 = j
            col.append(PhyloTree.robinson_foulds(t1, t2, unrooted_trees=True)[0])
        rez.append(col)
    axes=[ ''.join(["t",str(i)]) for i in range(len(trees)) ]
    dist_mat=pd.DataFrame(np.array(rez),index=axes,columns=axes)
    return(dist_mat)

def calculate_distances(tree):
    terminals = tree.get_terminals()
    num_terminals = len(terminals)
    distances = np.zeros((num_terminals, num_terminals))
    for i in range(num_terminals):
        for j in range(i+1, num_terminals):
            distance = tree.distance(terminals[i], terminals[j])
            distances[i][j] = distance
            distances[j][i] = distance
    return distances

def scale_tree_to_max_dist(tree):
    distances = calculate_distances(tree)
    scale = distances.max()
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = clade.branch_length/scale
    return(tree)