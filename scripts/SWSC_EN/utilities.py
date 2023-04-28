import os
from pathlib2 import Path
from Bio import AlignIO, SeqIO, SeqUtils
from itertools import combinations
import numpy as np
from math import factorial
from Bio.Nexus import Nexus

def check_taxa(matrices):
    '''Checks that nexus instances in a list [(name, instance)...] have
        the same taxa, provides useful error if not and returns None if
        everything matches
        From: http://biopython.org/wiki/Concatenate_nexus
    '''
    first_taxa = matrices[0][1].taxlabels
    for name, matrix in matrices[1:]:
        first_only = [t for t in first_taxa if t not in matrix.taxlabels]
        new_only = [t for t in matrix.taxlabels if t not in first_taxa]

        if first_only:
            missing = ', '.join(first_only)
            msg = '%s taxa %s not in martix %s' % (matrices[0][0], missing, name)
            raise Nexus.NexusError(msg)

        elif new_only:
            missing = ', '.join(new_only)
            msg = '%s taxa %s not in all matrices'  % (name, missing)
            raise Nexus.NexusError(msg)

    return None # will only get here if it hasn't thrown an exception


def concat(mypath, same_taxa):
    ''' Combine multiple nexus data matrices in one partitioned file.
        By default this will work if the same taxa are present in all files
        
        use  same_taxa=False if you are not concerned by this
        From: http://biopython.org/wiki/Concatenate_nexus
        small change: added onlyfiles block to remove hidden files
    '''

    onlyfiles = []
    for item in os.listdir(mypath):
        if not item.startswith('.') and os.path.isfile(os.path.join(mypath, item)):
            onlyfiles.append(item)
    
    nexi = []
    for nex in onlyfiles:
        nex_open = open(nex, 'r')
        nex_save = Nexus.Nexus(nex_open)
        nexi.append((nex, nex_save))
         
    if same_taxa:
        if not check_taxa(nexi):
            return Nexus.combine(nexi)
    else:
        return Nexus.combine(nexi)


def output_conc_nex(mypath, outfilename, same_taxa=False):

    os.chdir(mypath)
    combined = concat(mypath, same_taxa)
    combined.write_nexus_data(filename=open('%s.nex' % (outfilename), 'w'))

    return None


def blocks_pfinder_config(best_window, name, start, stop, uce_aln):

    # sometimes we couldn't split the window so it's all together
    if(best_window[1]-best_window[0] == stop-start):
        whole_UCE = '%s_all = %s-%s;\n' % (name, start+1, stop)
        return (whole_UCE)

    else:
        # left flank
        left_start = start + 1
        left_end = start + best_window[0]

        # core UCE
        core_start = left_end + 1
        core_end = start + best_window[1]

        #right UCE
        right_start = core_end + 1
        right_end = stop

    # do not output any undetermined blocks - if this happens, just output the whole UCE
    if(any_undetermined_blocks(best_window, uce_aln)==True or any_blocks_without_all_sites(best_window, uce_aln)==True):
        whole_UCE = '%s_all = %s-%s;\n' % (name, start+1, stop)
        return (whole_UCE)
    else:
        core_UCE = '%s_core = %s-%s;\n' % (name, core_start, core_end)
        left_UCE = '%s_left = %s-%s;\n' % (name, left_start, left_end)
        right_UCE = '%s_right = %s-%s;\n' % (name, right_start, right_end)

        return (left_UCE + core_UCE + right_UCE)


def any_undetermined_blocks(best_window, uce_aln):
    # Return TRUE if there are any blocks with only undeteremined characters
    # Defined as anything other than ACGT

    left_aln = uce_aln[:, 0 : best_window[0]]
    core_aln = uce_aln[:, best_window[0] : best_window[1]]
    right_aln = uce_aln[:, best_window[1] : uce_aln.get_alignment_length()]

    l_freq = bp_freqs_calc(left_aln)
    c_freq = bp_freqs_calc(core_aln)
    r_freq = bp_freqs_calc(right_aln)

    if(np.isnan(l_freq.max()) or np.isnan(c_freq.max()) or np.isnan(r_freq.max())):
        return(True)
    else:
        return(False)


def any_blocks_without_all_sites(best_window, uce_aln):
    # Return TRUE if there are any blocks with only undeteremined characters
    # Defined as anything other than ACGT

    left_aln = uce_aln[:, 0 : best_window[0]]
    core_aln = uce_aln[:, best_window[0] : best_window[1]]
    right_aln = uce_aln[:, best_window[1] : uce_aln.get_alignment_length()]

    l_counts = count_bases(left_aln)
    c_counts = count_bases(core_aln)
    r_counts = count_bases(right_aln)

    if(np.min(l_counts)==0 or np.min(c_counts)==0 or np.min(c_counts)==0):
        return(True)
    else:
        return(False)


def get_all_windows(aln, minimum_window_size):
    ''' aln: multiple sequence alignment
        minimum_window_size: smallest allowable window 
        
        return a list of all possible tuples [ (start : end) ]
    '''

    length = aln.get_alignment_length()

    keep_windows = []

    if length < 3 * minimum_window_size:
        # some things can't be split
        return ([(0, length)])

    for window in combinations(range(length), 2):
        start = window[0]
        stop = window[1]

        if start < minimum_window_size:
            continue
        if (length - stop) < minimum_window_size:
            continue
        if (stop - start) < minimum_window_size:
            continue

        keep_windows.append(window)

    return (keep_windows)


def output_paths(dataset_path):
    ''' dataset_path: path to a nexus alignment with UCE charsets 
    
        return a folder path with the name of the nexus UCE dataset  
    '''
    
    dataset_name = os.path.basename(dataset_path).rstrip(".nex")

    repository_dir      = Path(dataset_path).parents[1]
    processed_data_dir  = os.path.join(str(repository_dir), "PFinderUCE_output")

    output_path = os.path.join(processed_data_dir, dataset_name)
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    return (output_path)

def p_finder_start_block(dataset_name, branchlengths = 'linked', models = 'GTR+G', model_selection = 'aicc'):
    begin_block = str('## ALIGNMENT FILE ##\n' + 
                      'alignment = %s.phy;\n\n' % (dataset_name) +  


                      '## BRANCHLENGTHS: linked | unlinked ##\n' +
                      'branchlengths = %s;\n\n' % (branchlengths) +

                       '## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai <list> ##\n' +
                       'models = %s;\n\n' % (models) + 

                       '# MODEL SELECCTION: AIC | AICc | BIC #\n' +
                       'model_selection = %s;\n\n' % (model_selection) +

                       '## DATA BLOCKS: see manual for how to define ##\n' +
                       '[data_blocks]\n')

    return (begin_block)


def p_finder_end_block(dataset_name, search = 'rclusterf'):
    ''' dataset_name: name of the dataset
        search: pFinder input arguments
    
        returns str with information about the end block pFinder config block
    '''
    
    end_block = str('\n' +
                    '## SCHEMES, search: all | user | greedy | rcluster | hcluster | kmeans ##\n' +
                    '[schemes]\n' +
                    'search = %s;\n\n' % (search)
                    )

    return (end_block)

def sitewise_base_counts(aln):
    '''
    aln: biopython generic alignment
    
    returns a 4xN array of base counts (A,C,G,T) by site
    '''
    n_sites = aln.get_alignment_length()

    base_counts = np.zeros((4, n_sites))

    for i in range(n_sites):

        site_i = aln[:,i:i+1]
        counts_i = count_bases(site_i)
        base_counts[:,i] = counts_i

    return(base_counts)

def factorial_matrix(counts):
    '''
    input an array of integers
    
    returns an array of the colum-wise products of the factorials of those integers
    '''
    # NB: I used to do this like this, which seemed more sensible and almost certainly quicker,
    #fv = np.vectorize(factorial)
    #factorials = fv(counts)
    #f_product = factorials.prod(axis = 0)
    # however, I kept getting overflow errors that I couldn't figure out, so I gave up.

    f_product = []

    for c in counts.T:
        cf = []
        for i in c:
            cf.append(factorial(i))
        p = cf[0]*cf[1]*cf[2]*cf[3] #np has trouble with big numbers...
        f_product.append(p)

        # sanity check: numpy was occasionally returning negative products of factorials...
        if(p<0):
            raise ValueError('This is bad: you have a negative factorial, which should be impossible here')


    return(f_product)


def count_bases(aln):

    one_str = ""
    for seq in aln:
        one_str += seq

    seq = one_str.upper()

    A = seq.seq.count('A') 
    C = seq.seq.count('C')
    G = seq.seq.count('G')
    T = seq.seq.count('T')

    return(np.array([A, C, G, T]))

def bp_freqs_calc(aln):
    ''' aln: biopython generic alignment
     
    returns a 1-D array of base frequencies
    '''

    base_counts = count_bases(aln)

    sum_count = np.sum(base_counts)

    # fix so we get frequencies of zero when we don't count
    # any A, C, T or G
    if(sum_count == 0): sum_count = 1

    bp_freqs = base_counts/float(sum_count)
    
    return (bp_freqs)


def csv_col_to_plot_matrix(best_window, N):
    
    start = 0
    stop  = N

    wd_left  = np.repeat(-1, [start + best_window[0]])
    wd_core  = np.repeat(0, [best_window[1] - best_window[0]])
    wd_right = np.repeat(1, [stop - best_window[1]])

    concat = np.concatenate([wd_left,wd_core,wd_right])

    return(concat.tolist())

def get_min_var_window(windows, aln_length):
    '''
    input: a set of windows e.g. [(50, 100), (200, 400)]
           the length of the UCE alignment
    
    returns the window with the smallest variance in fragemeng length
    '''

    best_var = np.inf

    for w in windows:
        l1 = w[0]
        l2 = w[1] - w[0]
        l3 = aln_length - w[1]
        var = np.var([l1, l2, l3])

        if var < best_var:
            best_var = var
            best_window = w

    return(best_window)

def alignment_entropy(aln):
    ''' aln: biopython generic alignment
    
    returns an array with values of entropies
    '''

    bp_freqs = bp_freqs_calc(aln)
    entropy = entropy_calc(bp_freqs)
    
    return (entropy)


def entropy_calc(p):
    ''' p: 1D array of base frequencies
    
    returns a estimates of entropies 
    copied from: http://nbviewer.ipython.org/url/atwallab.cshl.edu/teaching/QBbootcamp3.ipynb
    '''
    p = p[p!=0] # modify p to include only those elements that are not equal to 0

    return np.dot(-p,np.log2(p)) # the function returns the entropy result

def sitewise_entropies(aln):
    ''' aln: biopython generic alignment
    
    returns an array with values of entropies per site
    '''

    entropies = []
    for i in range(aln.get_alignment_length()):

        site_i = aln[:,i:i+1]
        ent_i = alignment_entropy(site_i)
        entropies.append(ent_i)

    entropies = np.array(entropies)

    return (entropies)

def all_invariant_sites(sitevar):
    # return TRUE if aln has all invariant sites
    # return FALSE otherwise

    var = np.sum(sitevar)

    if(var==0):
        return True
    else:
        return False

def invariant_sites(aln):
    # output a list of whether or not sites are invariant
    # variable = 1
    # invariant = 0

    entropies = sitewise_entropies(aln)
    entropies[entropies > 0] = 1
    return(entropies)



def split_nexus_by_charsets(aln, charset_name):

    dat = Nexus.Nexus()
    dat.read(aln)
    aln_name = aln.rstrip(".nex")

    for name in charset_name:
        del dat.charsets[name]

    charset_list = []
    for name in tqdm(dat.charsets):
        charset_list.append('%s_%s.data' %  (aln_name, name))

    dat.write_nexus_data_partitions(charpartition = dat.charsets) # tried to use StringIO (doesn't work)

    return(charset_list)

def export_nexus(aln, charset_name):
    nexus_list_names = split_nexus_by_charsets(aln, charset_name)

    nexus_tuples = []
    for name in nexus_list_names:
        nexus_tuples.append((name, Nexus.Nexus(name)))

    concat = Nexus.combine(nexus_tuples)
    concat.write_nexus_data('%s_concat.nex' % (aln.rstrip(".nex")))
