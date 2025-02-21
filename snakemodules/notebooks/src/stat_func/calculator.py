from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from scipy.stats import entropy
import os
import time
import collections
import math
import subprocess as sp
import pandas as pd
import numpy as np

def double(lst):
    return [i*2 for i in lst]

def cal_entropy(sequence):
    """Calculate entropy for a sequence of numbers
 
    Args:
        sequence : contig sequence
 
    Returns:
        float: entropy
 
    """
    d = dict()
    counted = 0

    for base_i in range(len(sequence)-1):
        dimer = sequence[base_i:base_i+2]
        if (dimer[0]!='N' and dimer[1]!='N'):
            # if dimer is new to dict
            if not dimer in d:
                d[dimer] = 1
                counted+=1
            # otherwise
            else:
                d[dimer] += 1
                counted+=1

    # calculate the entropy for dinucleotide counts
    entropy = 0
    for key, value in d.items():
        if value == 0:
            continue
        p = float(value) / counted
        entropy -= p * math.log(p) / math.log(2)

    return entropy / 4

def estimate_shannon_entropy(dna_sequence):
    """Calculate entropy for each sequence
 
    Returns:
        float: entropy value
 
    """
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
    # define distribution
    dist = [x/sum(bases.values()) for x in bases.values()]
 
    # use scipy to calculate entropy
    entropy_value = entropy(dist, base=2)
 
    return entropy_value, bases

def entropy_dimer(sequence):
        d = dict()
        counted = 0

        for base_i in range(len(sequence)-1):
            dimer = sequence[base_i:base_i+2]
            if (dimer[0]!='N' and dimer[1]!='N'):
                # if dimer is new to dict
                if not dimer in d:
                    d[dimer] = 1
                    counted+=1
                # otherwise
                else:
                    d[dimer] += 1
                    counted+=1

        # calculate the entropy for dinucleotide counts
        entropy = 0
        for key, value in d.items():
            if value == 0:
                continue
            p = float(value) / counted
            entropy -= p * math.log(p) / math.log(2)

        return entropy / 4

def cal_N50(p):
    """Calculate N50 for a sequence of numbers
 
    Args:
        p : filepath
        len_of_sequences (list): List of sequence lengths
 
    Returns:
        float: median, N50 value.
        float: contigNum, number of contigs in the sample
 
    """
    seq_lengths=[]

    for seq_record in SeqIO.parse(p,'fasta'):  
        length=len(seq_record.seq)
        seq_lengths.append(length)
        contigNum = len(seq_lengths)


    all_len=sorted(seq_lengths, reverse=True)
    csum=np.cumsum(all_len)

    # n = int(sum(seq_lengths))
    n2=int(sum(seq_lengths)/2)

    # get index for cumsum >= N/2
    csumn2=min(csum[csum >= n2])
    ind=np.where(csum == csumn2)

    N50 = all_len[int(ind[0])]

 
    return N50, contigNum

def calculate_N50(len_of_sequences):
    """Calculate N50 for a sequence of numbers
 
    Args:
        len_of_sequences (list): List of sequence lengths
 
    Returns:
        float: N50 value.
 
    """
    tmp = []
    for tmp_num in set(len_of_sequences):
            tmp += [tmp_num] * len_of_sequences.count(tmp_num) * tmp_num
    tmp.sort()
 
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
 
    return median

def cal_gctable(p):

    """Calculate total GC percentage for assembled fasta file
 
    Args:
        p : the full filepath
 
    Returns:
        float: total GC% for the whole assambled sample
 
    """
    A=[]
    C=[]
    G=[]
    T=[]

    
    for seq_record in SeqIO.parse(p,'fasta'):
        seq = repr(seq_record.seq)
        A += [seq.count("A")]
        C += [seq.count("C")]
        G += [seq.count("G")]
        T += [seq.count("T")]

    A_T=sum(A)+sum(T)
    G_C=sum(G)+sum(C)

    sample_gc = round(G_C/(A_T + G_C)*100)

    return sample_gc


def process_sequences(sample_seq,sample_names):
    """
    Processes a list of sequence groups and computes statistics for each group.

    Args:
        sample_seq (list): A list where each element is a list of sequence records.
                           Each record should have attributes 'id' and 'seq'.

    Returns:
        list: A list of dictionaries containing computed statistics for each group of sequences.
              Each dictionary has the keys 'Seq_ID', 'Seq_LEN', 'GC_Content', and 'Entropy'.
    """
    sample_list = []
    sample_dict = {}

    for sequences in sample_seq:
        seq_ids = []
        seq_lengths = []
        gc_contents = []
        entropy_values = []

        for record in sequences:
            seq_id = record.id
            sequence = record.seq
            length = len(sequence)
            gc_content = gc_fraction(sequence)
            gc_content = round(gc_content, 2)
            entropy_value = entropy_dimer(sequence)

            seq_ids.append(seq_id)
            seq_lengths.append(length)
            gc_contents.append(gc_content)
            entropy_values.append(entropy_value)

        d = {
            'Seq_ID': seq_ids,
            'Seq_LEN': seq_lengths,
            'GC_Content': gc_contents,
            'Entropy': entropy_values
        }
        sample_list.append(d)

    for n,m in zip(sample_names, range(len(sample_list))):
        sample_dict[n] = pd.DataFrame(sample_list[m])

    return sample_list, sample_dict



# Function to replace outliers in each column with a set maximum value
def replace_column_outliers_with_max(df, column):
    # for column in df.columns:
        # Calculate the 25th and 75th percentiles
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1

    # Define outlier as outside of Q1 - 1.5 * IQR and Q3 + 1.5 * IQR
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    print(f'upper_bound {upper_bound}\n IQR {IQR}')
    # Replace outliers with the set maximum value
    outliers = ((df[column] > upper_bound))
    # outliers = df[column] > upper_bound
    print(df.loc[outliers, column])
    # df.loc[outliers] = upper_bound
    # print(df[column][outliers])
    df.loc[outliers, column] = upper_bound
    
    return df

def remove_outliers_replace_with_bounds(df):
    """
    Removes outliers from the 1st and 4th quartiles and replaces them with the lower and upper bounds.
    
    Parameters:
    df: DataFrame where the first column is the index and the next columns are the contig coverages for samples.
    
    Returns:
    df: DataFrame with outliers replaced by the respective bounds.
    """
    
    # Apply the process to each column starting from the second (ignoring index column)
    for col in df.columns:
        # Calculate Q1 (25th percentile) and Q3 (75th percentile)
        Q1 = df[col].quantile(0.25)
        Q3 = df[col].quantile(0.75)
        
        # Calculate IQR
        IQR = Q3 - Q1
        
        # Define bounds for outliers
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        
        # Replace values below the lower bound with the lower bound
        df[col] = df[col].apply(lambda x: lower_bound if x < lower_bound else x)
        
        # Replace values above the upper bound with the upper bound
        df[col] = df[col].apply(lambda x: upper_bound if x > upper_bound else x)
    
    return df


#quick function for getting an average of a list of numbers as a string    
def quickstats(x):
    x_array = np.asarray(x,dtype=np.float64)
    x_avg = np.average(x_array)
    x_avg = str(np.around(x_avg, decimals = 1))
    return x_avg


def calculate_coverage(bamfiles,out_folder):
    """
    Calculates coverage for a list of BAM files using samtools depth.

    For each BAM file provided, this function checks if a coverage file already exists.
    If not, it calls samtools depth to calculate the coverage and saves the output.
    It also prints the progress of the calculations.

    Args:
        bamfiles (list): List of file paths to BAM files.

    Returns:
        None: Coverage files are generated and saved in the current directory.
    """

    # Ensure the output folder exists
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Begin sample counter
    x = 1

    # Iterate across all sample names
    sample_names = []

    for fp in bamfiles:
        folder_path = os.path.dirname(fp)
        folder = os.path.basename(folder_path)
        sample_names.append(folder)
        out_name = os.path.join(out_folder, f"{folder}.coverage.txt")

        print("Sample {0}, {1}: Now calculating coverage for {2}".format(x, out_name, fp))

        # System call to samtools depth function, write to output designated above
        cov_string = "samtools depth {0} > {1}".format(fp, out_name)
        proc_cov = sp.call(cov_string, shell=True)

        x += 1

        print(time.asctime(time.localtime(time.time())))


def compile_coverage_results(output_file, coverage_dir='.', quickstats_func=None):
    """
    Compiles average coverage per contig from individual coverage files into a master output file.

    Args:
        output_file (str): The name of the master output file to be created.
        coverage_dir (str): The directory where the coverage files are located.
        quickstats_func (callable): A function to calculate the average coverage from a list of coverages.
                                    If None, the built-in average calculation will be used.

    Returns:
        None: The function writes the output to the specified file.
    """

    # Open the output file in write mode to override if it exists
    with open(output_file, 'w') as fh_out:
        # For all files in the coverage directory
        for fl in os.listdir(coverage_dir):
            # Find coverage output files
            if fl.endswith('.coverage.txt'):
                coverage_file = os.path.join(coverage_dir, fl)
                print("Starting calculations at", time.asctime(time.localtime(time.time())))
                print("Working to calculate averages for all contigs in {}:".format(fl))

                # Create a set to hold unique contig names
                temp_set = set()

                # Open the coverage file
                with open(coverage_file, 'r') as fh_temp:
                    # Read all lines
                    lines = fh_temp.readlines()
                    # Extract sample name from file name
                    sample_name = fl.split('.')[0]

                    # Extract contig names
                    for line in lines:
                        line = line.strip().split('\t')
                        temp_set.add(line[0])

                # Sort contig names
                temp_list = sorted(list(temp_set))
                number = len(temp_list)
                print('\t', "There are {} contigs to sift through...".format(number))

                # Process each contig
                for contig in temp_list:
                    print('\t', '\t', "- Calculating average for {}".format(contig))
                    concovlist = []

                    # Collect coverage values for the contig
                    for line in lines:
                        line = line.strip().split('\t')
                        if line[0] == contig:
                            concovlist.append(float(line[2]))

                    # Perform averaging for this contig
                    contig_avg = quickstats_func(concovlist)
                    

                    # Write sample name and average coverage for each contig
                    fh_out.write(f"{sample_name}\t{contig}\t{contig_avg}\n")
                    
                print('-----------------------------------------------------------------', '\n')
