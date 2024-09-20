# -*- coding: utf-8 -*-
# @Time    : 2024/8/21 19:24
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gff.py

import re
import gzip
import statistics
from collections import defaultdict, Counter

def extract_features_coord(file_path, file_format='gff3', feature='gene', key='ID'):
    """
    Extract the coordinates of the specified feature from the gff file
    :param file_path: str, the path of the gff file
    :param file_format: str, the format of the gff file, default is 'gff3', can be 'gff3' or 'gtf2'
    :param feature: str, the feature to be extracted, default is 'gene'
    :param key: str, the key of the feature, default is 'ID'
    :return: dict, the coordinates of the specified feature
    """
    coord = defaultdict(list)
    patterns = {
        'gff3': re.compile(r'(?P<key>\w+)=(?P<value>[^;]+)'),
        'gtf2': re.compile(r'(?P<key>\w+) "(?P<value>[^"]+)"')
    }
    pattern = patterns[file_format]
    # whether the gff/gtf file is compressed or not
    if file_path.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(file_path, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            line = line.strip().split('\t')
            if line[2] == feature:
                attributes = line[8].split(';')
                attributes = {pattern.match(attr).group('key'): pattern.match(attr).group('value') for attr in attributes}
                if key in attributes:
                    coord[attributes[key]].append([line[0], int(line[3]), int(line[4]), line[6]])
                else:
                    print(attributes)
                    raise ValueError('The key of the feature is not found in the attributes')
    return coord

def sort_features(file_path, file_format='gff3', output_file='sorted.gff'):
    """
    Sort the features in the gff file based on the order of features
    :param file_path: str, the path of the gff file
    :param file_format: str, the format of the gff file, default is 'gff3', can be 'gff3' or 'gtf2'
    :param output_file: str, the path of the output file, default is 'sorted.gff'
    :return: None
    """
    try:
        # Read the input file
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Filter out empty lines and comment lines
        lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]

        # Define the order of features for gff3 format and gtf2 format
        feature_order = {
            'gff3': {
                'gene': 1,
                'mRNA': 2,
                'exon': 3,
                'CDS': 4,
                'five_prime_UTR': 5,
                'three_prime_UTR': 6
            },
            'gtf2': {
                'gene': 1,
                'transcript': 2,
                'exon': 3,
                'CDS': 4,
                'start_codon': 5,
                'stop_codon': 6,
                'UTR': 7
            }
        }.get(file_format)

        # Custom sorting key
        def sort_key(line):
            cols = line.split('\t')
            chromosome = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            feature = cols[2]
            return chromosome, start, feature_order.get(feature, 99), end

        # Sort the lines based on the custom key
        sorted_lines = sorted(lines, key=sort_key)

        # Write the sorted lines to the output file
        with open(output_file, 'w') as f:
            f.write('\n'.join(sorted_lines))
            f.write('\n')

    except FileNotFoundError:
        print("Error: File not found.")
    except IOError:
        print("Error: Input/output error.")
    except Exception as e:
        print("An error occurred:", str(e))

def stat_features(file_path):
    """
    Calculate the statistics of the features in the gff file
    :param file_path: str, the path of the gff file
    :return:
    """
    # Initialize variables
    gene_lengths = []
    mRNA_lengths = []
    mRNA_count_per_gene = Counter()
    CDS_lengths_per_mRNA = Counter()
    CDS_count_per_mRNA = Counter()
    exon_lengths_per_mRNA = Counter()
    exon_count_per_mRNA = Counter()
    exon_end = {}
    intron_lengths_per_mRNA = Counter()
    intron_count_per_mRNA = Counter()
    five_prime_utr_lengths_per_mRNA = []
    three_prime_utr_lengths_per_mRNA = []

    try:
        # Open and read the GFF3 file
        with open(file_path, 'r') as gff:
            for line in gff:
                if not line.startswith('#') and line.strip():
                    # Split the line into a list
                    line = line.strip().split('\t')
                    # Get the feature type
                    feature = line[2]
                    # Get the attributes
                    attributes = line[-1]
                    # Get the feature coordinates
                    start = int(line[3])
                    end = int(line[4])
                    # Calculate the feature length
                    length = end - start + 1
                    # Calculate the gene length
                    if feature == 'gene':
                        gene_lengths.append(length)
                        continue
                    else:
                        # Get the parent ID
                        parent = re.search(r'Parent=([^;]+)', attributes).group(1)

                    # Calculate the mRNA length
                    if feature == 'mRNA':
                        mRNA_lengths.append(length)
                        mRNA_count_per_gene[parent] += 1
                    # Calculate the CDS length
                    elif feature == 'CDS':
                        CDS_lengths_per_mRNA[parent] += length
                        CDS_count_per_mRNA[parent] += 1
                    # Calculate the exon and intron length
                    elif feature == 'exon':
                        if exon_count_per_mRNA[parent] > 0:
                            intron_lengths_per_mRNA[parent] += start - exon_end[parent] - 1
                            intron_count_per_mRNA[parent] += 1
                        exon_lengths_per_mRNA[parent] += length
                        exon_count_per_mRNA[parent] += 1
                        exon_end[parent] = end
                    # Calculate the 5' UTR length
                    elif feature == 'five_prime_UTR':
                        five_prime_utr_lengths_per_mRNA.append(length)
                    # Calculate the 3' UTR length
                    elif feature == 'three_prime_UTR':
                        three_prime_utr_lengths_per_mRNA.append(length)
                    else:
                        continue
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return
    except Exception as e:
        print(f"Error: {e}")
        return

    def print_stats(name, lengths):
        count = len(lengths)
        if count == 0:
            print(f"{name}_count\t{count}")
            return

        max_len = max(lengths)
        min_len = min(lengths)
        mean_len = statistics.mean(lengths)
        median_len = statistics.median(lengths)

        print(f"{name}_count\t{count}")
        print(f"{name}_max\t{max_len}")
        print(f"{name}_min\t{min_len}")
        print(f"{name}_Median/AVG\t{median_len}/{mean_len:.2f}")

    # Calculate the statistics
    print_stats("gene", gene_lengths)
    print_stats("mRNA", mRNA_lengths)
    print_stats("CDS", CDS_lengths_per_mRNA.values())
    print_stats("exon", exon_lengths_per_mRNA.values())
    print_stats("intron", intron_lengths_per_mRNA.values())
    print_stats("five_prime_utr", five_prime_utr_lengths_per_mRNA)
    print_stats("three_prime_utr", three_prime_utr_lengths_per_mRNA)
