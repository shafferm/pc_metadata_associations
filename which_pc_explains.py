import pandas as pd
from scipy.stats import spearmanr
import argparse
import numpy as np

#TODO: implement use with categorical metadata labels
#TODO: for continuous data allow pearson correlation/use correl_method arg
#TODO: for continuous data sort correls

def read_pc(pc_fp, min_percent_explained):
    f = open(pc_fp)
    f = f.readlines()
    f = [i.strip().split('\t') for i in f]

    eig_index = -1
    eigvals = None
    percent_explained_line = -1
    percent_explained = None
    pc_start = -1
    pc_end = -1
    pcs = None
    for index, line in enumerate(f):
        if line == ['']:
            continue
        elif line[0] == "Eigvals":
            eig_index = index + 1
        elif index == eig_index:
            eigvals = line
        elif line[0] == "Proportion explained":
            percent_explained_line = index + 1
        elif index == percent_explained_line:
            percent_explained = [float(i) for i in line]
        elif line[0] == "Site":
            pc_start = index + 1
            pc_end = pc_start + int(line[1])
            pcs = pd.DataFrame(dtype=float)
        elif index >= pc_start and index < pc_end:
            pcs[line[0]] = [float(i) for i in line[1:]]
        else:
            continue
    pcs = pcs.transpose()

    # get rid of zeroes
    not_zeroes = pcs.columns[pcs.sum() != 0]
    pcs = pcs[not_zeroes]
    percent_explained = percent_explained[:len(not_zeroes)]

    # filter by percent explained
    last_percent_explained = None
    for index, value in enumerate(percent_explained):
        if value < .001:
            last_percent_explained = index-1
            break
    percent_explained = percent_explained[:last_percent_explained]
    pcs = pcs[pcs.columns[:last_percent_explained]]
    pcs.columns = ["PC"+str(i+1) for i in pcs.columns]

    return pcs


def bh_adjust(pvalues):
    """benjamini-hochberg p-value adjustment stolen from
    http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
    """
    pvalues = np.array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in xrange(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def main(args):
    pcs = read_pc(args.input_fp, args.min_explained)
    meta = pd.read_table(args.mapping_fp, index_col=0)

    if args.is_continuous:
        correls = pd.DataFrame(index=["R value", "p"])
        for pc_axis in pcs:
            corr = spearmanr(meta[args.metadata_label], pcs[pc_axis])
            correls[pc_axis] = corr
        correls = correls.transpose()
        correls['p_adj'] = bh_adjust(correls['p'])
        correls.to_csv(args.output_fp, sep='\t')

    else:
        raise NotImplementedError()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_fp", help="table of principle components only", required=True)
    parser.add_argument("-m", "--mapping_fp", help="mapping file location", required=True)
    parser.add_argument("-o", "--output_fp", help="output file location", default="pc_correlations.txt")
    parser.add_argument("-c", "--metadata_label", help="minimum number of samples present in", required=True)
    parser.add_argument("--correl_method", help="correlation method", default="spearman")
    parser.add_argument("--is_continuous", help="minimum p-value to determine edges", default=False, action="store_true")
    parser.add_argument("--min_explained", help="minimum percent explained to be included in analysis", default=0., type=float)

    args = parser.parse_args()

    main(args)
