#!/usr/bin/env python
from operator import itemgetter
from collections import Counter
from collections import OrderedDict
import scipy.stats as stats
import os
import numpy as np
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter


# plt.rcParams['agg.path.chunksize'] = 20000


def make_histogram(array, bin_size):
    binned_array = [value / bin_size * bin_size for value in array]
    counter = Counter(binned_array)
    counter = sorted(counter.items(), key=itemgetter(0))
    x = np.array([value[0] for value in counter])
    y = np.array([value[1] for value in counter])
    return x, y


def plot_fragment_length_distribution(outpath, lengths):
    lengths = lengths[lengths > 500]
    lengths = lengths[lengths < 90000]
    weights = np.ones_like(lengths)/float(len(lengths))
    plt.hist(lengths, bins=50, weights=weights)
    plt.xlabel("Length of the fragment")
    plt.ylabel("Probability of given fragment length")
    plt.title("Fragment length distribution")
    plt.savefig(os.path.join(outpath, "length_distribution"))
    plt.clf()


def plot_pool_size_distribution(outpath, fragment_to_reads):
    left_size = 100
    plt.hist(fragment_to_reads, bins=100)
    # plt.xlim(left_size, max(fragment_to_reads))
    plt.ylabel("Number of fragments")
    plt.xlabel("Number of barcoded reads in a fragment (normalized by length)")
    plt.yscale('log')
    plt.title("Fragment coverage distribution")
    plt.savefig(os.path.join(outpath, "read_distribution"))
    plt.clf()


def plot_gap_distribution(outpath, gaps):
    gaps = [gap for gap in gaps if gap > 1000]
    x, y = make_histogram(gaps, 100)
    y = np.cumsum(y[::-1])[::-1]
    y = y / float(y[0])
    plt.xlabel("Distance between consecutive reads inside a fragment")
    plt.ylabel("Probability of gap being larger than given value")
    plt.title("Internal gap length cumulative distribution (reference estimation)")
    plt.plot(x, y)
    plt.xlim(-20, np.percentile(gaps, 99))
    plt.savefig(os.path.join(outpath, "gap_distribution"))
    plt.clf()


def plot_internal_coverage_distribution(outpath, length_to_count):
    read_length = 140
    coverages = [tup[1] * read_length / float(tup[0]) for tup in length_to_count if tup[0] > 5000]
    # print(coverages) 
    print(min(coverages))
    print(len(coverages))
    print(max(coverages))
    weights = np.ones_like(coverages)/float(len(coverages))
    plt.xlabel("Internal coverage of a fragment")
    plt.ylabel("Probability of given internal coverage")
    plt.title("Internal coverage distribution")
    plt.hist(coverages, 50, range=(0, np.percentile(coverages, 95)), weights=weights)
    plt.savefig(os.path.join(outpath, "internal_coverage"))
    plt.clf()

def draw_pie_edges(plot_path, edgestats):
    # mpl.rcParams['font.size'] = 60
    # mpl.rcParams['figure.figsize'] = 30, 22
    # labels = 'Single extension', 'Multiple extensions', 'No extensions'
    sizes = [edgestats['Single'], edgestats['Multiple'],
             edgestats['Nonextendable'], edgestats['No candidates']]
    colors = ['green', 'yellow', 'red', 'gray']
    explode = (0, 0.1, 0.1, 0)  # explode 1st slice
    pieWedgesCollection = plt.pie(sizes, explode=explode, colors=colors,
                                  shadow=True, startangle=90, radius=1)
    plt.axis('equal')
    out_path = os.path.join(plot_path, "edges_pie")
    plt.savefig(out_path)
    plt.clf()

def plot_score_along_path(outpath, scores):
    plt.hist(scores, bins=50)
    plt.xlabel("Fraction of shared barcodes between adjacent edges in genome path")
    plt.ylabel("Number of edges")
    plt.title("Score distribution")
    plt.savefig(os.path.join(outpath, "genome_path_score"))
    plt.clf()


def plot_score_with_gaps(outpath, scores, gaps):
    scores_with_gaps = sorted(zip(gaps, scores), key=lambda tup: tup[0])
    x, y = zip(*scores_with_gaps)
    x, y = x[:-10], y[:-10]
    bins = 30
    bin_means, bin_edges, binnumber = stats.binned_statistic(x, y, statistic='mean', bins=bins)
    plt.xlabel("Gap between adjacent edges in genome path")
    plt.ylabel("Mean score")
    plt.title("Relation between gap and score")
    plt.bar(bin_edges[:-1], bin_means, width=max(x) / bins, color='green')
    plt.savefig(os.path.join(outpath, "scores_with_gaps"))
    plt.clf()


def plot_barcode_to_len(outpath, barcode_to_len):
    plt.hist(barcode_to_len, bins=100)
    plt.xlabel("Number of fragments in a barcode")
    plt.ylabel("Number of barcodes")
    plt.title("Distribution of number of long fragments in a single barcode")
    plt.xlim(0, 20)
    plt.savefig(os.path.join(outpath, "fragments_in_barcode"))
    plt.clf()


def plot_graph_gap_distribution(distribution, outpath):
    gap_threshold = 500
    new_gaps = OrderedDict((gap, distribution[gap]) for gap in distribution.keys() if gap >= gap_threshold)
    x, y = new_gaps.keys(), new_gaps.values()
    y = np.cumsum(y[::-1])[::-1]
    y = y / float(y[0])
    plt.xlabel("Distance between consecutive reads inside a fragment")
    plt.ylabel("Probability of gap being larger than given value")
    plt.title("Internal gap length cumulative distribution (long edge estimation)")
    plt.plot(x, y)
    plt.xlim(-20, 20000)
    plt.savefig(os.path.join(outpath, "gap_distribution"))
    plt.clf()



def plot_cov_gap_distribution(x, y, outpath):
    print(y.shape)
    y_medians = np.array([np.median(distribution) for distribution in y])
    bins = 50
    bin_means, bin_edges, binnumber = stats.binned_statistic(x, y_medians, statistic = 'mean', bins=bins)
    widths = [bin_edges[i + 1] - bin_edges[i] for i in range(len(bin_means))]
    print(widths)
    plt.xlabel("Number of reads in the fragment normalized by length of the covered part")
    plt.ylabel("Median of internal gap distribution")
    plt.title("Median gap with respect to coverage")
    plt.bar(bin_edges[:-1], bin_means, width=widths, color='green')
    plt.savefig(os.path.join(outpath, "cov_gap_distribution"))


def plot_gap_to_reliable_coverage(gaps, coverages, outpath):
    plt.xlabel("Length of reliable cloud")
    plt.ylabel("Covered fraction of long edges")
    plt.title("Estimated fraction of genome covered by reliable clouds")
    plt.plot(gaps, coverages)
    plt.savefig(os.path.join(outpath, "reliable_distribution"))
    plt.clf()