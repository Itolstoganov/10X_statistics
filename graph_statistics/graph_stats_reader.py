from collections import OrderedDict
import numpy as np
import matplotlib as mpl
mpl.use('Agg')


def read_edges(edges_file):
    edgestats = {'Nonextendable': 0, 'Single': 0, 'Multiple': 0, 'No candidates': 0}
    with open(edges_file, 'r') as f:
        for line in f:
            name, value = line.split(':')[0], int(line.split(':')[1])
            if name in edgestats:
                edgestats[name] = value
    return edgestats


def read_graph_gap_distribution(gaps_file):
    distribution = OrderedDict()
    with open(gaps_file, 'r') as f:
        n = f.readline()
        for line in f:
            split = line.split()
            gap, count = int(split[0]), int(split[1])
            # if gap < 20000:
            distribution[gap] = count
    return distribution


def read_coverage_to_gap_distribution(covgaps_file):
    coverages = []
    gap_distributions = []
    with open(covgaps_file, 'r') as f:
        for line in f:
            coverage = float(line)
            next_line = next(f)
            gap_distribution = [int(word) for word in next_line.split() if int(word) < 10000]
            if len(gap_distribution) > 0 and coverage < 0.075:
                coverages.append(coverage)
                gap_distributions.append(gap_distribution)
    return np.array(coverages), np.array([np.array(distribution) for distribution in gap_distributions])


def read_score_along_path(score_file):
    scores = []
    gaps = []
    with open(score_file, 'r') as f:
        for line in f:
            score_and_gap = line.strip().split()
            # print(score_and_gap)
            scores.append(float(score_and_gap[0]))
            gaps.append(int(score_and_gap[1]))
    return scores, gaps


def read_reliable_barcode_distribution_coverage(input_path):
    gaps = []
    coverages = []
    with open(input_path, 'r') as f:
        for line in f:
            gap_and_coverage = line.strip().split()
            gaps.append(int(gap_and_coverage[0]))
            coverages.append(float(gap_and_coverage[1]))
    return gaps, coverages
