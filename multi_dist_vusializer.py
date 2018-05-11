#!/usr/bin/env python3
import pandas
import os
import shutil
import argparse
from collections import Counter

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from bam_statistics import bam_stats_io


def createparser():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', help="Path to directories with statistics", required=True)
    parser.add_argument('-o', '--output', help="Path with plots (in the input directory by default)")
    return parser


class MultiStats(object):
    def __init__(self, dist_to_stats):
        self.dist_to_stats = dist_to_stats


class MultiStatsReader(object):
    def __init__(self, stats_paths):
        self.readers = []
        for path in stats_paths:
            stats_reader = bam_stats_io.BamStatsReader(stats_path=path)
            self.readers.append(stats_reader)

    def extract_dist(self, path):
        path_string = os.path.normpath(path)
        dist_name = path_string.split(os.sep)[-2]
        return dist_name.split("_")[1]

    def read_dist_to_stats(self):
        dist_to_stats = {}
        for reader in self.readers:
            dist = self.extract_dist(reader.stats_path)
            dist_to_stats[dist] = reader.read_bam_stats()
        return dist_to_stats


class MultiDistVisualizer(object):
    @classmethod
    def get_name_to_method(cls):
        return {
            "fragment_length_distribution": getattr(cls, "draw_fragment_length_distributions"),
            "fragments_per_container": getattr(cls, "draw_fragments_per_container"),
            "internal_coverage_distribution": getattr(cls, "draw_initial_coverage_distributions")
        }

    def __init__(self, output_path, dist_range):
        self.output_path = output_path
        self.dist_range = dist_range

    def filter_dist_list(self, dist_list):
        return [dist for dist in dist_list if self.dist_range[0] <= dist <= self.dist_range[1]]

    @classmethod
    def check_dist_to_bam_stat(self, dist_to_bam_stat):
        prev_name = next(iter(dist_to_bam_stat.values())).get_name()
        for stat in dist_to_bam_stat.values():
            name = stat.get_name()
            if not name == prev_name:
                return False
        return True

    @classmethod
    def check_dist_to_stats(cls, dist_to_stats):
        number_of_stats = len(next(iter(dist_to_stats.values())))
        for value in dist_to_stats.values():
            if len(value) != number_of_stats:
                raise IndexError("Statistics are incomplete")

    def draw_multi_stats(self, dist_to_stats, output_suffix):
        MultiDistVisualizer.check_dist_to_stats(dist_to_stats)

        number_of_stats = len(next(iter(dist_to_stats.values())))
        for i in range(number_of_stats):
            dist_to_bam_stat = {int(dist): dist_to_stats[dist][i] for dist in dist_to_stats.keys()}
            assert MultiDistVisualizer.check_dist_to_bam_stat(dist_to_bam_stat)
            name = next(iter(dist_to_bam_stat.values())).get_name()
            name_to_method = MultiDistVisualizer.get_name_to_method()
            if name in name_to_method:
                distribution = {dist: stat.get_data() for (dist, stat) in dist_to_bam_stat.items()}
                name_to_method[name](self, distribution, output_suffix)
                if name == "fragment_length_distribution":
                    self.draw_mean_lengths_for_multiple_dists(distribution, output_suffix)

    def get_filtered_list(self, dictionary):
        raw_list = list(Counter(dictionary).elements())
        upper_bound = np.percentile(raw_list, 95)
        return [length for length in raw_list if length < upper_bound]

    def draw_fragment_length_distributions(self, dist_to_length_distributions, output_suffix):
        dist_list = sorted(self.filter_dist_list(list(dist_to_length_distributions.keys())))
        dist_and_len_list = [(dist, length) for dist in dist_list
                             for length in self.get_filtered_list(dist_to_length_distributions[dist])]
        labels = ["Distance", "Length"]
        len_data_frame = pandas.DataFrame.from_records(dist_and_len_list, columns=labels)
        sns.violinplot(x="Distance", y="Length", data=len_data_frame)
        plt.xlabel("Distance")
        plt.ylabel("Length of the cluster")
        plt.title("Cluster length distribution")
        plt.savefig(os.path.join(self.output_path, "multi_length_distribution" + output_suffix))
        plt.clf()

    def draw_fragments_per_container(self, dist_to_fragments_per_container, output_suffix):
        dist_list = sorted(self.filter_dist_list(list(dist_to_fragments_per_container.keys())))
        dist_and_cov_list = [(dist, length) for dist in dist_list
                             for length in self.get_filtered_list(dist_to_fragments_per_container[dist])]
        labels = ["Distance", "Clusters"]
        cov_data_frame = pandas.DataFrame.from_records(dist_and_cov_list, columns=labels)
        sns.boxplot(x="Distance", y="Clusters", data=cov_data_frame)
        plt.xlabel("Distance")
        plt.ylabel("Number of clusters per container")
        plt.title("Number of clusters per container distribution")
        plt.savefig(os.path.join(self.output_path, "multi_fragment_distribution" + output_suffix))
        plt.clf()

    def draw_initial_coverage_distributions(self, dist_to_initial_coverage_distributions, output_suffix):
        dist_list = sorted(self.filter_dist_list(list(dist_to_initial_coverage_distributions.keys())))
        dist_and_cov_list = [(dist, length) for dist in dist_list for length in
                             self.get_filtered_list(dist_to_initial_coverage_distributions[dist])]
        labels = ["Distance", "Coverage"]
        cov_data_frame = pandas.DataFrame.from_records(dist_and_cov_list, columns=labels)
        sns.violinplot(x="Distance", y="Coverage", data=cov_data_frame)
        plt.xlabel("Distance")
        plt.ylabel("Coverage of the cluster")
        plt.title("Cluster coverage distribution")
        plt.savefig(os.path.join(self.output_path, "multi_coverage_distribution" + output_suffix))
        plt.clf()

    def draw_mean_lengths_for_multiple_dists(self, dist_to_length_distributions, output_suffix):
        dist_list = sorted(self.filter_dist_list(list(dist_to_length_distributions.keys())))
        dist_to_mean = sorted({dist: np.mean(list(Counter(dist_to_length_distributions[dist]).elements()))
                               for dist in dist_list}.items())
        x, y = zip(*dist_to_mean)
        plt.plot(x, y)
        plt.xlabel("Distance")
        plt.ylabel("Mean length of the cluster")
        plt.title("Mean cluster length distribution")
        plt.savefig(os.path.join(self.output_path, "mean_length_distribution" + output_suffix))
        plt.clf()


def select_stats(base_input_path, prefix, stats_prefix):
    subpaths = [path for path in os.listdir(base_input_path) if path.startswith(prefix)]
    print("Subpaths: ", subpaths)
    stats_paths = [os.path.join(base_input_path, subpath, stats_prefix) for subpath in subpaths]
    return stats_paths


if __name__ == "__main__":
    parser = createparser()
    args = parser.parse_args()
    base_input_path = args.input
    output_path = args.output
    print(output_path)

    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)

    print("Input", base_input_path)
    stats_paths = select_stats(base_input_path=base_input_path, prefix="dist", stats_prefix="stats")
    multireader = MultiStatsReader(stats_paths)
    dist_to_stats = multireader.read_dist_to_stats()
    print(len(dist_to_stats))

    very_large_dist_range = (100000, 500000)
    large_dist_range = (10000, 70000)
    small_dist_range = (500, 5000)
    print("Drawing plots for very large distances")
    multivisualizer = MultiDistVisualizer(output_path=output_path, dist_range=very_large_dist_range)
    multivisualizer.draw_multi_stats(dist_to_stats, output_suffix="_very_large")
    print("Drawing plots for large distances")
    multivisualizer = MultiDistVisualizer(output_path=output_path, dist_range=large_dist_range)
    multivisualizer.draw_multi_stats(dist_to_stats, output_suffix="_large")
    print("Drawing plots for small distances")
    multivisualizer = MultiDistVisualizer(output_path=output_path, dist_range=small_dist_range)
    multivisualizer.draw_multi_stats(dist_to_stats, output_suffix="_small")
