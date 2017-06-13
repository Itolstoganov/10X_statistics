#!/usr/bin/env python3
import argparse
import logging
import math
import os
import shutil
from collections import Counter
from abc import ABCMeta, abstractmethod

import matplotlib as mpl
import numpy as np

import bam_stats_io

mpl.use('Agg')
import matplotlib.pyplot as plt


def createparser():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', help="Path to directory with statistics", required=True)
    parser.add_argument('-o', '--output', help="Path with plots (next to input directory by default)")
    return parser


class BaseBamStatsVisualizer(metaclass=ABCMeta):
    @classmethod
    @abstractmethod
    def get_supported_names(cls):
        raise NotImplementedError()

    @classmethod
    def get_name_to_method_dict(cls):
        name_to_method_dict = {}
        for name in cls.get_supported_names():
            method_name = "draw_" + name
            try:
                method = getattr(cls, method_name)
                name_to_method_dict[name] = method
            except AttributeError:
                raise NotImplementedError(
                    "Class `{}` does not implement `{}`".format(cls.__class__.__name__, method_name))
        return name_to_method_dict

    def draw_statistics(self, statistics):
        name_to_method_dict = self.get_name_to_method_dict()
        for stat in statistics:
            self.draw_stat(bam_stat=stat, name_to_method_dict=name_to_method_dict)

    def draw_stat(self, bam_stat, name_to_method_dict):
        if bam_stat.get_name() in name_to_method_dict:
            method = name_to_method_dict[bam_stat.get_name()]
            method(self, bam_stat.get_data())
        # else:
        #     logging.warning('Stat {} is not drawable'.format(bam_stat.get_name()))


class OldBamStatsVisualizer(BaseBamStatsVisualizer):
    @classmethod
    def get_supported_names(cls):
        return ["fragment_length_distribution", "internal_coverage_distribution",
                "fragments_per_container", "reads_per_container", "read_reference_coverage",
                "fragment_reference_coverage"]

    def __init__(self, color, output_path):
        self.color = color
        self.output_path = output_path

        self.output_statistics_filename = os.path.join(output_path, "distribution_statistics")
        open(self.output_statistics_filename, "w")

    def print_distribution_params(self, name, distribution):
        with open(self.output_statistics_filename, "a") as fout:
            fout.write(name + ": \n")
            fout.write("Mean: ")
            fout.write(str(sum(distribution) / float(len(distribution))))
            fout.write('\n')
            fout.write("Median: ")
            fout.write(str(distribution[len(distribution) // 2]))
            fout.write('\n')
            fout.write("Variance: ")
            fout.write(str(np.var(distribution)))
            fout.write("\n\n")

    def draw_fragment_length_distribution(self, lengths_dict):
        length_counter = Counter(lengths_dict)
        lengths = list(length_counter.elements())
        lengths = [length for length in lengths if length > 1000]
        weights = np.ones_like(lengths) / float(len(lengths))
        self.print_distribution_params(name="Fragment length distribution", distribution=lengths)
        plt.hist(lengths, bins=50,
                 range=(0, np.percentile(lengths, 95)),
                 weights=weights, color=self.color)
        plt.xlabel("Length of the cluster")
        plt.ylabel("Probability of given cluster length")
        plt.title("Cluster length distribution")
        plt.savefig(os.path.join(self.output_path, "length_distribution"))
        plt.clf()

    def draw_internal_coverage_distribution(self, coverages_dict):
        coverage_counter = Counter(coverages_dict)
        coverages = list(coverage_counter.elements())
        weights = np.ones_like(coverages) / float(len(coverages))
        self.print_distribution_params(name="Internal coverage distribution", distribution=coverages)
        plt.xlabel("Internal coverage of a cluster")
        plt.ylabel("Probability of given internal coverage")
        plt.title("Distribution of the coverage of clusters")
        plt.hist(coverages, bins=50,
                 range=(0, np.percentile(coverages, 97)),
                 weights=weights, color=self.color)
        plt.savefig(os.path.join(self.output_path, "internal_coverage"))
        plt.clf()

    def draw_fragments_per_container(self, fragments_per_barcode):
        fragments_counter = Counter(fragments_per_barcode)
        fragments_list = list(fragments_counter.elements())
        self.print_distribution_params(name="Fragment per container", distribution=fragments_list)
        weights = np.ones_like(fragments_list) / float(len(fragments_list))
        plt.xlabel("Number of clusters")
        plt.ylabel("Fraction of containers with given number of clusters")
        plt.title("Distribution of the number of clusters per container")
        x_right_lim = math.ceil(np.percentile(fragments_list, 95))
        xrange = (1, x_right_lim)

        plt.xlim(xrange)
        plt.xticks(np.arange(1, x_right_lim + 1, 2))
        plt.hist(fragments_list, weights=weights, bins=np.arange(x_right_lim + 2), range=xrange,
                 color=self.color)
        plt.savefig(os.path.join(self.output_path, "fragments_per_barcode"))
        plt.clf()

    def draw_reads_per_container(self, reads_per_container):
        reads_counter = Counter(reads_per_container)
        reads_list = list(reads_counter.elements())
        self.print_distribution_params(name="Fragment per container", distribution=reads_list)
        weights = np.ones_like(reads_list) / float(len(reads_list))
        plt.xlabel("Number of reads")
        plt.ylabel("Fraction of containers with given number of reads")
        plt.title("Distribution of the number of reads per container")
        plt.xlim(0, 51)
        plt.hist(reads_list, bins=50, range=(0, 51), weights=weights, color=self.color)
        plt.savefig(os.path.join(self.output_path, "reads_per_container"))
        plt.clf()

    def draw_read_reference_coverage(self, reads_per_reference):
        with open(self.output_statistics_filename, "a") as fout:
            fout.write("Read reference coverages\n")
            for (key, value) in reads_per_reference.items():
                fout.write('{}\t{}\n'.format(key, value))
            fout.write('\n')

    def draw_fragment_reference_coverage(self, fragments_per_reference):
        with open(self.output_statistics_filename, "a") as fout:
            fout.write("Fragment reference coverages\n")
            for (key, value) in fragments_per_reference.items():
                fout.write('{}\t{}\n'.format(key, value))
            fout.write('\n')


class NewBamStatsVisualizer(BaseBamStatsVisualizer):
    @classmethod
    def get_supported_names(cls):
        return ["fragment_length_distribution", "internal_coverage_distribution",
                "fragments_per_container", "reads_per_container"]

    def __init__(self, output_path):
        """
        @param output_path: Where pictures should be stored
        """
        self.output_path = output_path

    def draw_fragment_length_distribution(self, length_distribution):
        """
        Draws distribution of fragment lengths.
        @type  length_distribution: dict
        @param length_distribution: Map <length of fragment>: <number of such fragments>
        """
        pass

    def draw_internal_coverage_distribution(self, coverages_dict):
        """
        Draws distribution of internal coverages.
        @type  length_distribution: dict
        @param length_distribution: Map <internal coverage of a fragment>: <number of such fragments>
        """
        pass

    def draw_fragments_per_container(self, fragments_per_barcode):
        """
        Draws distribution of number of fragments per container(barcode).
        @type  length_distribution: dict
        @param length_distribution: Map <number of fragments>: <number of such containers>
        """
        pass

    def draw_reads_per_container(self, reads_per_container):
        """
        Draws distribution of number of reads per container(barcode).
        @type  length_distribution: dict
        @param length_distribution: Map <number of reads>: <number of such containers>
        """
        pass


if __name__ == "__main__":
    parser = createparser()
    args = parser.parse_args()
    output_path = args.output
    if not args.output:
        output_path = os.path.join(args.input, "..", "pictures")

    print("Output path: ", output_path)
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)

    bam_stats_reader = bam_stats_io.BamStatsReader(args.input)
    bam_stats = bam_stats_reader.read_bam_stats()

    # bam_stats_visualizer = OldBamStatsVisualizer(color="green", output_path=output_path)
    bam_stats_visualizer = NewBamStatsVisualizer(output_path=output_path)
    bam_stats_visualizer.draw_statistics(bam_stats)
