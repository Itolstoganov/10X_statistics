#!/usr/bin/env python3
import argparse
import logging
import math
import os
import sys
import shutil
from collections import Counter
from abc import ABCMeta, abstractmethod
import pandas
import matplotlib as mpl
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from bam_statistics import bam_stats_io


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
        if bam_stat.get_name() in name_to_method_dict and len(bam_stat.get_data()) != 0:
            print(bam_stat.get_name())
            print(len(bam_stat.get_data()))
            method = name_to_method_dict[bam_stat.get_name()]
            method(self, bam_stat.get_data())
        # else:
        #     logging.warning('Stat {} is not drawable'.format(bam_stat.get_name()))


class OldBamStatsVisualizer(BaseBamStatsVisualizer):
    @classmethod
    def get_supported_names(cls):
        return ["fragment_length_distribution", "internal_coverage_distribution",
                "fragments_per_container", "reads_per_container", "read_reference_coverage",
                "fragment_reference_coverage", "covered_references_per_barcode", 
                "barcode_to_reads", "score_distributions", "gap_distribution"]

    def __init__(self, color, output_path):
        self.color = color
        self.output_path = output_path

        self.output_statistics_filename = os.path.join(output_path, "distribution_statistics")
        open(self.output_statistics_filename, "w")

    def print_distribution_params(self, name, distribution):
        percentile = np.percentile(distribution, 95)
        distribution = [element for element in distribution if element < percentile]
        with open(self.output_statistics_filename, "a") as fout:
            fout.write(name + ": \n")
            fout.write("Mean=")
            fout.write(str(np.mean(distribution)))
            fout.write(", ")
            fout.write("Median=")
            fout.write(str(np.median(distribution)))
            fout.write(", ")
            fout.write("Standard deviation=")
            fout.write(str(np.std(distribution)))
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
        # plt.plot(lengths)
        # plt.xlim(0, np.percentile(lengths, 95))
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
        xrange = (0, x_right_lim)
        plt.xlim(xrange)
        plt.xticks(np.arange(1, x_right_lim + 1, 2))
        print(xrange)
        plt.hist(fragments_list, weights=weights, bins=np.arange(x_right_lim + 2), range=xrange,
                 color=self.color)
        plt.savefig(os.path.join(self.output_path, "fragments_per_barcode"))
        plt.clf()

    def draw_reads_per_container(self, reads_per_container):
        reads_counter = Counter(reads_per_container)
        reads_list = list(reads_counter.elements())
        self.print_distribution_params(name="Reads   per container", distribution=reads_list)
        weights = np.ones_like(reads_list) / float(len(reads_list))
        plt.xlabel("Number of reads")
        plt.ylabel("Fraction of containers with given number of reads")
        plt.title("Distribution of the number of reads per container")
        # plt.xlim(0, 51)
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

    def draw_covered_references_per_barcode(self, references_per_barcode):
        refs = len(references_per_barcode.keys())
        counter = Counter(references_per_barcode)
        ref_list = list(counter.elements())
        self.print_distribution_params(name="Covered genomes per barcode", distribution=ref_list)
        weights = np.ones_like(ref_list) / float(len(ref_list))
        plt.xlabel("Number of covered genomes")
        plt.ylabel("Fraction of barcodes covering given number of genomes")
        plt.title("Distribution of covered genomes per barcode")
        plt.hist(ref_list, bins=refs, weights=weights, color=self.color)
        plt.xticks(np.arange(1, refs))
        plt.savefig(os.path.join(self.output_path, "references_per_container"))
        plt.clf()

    def draw_gap_distribution(self, gap_distribution):
        left_bound = 1000
        right_bound = 200000
        gaps = sorted([gap for gap in list(Counter(gap_distribution).elements())])
        print(max(gaps))
        gaps = sorted([gap for gap in gaps if right_bound > gap > left_bound])
        print(max(gaps))
        if len(gaps) > 0:
            self.print_distribution_params(name="Gap distribution", distribution=gaps)
        weights = np.ones_like(gaps) / float(len(gaps))
        plt.xlabel("Length of the gap")
        plt.ylabel("Fraction of gaps with given length")
        plt.title("Distribution of the gap length")
        plt.hist(gaps, bins=100, weights=weights, color=self.color)
        plt.savefig(os.path.join(self.output_path, "gap_distribution"))
        plt.clf()

    def draw_barcode_to_reads(self, barcode_to_reads):
        output_name = os.path.join(self.output_path, "barcode_to_reads")
        with open(output_name, "w") as f:
            for key, value in barcode_to_reads.items():
                f.write(key[:16] + " " + str(value) + "\n")


    def draw_score_distributions(self, score_distributions):
        close_distribution = score_distributions[0]
        distant_distribution = score_distributions[1]
        cov_bin_size = 2
        actual_close_scores = [entry for entry in close_distribution if entry[1] > 0.001]
        actual_distant_scores = [entry for entry in distant_distribution if entry[1] > 0.001]

        data = []
        data += [(entry[0], int(entry[1] / cov_bin_size) * cov_bin_size, "next") for entry in actual_close_scores]
        data += [(entry[0], int(entry[1] / cov_bin_size) * cov_bin_size, "close") for entry in actual_distant_scores]
        coverages = [entry[1] for entry in data]
        min_bin = np.percentile(coverages, 5)
        max_bin = np.percentile(coverages, 95)
        data = [entry for entry in data if min_bin <= entry[1] <= max_bin]

        labels = ["Score", "Coverage", "State"]
        cov_data_frame = pandas.DataFrame.from_records(data, columns=labels)
        sns.violinplot(x="Coverage", y="Score", hue="State", data=cov_data_frame)
        plt.title("Conditional distribution of score given coverage")
        plt.savefig(os.path.join(self.output_path, "score_distributions"))
        plt.clf()


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
        @type  coverages_dict: dict
        @param coverages_dict: Map <internal coverage of a fragment>: <number of such fragments>
        """
        pass

    def draw_fragments_per_container(self, fragments_per_barcode):
        """
        Draws distribution of number of fragments per container(barcode).
        @type  fragments_per_barcode: dict
        @param fragments_per_barcode: Map <number of fragments>: <number of such containers>
        """
        pass

    def draw_reads_per_container(self, reads_per_container):
        """
        Draws distribution of number of reads per container(barcode).
        @type  reads_per_container: dict
        @param reads_per_container: Map <number of reads>: <number of such containers>
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

    bam_stats_visualizer = OldBamStatsVisualizer(color="green", output_path=output_path)
    # bam_stats_visualizer = NewBamStatsVisualizer(output_path=output_path)
    bam_stats_visualizer.draw_statistics(bam_stats)
