#!/usr/bin/env python3
from operator import itemgetter
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter
from graph_statistics import initial_filter_stats_processor
from graph_statistics import short_edge_dataset_processor
from graph_statistics import long_edge_dataset_processor
from collections import OrderedDict
import os
import numpy as np
import matplotlib as mpl
import pandas
import seaborn as sns
import operator
import logging

mpl.use('Agg')
import matplotlib.pyplot as plt


def make_histogram(array, bin_size):
    binned_array = [value / bin_size * bin_size for value in array]
    counter = Counter(binned_array)
    counter = sorted(counter.items(), key=itemgetter(0))
    x = np.array([value[0] for value in counter])
    y = np.array([value[1] for value in counter])
    return x, y


class ClusterStatsMultiVisualizer(object):
    @classmethod
    def get_name_to_method(cls, name):
        return getattr(cls, "draw_" + name)

    @classmethod
    def get_permitted_names(cls):
        permitted_names = frozenset(["summary_cluster_statistics",
                                     "edge_to_clusters_distribution",
                                     "loop_path_cluster_statistic"])
        return permitted_names

    def __init__(self, output_path, dist_range):
        self.output_path = output_path
        self.dist_range = dist_range

    def filter_dist_list(self, dist_list):
        return [dist for dist in dist_list if self.dist_range[0] <= dist <= self.dist_range[1]]

    @classmethod
    def check_dist_to_stat(self, dist_to_stat):
        prev_name = next(iter(dist_to_stat.values())).get_name()
        for stat in dist_to_stat.values():
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
        ClusterStatsMultiVisualizer.check_dist_to_stats(dist_to_stats)

        number_of_stats = len(next(iter(dist_to_stats.values())))
        logging.info("Number of stats: {}".format(number_of_stats))
        for i in range(number_of_stats):
            dist_to_bam_stat = {int(dist): dist_to_stats[dist][i] for dist in dist_to_stats.keys()}
            assert ClusterStatsMultiVisualizer.check_dist_to_stat(dist_to_bam_stat)
            name = next(iter(dist_to_bam_stat.values())).get_name()
            logging.info("Trying to draw {}".format(name))
            if name in ClusterStatsMultiVisualizer.get_permitted_names():
                logging.info("Drawing {}".format(name))
                distribution = {dist: stat.get_data() for (dist, stat) in dist_to_bam_stat.items()}
                method = ClusterStatsMultiVisualizer.get_name_to_method(name)
                method(self, distribution, output_suffix)

    def get_filtered_list(self, distribution, percent):
        raw_list = sorted(list(Counter(distribution).elements()))
        upper_bound = np.percentile(raw_list, percent)
        return [length for length in raw_list if length < upper_bound]

    def print_stats(self, distribution, filename, percent):
        with open(filename, "w") as fout:
            for distance in sorted(distribution.keys()):
                fout.write('Distance: {}\n'.format(distance))
                list_distr = self.get_filtered_list(distribution[distance], percent)
                fout.write('Mean={}, '.format(np.mean(list_distr)))
                fout.write('Median={}, '.format(np.median(list_distr)))
                fout.write('Standard deviation={}\n'.format(np.std(list_distr)))

    def draw_cluster_span_distribution(self, dist_to_size_dict, output_suffix):
        dist_list = sorted(self.filter_dist_list(list(dist_to_size_dict.keys())))
        percent = 85
        dist_and_len_list = [(dist, length) for dist in dist_list
                             for length in self.get_filtered_list(dist_to_size_dict[dist], percent) if length > 1750]
        labels = ["Distance", "Span"]
        data_frame = pandas.DataFrame.from_records(dist_and_len_list, columns=labels)
        sns.violinplot(x="Distance", y="Span", data=data_frame)
        plt.xlabel("Distance")
        plt.ylabel("Span of the cluster")
        plt.title("Cluster span distribution")
        stat_filename = os.path.join(self.output_path, "cluster_span_distribution_text_" + output_suffix)
        self.print_stats(dist_to_size_dict, stat_filename, 95)
        plt.savefig(os.path.join(self.output_path, "cluster_span_distribution_" + output_suffix))
        plt.clf()

    def draw_cluster_coverage_distribution(self, dist_to_coverage_list, output_suffix):
        dist_list = sorted(self.filter_dist_list(list(dist_to_coverage_list.keys())))
        percent = 90
        dist_and_len_list = [(dist, coverage * 80) for dist in dist_list
                             for coverage in self.get_filtered_list(dist_to_coverage_list[dist], percent)]
        labels = ["Distance", "Coverage"]
        data_frame = pandas.DataFrame.from_records(dist_and_len_list, columns=labels)
        sns.violinplot(x="Distance", y="Coverage", data=data_frame)
        plt.xlabel("Distance")
        plt.ylabel("Coverage of the cluster")
        plt.title("Cluster coverage distribution")
        stat_filename = os.path.join(self.output_path, "cluster_coverage_distribution_text_" + output_suffix)
        dist_to_nucl_coverage = {dist: [el * 100 for el in cov_list] for (dist, cov_list)
                                 in dist_to_coverage_list.items()}
        self.print_stats(dist_to_nucl_coverage, stat_filename, percent)
        plt.savefig(os.path.join(self.output_path, "cluster_coverage_distribution_" + output_suffix))
        plt.clf()

    def draw_barcode_to_cluster_distribution(self, dist_to_size_list, output_suffix):
        dist_list = sorted(self.filter_dist_list(list(dist_to_size_list.keys())))
        percent = 95
        dist_and_len_list = [(dist, length) for dist in dist_list
                             for length in self.get_filtered_list(dist_to_size_list[dist], percent)]
        labels = ["Distance", "Size"]
        data_frame = pandas.DataFrame.from_records(dist_and_len_list, columns=labels)
        sns.violinplot(x="Distance", y="Size", data=data_frame)
        plt.xlabel("Distance")
        plt.ylabel("Number of clusters per container")
        plt.title("Distribution of number of clusters per container")
        stat_filename = os.path.join(self.output_path, "number_of_clusters_distribution_text_" + output_suffix)
        self.print_stats(dist_to_size_list, stat_filename, percent)
        plt.savefig(os.path.join(self.output_path, "number_of_clusters_distribution_" + output_suffix))
        plt.clf()

    def draw_edges_to_clusters_distribution(self, dist_to_distr_dict, output_suffix):
        output_path = os.path.join(self.output_path, "edges_to_cluster_statistics")
        keys = sorted(dist_to_distr_dict.keys())
        threshold = 4
        for key in keys:
            sum_tuple = (0, 0, 0)
            for edge in dist_to_distr_dict[key].keys():
                if int(edge) >= threshold:
                    sum_tuple = tuple(map(operator.add, dist_to_distr_dict[key][edge], sum_tuple))
            dist_to_distr_dict[key][threshold] = sum_tuple

        with open(output_path, "w") as f:
            first_string = "* " + ' '.join([str(el) for el in range(1, threshold)]) + " >" + str(threshold) + "\n"
            f.write(first_string)
            for key in keys:
                f.write(str(key) + " ")
                for elem in dist_to_distr_dict[key].keys():
                    if elem <= threshold:
                        f.write(';'.join([str(x) for x in dist_to_distr_dict[key][elem]]) + " ")
                f.write("\n")

    def draw_edges_to_paths_distribution(self, dist_to_distr_dict, output_suffix):
        output_path = os.path.join(self.output_path, "edges_to_path_statistics")
        keys = sorted(dist_to_distr_dict.keys())
        threshold = 4
        for key in keys:
            sum_tuple = (0, 0)
            for edge in dist_to_distr_dict[key].keys():
                if int(edge) >= threshold:
                    sum_tuple = tuple(map(operator.add, dist_to_distr_dict[key][edge], sum_tuple))
            dist_to_distr_dict[key][threshold] = sum_tuple

        with open(output_path, "w") as f:
            first_string = "* " + ' '.join([str(el) for el in range(1, threshold)]) + " >" + str(threshold) + "\n"
            f.write(first_string)
            for key in keys:
                f.write(str(key) + " ")
                for elem in dist_to_distr_dict[key].keys():
                    if elem <= threshold:
                        f.write(';'.join([str(x) for x in dist_to_distr_dict[key][elem]]) + " ")
                f.write("\n")

    def draw_mean_lengths_for_multiple_dists(self, dist_to_size_dict, output_suffix):
        print(len(dist_to_size_dict))
        dist_to_mean = sorted({dist: np.mean(list(Counter(dist_to_size_dict[dist]).elements()))
                               for dist in dist_to_size_dict.keys()}.items())
        x, y = zip(*dist_to_mean)
        plt.plot(x, y)
        plt.xlabel("Distance")
        plt.ylabel("Mean length of the cluster")
        plt.title("Mean cluster length distribution")
        plt.savefig(os.path.join(self.output_path, "mean_length_distribution" + output_suffix))
        plt.clf()

    def draw_summary_cluster_statistics(self, dist_to_stats, output_suffix):
        sep = "\t"
        names_translator = {"true_pos_rate": "True positive rate",
                            "avg_true_coverage": "Average coverage of true transitions",
                            "average_false_transition_coverage": "Average coverage of false tranisitons",
                            "false_transitions": "False transitions", "true_transitions": "True transitions",
                            "false_transitions_rate": "False transitions rate"}
        output_path = os.path.join(self.output_path, "summary_cluster_statistics" + output_suffix)
        sorted_dists = sorted(dist_to_stats.keys())
        dist_to_record = self._get_dist_to_record(dist_to_stats, sep, names_translator, sorted_dists)
        first_string = "*" + sep + sep.join(sorted(names_translator.values()))
        with open(output_path, "w") as fout:
            fout.write(first_string + "\n")
        self._print_dist_to_record(dist_to_record, sorted_dists, output_path, sep)

    def _get_dist_to_record(self, dist_to_stats, sep, names_translator, sorted_dists):
        dist_to_record = {}
        for dist in sorted_dists:
            stats = dist_to_stats[dist]
            template = stats[0]
            dist_record = [(names_translator[name], [stat[name] for stat in stats])
                           for name in template.keys() if name in names_translator]
            dist_to_record[dist] = sorted(dist_record, key=lambda record: record[0])
        for dist in sorted_dists:
            record = dist_to_record[dist]
        return dist_to_record

    def _print_dist_to_record(self, dist_to_record, sorted_dists, output_path, sep):
        with open(output_path, "w") as fout:
            for dist in sorted_dists:
                record = dist_to_record[dist]
                print("Record: {}".format(record))
                fout.write(sep.join([str(dist)] + [";".join(self._formatize_stats(name_with_stats[1]))
                                                   for name_with_stats in record]) + "\n")

    def _convert_int_to_string(self, num):
        if num == 0:
            return "0"
        current_num = num
        block_size = 3
        module = pow(10, block_size)
        current_list = []
        while current_num > block_size:
            block_value = current_num % module
            block_string = str(block_value)
            for i in range(1, block_size):
                curr_mod = pow(10, i)
                if block_value < curr_mod:
                    block_string = "0" + block_string
            current_list.append(block_string)
            current_num //= module
        if current_num != 0:
            current_list.append(str(current_num))
        head = current_list[len(current_list) - 1]
        curr_idx = 0
        while head[curr_idx] == '0':
            curr_idx += 1
        current_list[len(current_list) - 1] = head[curr_idx:]
        print(current_list[::-1])
        return ",".join(current_list[::-1])

    def _convert_float_to_string(self, float_num):
        return "{0:.2f}".format(float_num)

    def _formatize_string(self, string):
        print("String: {}".format(string))
        if "." in string:
            return self._convert_float_to_string(float(string))
        else:
            return self._convert_int_to_string(int(string))

    def _formatize_stats(self, stats):
        return [self._formatize_string(stat) for stat in stats]

    def draw_loop_path_cluster_statistic(self, dist_to_stats, output_suffix):
        sep = "\t"
        output_path = os.path.join(self.output_path, "loop_path_cluster_statistics" + output_suffix)
        sorted_dists = sorted(dist_to_stats.keys())
        dist_to_record = {}
        print(dist_to_stats)
        for dist, stat in dist_to_stats.items():
            record = stat.items()
            dist_to_record[dist] = sorted(record, key=lambda entry: entry[0])
        print(list(dist_to_record.values())[0])
        sorted_names = sorted([pair[0] for pair in list(dist_to_record.values())[0]])
        with open(output_path, "w") as fout:
            fout.write("*" + sep + sep.join(sorted_names) + "\n")
            for dist in sorted_dists:
                record = dist_to_record[dist]
                fout.write(str(dist) + sep + sep.join([self._formatize_string(name_with_stat[1])
                                                       for name_with_stat in record]) + "\n")


class GraphStatsVisualizer(object):
    def __init__(self, output_path):
        self.output_path_ = output_path

    @classmethod
    def get_permitted_names(cls):
        # score_distribution_info
        permitted_names = frozenset(["initial_filter_stats", "threshold_distance",
                                     "short_edge_dataset", "long_edge_dataset"])
        return permitted_names

    @classmethod
    def get_method_from_name(cls, name):
        return getattr(cls, "draw_" + name)

    def draw_stats(self, stats):
        permitted_names = GraphStatsVisualizer.get_permitted_names()

        for stat in stats.get_leaves():
            stat_name = stat.get_name()
            print("Stat name: " + stat_name)
            if stat_name in permitted_names:
                method = GraphStatsVisualizer.get_method_from_name(name=stat_name)
                method(self, stat.get_data(), stat.get_name())

    def draw_initial_filter_stats(self, stats, output_suffix):
        output_path = os.path.join(self.output_path_, output_suffix)
        extractor = initial_filter_stats_processor.InitialFilterExtractor(stats=stats)
        processor = initial_filter_stats_processor.InitialFilterProcessor(extractor=extractor)
        next_map = processor.get_next_edge_map()
        score_error_statistics = processor.get_score_error_statistics()
        print(score_error_statistics)

        for top_threshold in range(5, 6):
            print("Top threshold: {}".format(top_threshold))
            top_candidate_statistics = processor.get_top_candidates_statistics(top_threshold=top_threshold)
            print(top_candidate_statistics)

        self.draw_next_distribution(processor.get_next_score_distribution(), self.output_path_)

    def draw_long_edge_dataset(self, long_edge_data, output_suffix):
        output_path = os.path.join(self.output_path_, output_suffix)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        long_edge_data = long_edge_dataset_processor.prepare_dataset(data=long_edge_data)
        # long_edge_dataset_processor.draw_2d_plots(data=long_edge_data, output_path=output_path)
        long_edge_dataset_processor.get_stats(data=long_edge_data)
        long_edge_dataset_processor.draw_cont_index_histogram(data=long_edge_data, output_path=output_path)
        long_edge_dataset_processor.draw_violin_plot(data=long_edge_data, output_path=output_path)

    def draw_short_edge_dataset(self, short_edge_data, output_suffix):
        output_path = os.path.join(self.output_path_, output_suffix)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        # short_data = short_edge_dataset_processor.clean_and_add_columns(short_data=short_edge_data)
        # short_edge_dataset_processor.draw_violin_plots(short_data=short_data, output_path=output_path)
        # short_edge_dataset_processor.draw_score_to_params(short_data=short_data, output_path=output_path)

        # short_edge_dataset_processor.test_criterias(data=short_data)
        # short_edge_dataset_processor.draw_3d_plot(short_data=short_data, output_path=output_path)
        # short_edge_dataset_processor.launch_svm(data=short_data, output_path=output_path)

    def draw_threshold_distance(self, stats, output_suffix):
        output_path = os.path.join(self.output_path_, output_suffix);
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        distances, thresholds, data = stats
        sns.heatmap(data=data)
        plt.savefig(os.path.join(output_path, "threshold_failed_heatmap"))

    def draw_score_distribution_info(self, stats, output_suffix):
        output_path = os.path.join(self.output_path_, output_suffix)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        for name, distr_pair in stats.items():
            close_scores, random_scores = distr_pair
            close_distr = np.array(close_scores)
            random_distr = np.array(random_scores)
            sns.distplot(close_distr, kde=False, hist_kws={"label": "Score between close edges", "cumulative": True},
                         color="green", bins=len(close_scores), norm_hist=True)
            sns.distplot(random_distr, kde=False, hist_kws={"label": "Score between distant edges", "cumulative": -1},
                         color="red", norm_hist=True)
            percentile = 0.2
            sns.plt.xlim(0, np.percentile(close_scores, int(percentile * 100)))
            sns.plt.ylim(0, percentile)
            plt.xlabel("Score")
            plt.ylabel("Probability of smaller/larger score for close/distant transitions")
            plt.title("Cumulative score distributions")
            plt.legend()
            plt.savefig(os.path.join(output_path, "score_distribution_info_" + name))
            plt.clf()

    def draw_next_distribution(self, distribution, output_base):
        distr = np.array(distribution)
        print([x for x in distribution if x < 5])
        print(len([x for x in distribution if x < 5]))
        print(type(distr))
        sns.distplot(distr)
        plt.savefig(os.path.join(output_base, "next_score_distribution"))
        plt.clf()
