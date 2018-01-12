from collections import OrderedDict
from collections import namedtuple
import numpy as np
import pandas as pd
import matplotlib as mpl
import os
import sys
import logging

mpl.use('Agg')


class AbstractGraphStats(object):
    def __init__(self, name):
        self.name_ = name
        self.sons = []

    def add_son(self, son):
        self.sons.append(son)

    def get_name(self):
        return self.name_

    def get_leaves(self):
        leaves = []
        self._leaves_search(leaves=leaves, node=self)
        return leaves

    def __len__(self):
        return len(self.get_leaves())

    def _leaves_search(self, leaves, node):
        for son in node.sons:
            logging.debug("Son: {}".format(son.get_name()))
            logging.debug("Grandsons: {}".format(son.sons))
            if isinstance(son, LeafGraphStats):
                logging.debug("Appending {}".format(son.name_))
                logging.debug("{} sons".format(len(son.sons)))
                leaves.append(son)
            else:
                self._leaves_search(leaves=leaves, node=son)


class LeafGraphStats(AbstractGraphStats):
    def __init__(self, name, data):
        super().__init__(name)
        self.data_ = data
        self.name_ = name

    def get_data(self):
        return self.data_

    def get_name(self):
        return self.name_


class ClusterSpanDistribution(LeafGraphStats):
    def __init__(self, name):
        super(ClusterSpanDistribution, self).__init__(name=name, data={})

    def load_data(self, path):
        size_to_number = {}
        with open(path, "r") as fin:
            for line in fin:
                size, number = [int(val) for val in line.split(" ")]
                size_to_number[size] = number
        self.data_ = size_to_number


class ClusterCoverageDistribution(LeafGraphStats):
    def __init__(self, name):
        super(ClusterCoverageDistribution, self).__init__(name=name, data=[])

    def load_data(self, path):
        coverages = []
        with open(path, "r") as fin:
            for line in fin:
                coverages.append(float(line))
        self.data_ = coverages


class BarcodeToClustersDistribution(LeafGraphStats):
    def __init__(self, name):
        super(BarcodeToClustersDistribution, self).__init__(name=name, data=[])

    def load_data(self, path):
        cluster_sizes = []
        with open(path, "r") as fin:
            for line in fin:
                cluster_sizes.append(int(line))
        self.data_ = cluster_sizes


class EdgesToClustersDistribution(LeafGraphStats):
    def __init__(self, name):
        super(EdgesToClustersDistribution, self).__init__(name=name, data={})

    def load_data(self, path):
        edge_number_to_paths = {}
        with open(path, "r") as fin:
            for line in fin:
                edges, total, single, correct = [int(val) for val in line.split(" ")]
                edge_number_to_paths[edges] = (total, single, correct)
        self.data_ = edge_number_to_paths


class EdgesToPathsDistribution(LeafGraphStats):
    def __init__(self, name):
        super(EdgesToPathsDistribution, self).__init__(name=name, data={})

    def load_data(self, path):
        edge_number_to_paths = {}
        with open(path, "r") as fin:
            for line in fin:
                edges, single, correct = [int(val) for val in line.split(" ")]
                edge_number_to_paths[edges] = (single, correct)
        self.data_ = edge_number_to_paths


class SummaryClusterStatistics(LeafGraphStats):
    def __init__(self, name):
        super(SummaryClusterStatistics, self).__init__(name=name, data={})

    def load_data(self, path):
        path_cluster_stats = {}
        nonpath_cluster_stats = {}
        with open(path, "r") as fin:
            names = fin.readline().split()
            path_cluster_stats = self._load_stats(line=fin.readline(), names=names)
            nonpath_cluster_stats = self._load_stats(line=fin.readline(), names=names)
        self.data_ = (path_cluster_stats, nonpath_cluster_stats)

    def _load_stats(self, line, names):
        stats_split = line.split()
        if len(names) != len(stats_split):
            raise ValueError("Number of stats does not correspond to the number of values")
        stats = {names[i]: stats_split[i] for i in range(len(names))}
        return stats


class LoopPathClusterStatistic(LeafGraphStats):
    def __init__(self, name):
        super(LoopPathClusterStatistic, self).__init__(name=name, data={})

    def load_data(self, path):
        sep = "\t"
        with open(path, "r") as fin:
            names = fin.readline().strip().split(sep)
            stats = self._load_stats(line=fin.readline(), names=names, sep=sep)
        self.data_ = stats

    def _load_stats(self, line, names, sep):
        stats_split = line.strip().split(sep)
        if len(names) != len(stats_split):
            raise ValueError("Number of stats does not correspond to the number of values")
        stats = {names[i]: stats_split[i] for i in range(len(names))}
        return stats


FirstEdgeStats = namedtuple("FirstEdgeStats", ["max_score", "cov", "barcodes", "next_id", "rc_id"])
EdgeInfo = namedtuple("EdgeInfo", ["first_edge_stats", "candidates"])
CandidateInfo = namedtuple("CandidateInfo", ["shared", "cov", "barcodes", "score", "distance"])


class InitialFilterStats(LeafGraphStats):
    def __init__(self, name):
        super(InitialFilterStats, self).__init__(name=name, data={})

    def load_data(self, path):
        with open(path, "r") as fin:
            number_of_edges = int(fin.readline())
            edge_info_dict = {}
            for _ in range(number_of_edges):
                edge_line = fin.readline()
                edge_id, edge_value = self.read_edge_info(edge_line, fin)
                edge_info_dict[edge_id] = edge_value
            self.data_ = edge_info_dict

    def read_edge_info(self, edge_line, file_stream):
        edge_stats_split = edge_line.split()
        # print("Edge split: {}".format(edge_stats_split))
        edge_id, first_edge_stats, number_of_candidates = self.read_edge_stats_split(edge_stats_split)
        candidates = {}
        for _ in range(number_of_candidates):
            cand_stats_split = file_stream.readline().split()
            cand_id, cand_value = self.read_candidate_stats_split(cand_stats_split=cand_stats_split)
            candidates[cand_id] = cand_value
        edge_info = EdgeInfo(first_edge_stats=first_edge_stats, candidates=candidates)
        file_stream.readline()
        return edge_id, edge_info

    def read_edge_stats_split(self, edge_stats_split):
        logging.debug("Edge stats: {}".format(edge_stats_split))
        edge_id = int(edge_stats_split[0])
        number_of_candidates = int(edge_stats_split[1])
        max_score = float(edge_stats_split[2])
        cov = float(edge_stats_split[3])
        barcodes = int(edge_stats_split[4])
        next_id = int(edge_stats_split[5])
        rc_id = int(edge_stats_split[6])
        edge_info = FirstEdgeStats(max_score=max_score, cov=cov, barcodes=barcodes, next_id=next_id, rc_id=rc_id)
        return edge_id, edge_info, number_of_candidates

    def read_candidate_stats_split(self, cand_stats_split):
        # print(cand_stats_split)
        logging.debug(cand_stats_split)
        id = int(cand_stats_split[0])
        shared = int(cand_stats_split[1])
        cov = float(cand_stats_split[2])
        barcodes = int(cand_stats_split[3])
        score = float(cand_stats_split[4])
        distance = int(cand_stats_split[5])
        cand_info = CandidateInfo(shared=shared, cov=cov, barcodes=barcodes, score=score,
                                  distance=distance)
        return id, cand_info


class ScoreDistributionInfo(LeafGraphStats):
    def __init__(self, name):
        super(ScoreDistributionInfo, self).__init__(name=name, data={})

    def load_data(self, path):
        with open(path, "r") as fin:
            number_of_score_functions = int(fin.readline())
            for _ in range(number_of_score_functions):
                name = fin.readline().strip()
                close_scores_string = fin.readline()
                random_scores_string = fin.readline()
                close_scores = [float(word) for word in close_scores_string.split()]
                random_scores = [float(word) for word in random_scores_string.split()]
                self.data_[name] = (close_scores, random_scores)


class ThresholdDistance(LeafGraphStats):
    def __init__(self, name):
        super(ThresholdDistance, self).__init__(name=name, data={})

    def load_data(self, path):
        with open(path, "r") as fin:
            sep = "\t"
            distances_len = int(fin.readline())
            thresholds_len = int(fin.readline())
            distances, thresholds = [], []
            thresholds_string = fin.readline()
            print(thresholds_string)
            thresholds = [int(word) for word in thresholds_string.strip().split(sep)[1:]]
            print(thresholds)
            assert (len(thresholds) == thresholds_len)
            failed_array = np.zeros(shape=(distances_len, thresholds_len))
            for i in range(distances_len):
                row_string = fin.readline()
                row_string_list = [int(word) for word in row_string.strip().split(sep)]
                distance = row_string_list[0]
                distances.append(distance)
                assert len(row_string_list) == len(thresholds) + 1
                for j in range(1, len(row_string_list)):
                    failed_array[i][j - 1] = row_string_list[j]
            self.data_ = (distances, thresholds, failed_array)


class ShortEdgeDataset(LeafGraphStats):
    def __init__(self, name):
        super(ShortEdgeDataset, self).__init__(name=name, data={})

    def load_data(self, path):
        self.data_ = pd.read_table(path, sep=',')


NameStructureNode = namedtuple("StructureInternalNode", ["name", "sons"])


def get_cluster_structure():
    return NameStructureNode(name="cluster_statistics",
                             sons=["barcode_to_clusters_distribution",
                                   "cluster_coverage_distribution",
                                   "cluster_span_distribution",
                                   "edges_to_clusters_distribution",
                                   "edges_to_paths_distribution",
                                   "summary_cluster_statistics"])


def get_contracted_cluster_structure():
    return NameStructureNode(name="contracted_cluster_statistics",
                             sons=["loop_path_cluster_statistic"])


def get_standard_structure():
    cluster_statistics = get_cluster_structure()
    contracted_cluster_statistics = get_contracted_cluster_structure()
    initial_filter_stats = NameStructureNode(name="initial_filter_stats",
                                             sons=["initial_filter_stats"])
    scaffolder_stats = NameStructureNode(name="scaffolder_statistics",
                                         sons=["score_distribution_info",
                                               "threshold_distance",
                                               "short_edge_dataset"])
    graph_stats = NameStructureNode(name="barcode_stats",
                                    sons=[cluster_statistics,
                                          contracted_cluster_statistics,
                                          initial_filter_stats,
                                          scaffolder_stats])
    return graph_stats


class GraphStatsReader:
    def __init__(self, base_path, structure=get_standard_structure()):
        self.path_ = base_path
        self.structure_ = structure

    def load_data(self):
        path = os.path.join(self.path_, self.structure_.name)
        return self.traverse_structure(self.structure_, path, self.structure_.name)

    def traverse_structure(self, structure, relative_path, name):
        node = AbstractGraphStats(name)
        for item in structure.sons:
            if isinstance(item, NameStructureNode):
                son = self.traverse_structure(item, os.path.join(relative_path, item.name), item.name)
                node.add_son(son)
            else:
                if not isinstance(item, str):
                    raise TypeError("Unknown entry type in the structure")
                stat_class = get_class_from_stat_name(item)
                stat = stat_class(item)
                stat_path = os.path.join(relative_path, item)
                logging.info("Loading from " + stat_path)
                if os.path.exists(stat_path):
                    stat.load_data(os.path.join(relative_path, item))
                    node.add_son(stat)
        return node

    def get_stats_path(self):
        return self.path_


class ClusterMultiStatsReader:
    def __init__(self, stats_paths, structure=get_cluster_structure()):
        self.readers = []
        for path in stats_paths:
            stats_reader = GraphStatsReader(base_path=path, structure=structure)
            self.readers.append(stats_reader)

    def extract_dist(self, path):
        path_string = os.path.normpath(path)
        dist_name = path_string.split(os.sep)[-1]
        return dist_name.split("_")[1]

    def read_dist_to_stats(self):
        dist_to_stats = {}
        logging.info("Distance to stats:")
        for reader in self.readers:
            logging.info(reader.get_stats_path())
            dist = self.extract_dist(reader.get_stats_path())
            logging.info(dist)
            data = reader.load_data()
            print(data)
            print(data.get_name())
            print(data.get_leaves())
            dist_to_stats[dist] = data.get_leaves()
        return dist_to_stats


def get_class_from_stat_name(name):
    camel_case_name = ''.join(word.capitalize() for word in name.split("_"))
    return getattr(sys.modules[__name__], camel_case_name)
