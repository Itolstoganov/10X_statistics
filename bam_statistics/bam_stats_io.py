import logging
import os
import pickle

import bam_stats


class BamStatsPrinter:
    def __init__(self, output_base):
        self.output_base = output_base

    def print_bam_stats(self, stats):
        filename = stats.get_name()
        output_path = os.path.join(self.output_base, filename)
        with open(output_path, 'wb') as fout:
            pickle.dump(stats.get_data(), fout)


class BamStatsReader:
    def __init__(self, stats_path):
        self.stats_path = stats_path

    def read_bam_stats(self):
        name_to_class = bam_stats.name_to_class
        files = os.listdir(self.stats_path)
        stats = []
        for filename in files:
            logging.info("Reading ", os.path.join(self.stats_path, filename))
            name = filename
            with open(os.path.join(self.stats_path, filename), "rb") as fin:
                data = pickle.load(fin)
            stats.append(name_to_class[name](data))
        return stats
