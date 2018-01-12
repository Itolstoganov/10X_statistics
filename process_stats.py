#! /usr/bin/env python3
import argparse
import logging
import os
import shutil

import pysam

from bam_statistics import bam_stats_visualizer, bam_stats_extractor, bam_stats_io, primary_statistics
from graph_statistics import graph_stats_reader, graph_stats_visualizer


def createparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile', help="Barcoded BAM File")
    parser.add_argument('-t', '--tag', help="Barcode tag name", default="BX:Z")
    parser.add_argument('--test', help="Test mode: use only first 1000000 reads from BAM", action='store_true')
    parser.add_argument('-d', '--distance', help="Distance at which clusters are splitted", type=int, default=10000)
    parser.add_argument('-o', '--output', help="Path with plots")
    parser.add_argument('-g', '--graph', help="Path to graph_statistics statistics")
    return parser


def prepare_output_dirs(output_name):
    if not os.path.exists(output_name):
        os.mkdir(output_name)
    log_outname = os.path.join(output_name, "log")
    saves_outname = os.path.join(output_name, "stats")
    pictures_outname = os.path.join(output_name, "pictures")

    os.mkdir(saves_outname)
    os.mkdir(pictures_outname)

    open(log_outname, "w")
    FORMAT = '[%(asctime)s] [%(levelname)8s] (%(filename)25s:%(lineno)4s) --- %(message)s'

    logging.basicConfig(format=FORMAT,
                        level=logging.INFO, filename=log_outname)
    return saves_outname, pictures_outname


def process_bamfile(bamname, tag, test, distance, output_name):
    saves_name, pictures_name = prepare_output_dirs(output_name=output_name)
    bamfile = pysam.Samfile(bamname, "rb")

    logging.info('Tag: {}'.format(tag))
    logging.info('Distance: {}'.format(distance))
    logging.info('Test mode: {}'.format(test))
    params = primary_statistics.ProcessorParameters(gap_threshold=distance,
                                                    tag=tag,
                                                    test_mode=test)
    processor = primary_statistics.BamStatisticsProcessor(bamfile, params)
    primary_storage = processor.get_primary_storage()

    min_reference_length = 200000
    min_fragment_length = 2000
    min_fragment_reads = 1
    min_gap_length = 1000

    logging.info("Min reference length: {}".format(min_reference_length))
    logging.info("Min fragment length: {}".format(min_fragment_length))
    logging.info("Min reads in a fragment: {}".format(min_fragment_reads))

    read_length = processor.stats.read_length
    logging.info("Reading stats")
    stats_extractor = bam_stats_extractor.BamStatsExtractor(min_reference_length=min_reference_length,
                                                            read_length=read_length,
                                                            min_fragment_length=min_fragment_length,
                                                            min_fragment_reads=min_fragment_reads,
                                                            min_gap_length=min_gap_length)
    stats_list = stats_extractor.extract_secondary_stats(primary_storage)

    logging.info("Saving stats")
    stats_printer = bam_stats_io.BamStatsPrinter(saves_name)
    for stats in stats_list:
        stats_printer.print_bam_stats(stats)

    logging.info("Drawing plots")

    stats_visualizer = bam_stats_visualizer.OldBamStatsVisualizer("green", pictures_name)
    stats_visualizer.draw_statistics(stats_list)


def get_graph_names(graphname, prefix):
    subpaths = [path for path in os.listdir(graphname) if path.startswith(prefix)]
    print("Subpaths: ", subpaths)
    stats_paths = [os.path.join(graphname, subpath) for subpath in subpaths]
    return stats_paths


def draw_multi_distance_stats(path_to_graph_stats, picture_path):
    cluster_pictures_output_path = os.path.join(picture_path, "cluster_statistics")
    stats_paths = get_graph_names(path_to_graph_stats, "distance")
    print(stats_paths)
    graph_reader = graph_stats_reader.ClusterMultiStatsReader(stats_paths=stats_paths,
                                                              structure=graph_stats_reader.get_contracted_cluster_structure())
    dist_to_stats = graph_reader.read_dist_to_stats()
    small_range = (2000, 11000)
    big_range = (20000, 51000)
    print("Drawing plots for small distances")
    multivisualizer = graph_stats_visualizer.ClusterStatsMultiVisualizer(output_path=picture_path, dist_range=small_range)
    multivisualizer.draw_multi_stats(dist_to_stats, output_suffix="_small")
    print("Drawing plots for large distances")
    multivisualizer = graph_stats_visualizer.ClusterStatsMultiVisualizer(output_path=picture_path, dist_range=big_range)
    multivisualizer.draw_multi_stats(dist_to_stats, output_suffix="_big")


def draw_graph_stats(graphname, picture_path):
    graph_stats_output_path = os.path.join(picture_path, "graph_stats")
    if not os.path.exists(graph_stats_output_path):
        os.mkdir(graph_stats_output_path)
    graph_reader = graph_stats_reader.GraphStatsReader(graphname)
    stats = graph_reader.load_data()
    for son in stats.get_leaves():
        print(son.get_name())
    initial_filter_visualizer = graph_stats_visualizer.GraphStatsVisualizer(output_path=graph_stats_output_path)
    initial_filter_visualizer.draw_stats(stats)


def process_graph_stats(path_to_graph_stats, outname):
    saves_path, picture_path = prepare_output_dirs(output_name=outname)
    # draw_multi_distance_stats(path_to_graph_stats, picture_path)
    draw_graph_stats(path_to_graph_stats, picture_path)


if __name__ == "__main__":
    parser = createparser()
    args = parser.parse_args()
    outname = str(args.output)
    if os.path.exists(outname):
        shutil.rmtree(outname)

    os.mkdir(outname)
    if args.bamfile:
        bamname = args.bamfile
        tag = args.tag
        test = args.test
        distance = args.distance
        process_bamfile(bamname=bamname, tag=tag, test=test, distance=distance, output_name=outname)
        logging.info("BAM statistics extracted")

    if args.graph:
        process_graph_stats(path_to_graph_stats=args.graph, outname=outname)
