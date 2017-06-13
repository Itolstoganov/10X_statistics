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
    parser.add_argument('-g', '--graph_statistics', help="Path to graph_statistics statistics")
    return parser


def prepare_output_dirs(output_name):
    if os.path.exists(output_name):
        shutil.rmtree(output_name)

    log_outname = os.path.join(output_name, "log")
    saves_outname = os.path.join(output_name, "stats")
    pictures_outname = os.path.join(output_name, "pictures")

    os.mkdir(output_name)
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
    min_fragment_reads = 15
    logging.info("Min reference length: {}".format(min_reference_length))
    logging.info("Min fragment length: {}".format(min_fragment_length))
    logging.info("Min reads in a fragment: {}".format(min_fragment_reads))

    read_length = processor.stats.read_length
    logging.info("Reading stats")
    stats_extractor = bam_stats_extractor.BamStatsExtractor(min_reference_length=min_reference_length,
                                                            read_length=read_length,
                                                            min_fragment_length=min_fragment_length,
                                                            min_fragment_reads=min_fragment_reads)
    stats_list = stats_extractor.extract_secondary_stats(primary_storage)

    logging.info("Saving stats")
    stats_printer = bam_stats_io.BamStatsPrinter(saves_name)
    for stats in stats_list:
        stats_printer.print_bam_stats(stats)

    logging.info("Drawing plots")

    stats_visualizer = bam_stats_visualizer.OldBamStatsVisualizer("green", pictures_name)
    stats_visualizer.draw_statistics(stats_list)


def process_graph_stats(graphname, outname):
    # scores, edge_gaps = tenx_stats_reader.read_score_along_path(os.path.join(graphname, "genome_path_statistics"))
    # edgestats = tenx_stats_reader.read_edges(os.path.join(graphname, "edge_statistics"))
    # covs, gaps = tenx_stats_reader.read_coverage_to_gap_distribution(graphname + "/length_based_gap_distribution")
    overall_gap_distribution = graph_stats_reader.read_graph_gap_distribution(graphname + "/overall_gap_distribution")
    gaps, coverages = graph_stats_reader.read_reliable_barcode_distribution_coverage(
        graphname + "/reliable_barcodes_coverage")

    # tenx_stats_visualizer.draw_pie_edges(outname, edgestats)
    # tenx_stats_visualizer.plot_cov_gap_distribution(covs, gaps, outname)
    graph_stats_visualizer.plot_graph_gap_distribution(overall_gap_distribution, outname)
    graph_stats_visualizer.plot_gap_to_reliable_coverage(gaps, coverages, outname)


if __name__ == "__main__":
    parser = createparser()
    args = parser.parse_args()
    outname = str(args.output)
    if args.bamfile:
        bamname = args.bamfile
        tag = args.tag
        test = args.test
        distance = args.distance
        process_bamfile(bamname=bamname, tag=tag, test=test, distance=distance, output_name=outname)

    if args.graph:
        process_graph_stats(graphname=args.graph, outname=args.outname)
