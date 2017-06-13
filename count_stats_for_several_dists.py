#!/usr/bin/env python3
from joblib import Parallel, delayed
import os
import argparse
import process_stats


def createparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile', help="Barcoded BAM File")
    parser.add_argument('-t', '--tag', help="Barcode tag name", default="BX:Z")
    parser.add_argument('-o', '--output', help="Output base directory")
    parser.add_argument('--test', help="Test mode: use only first 1000000 reads from BAM", action='store_true')
    return parser


def get_output_name(output_base, dist):
    return os.path.join(output_base, "dist_" + str(dist))


parser = createparser()
args = parser.parse_args()

bam_path = args.bamfile
output_base = args.output
tag = args.tag
dist_list = [500, 1000, 2500, 5000, 10000, 20000]
test = args.test
print("Getting stats for distances {}".format(dist_list))

Parallel(n_jobs=len(dist_list))(delayed(process_stats.process_bamfile)(bamname=bam_path, tag=tag,
                                                                       output_name=get_output_name(output_base, dist),
                                                                       test=test, distance=dist) for dist in dist_list)
