import logging
from collections import Counter
from collections import namedtuple
from operator import itemgetter

from bam_statistics.bam_stats import *


class FragmentFilter(object):
    def __init__(self, min_fragment_length, min_fragment_reads):
        self.min_fragment_length = min_fragment_length
        self.min_fragment_reads = min_fragment_reads

    def __call__(self, fragment):
        return fragment.len() >= self.min_fragment_length and fragment.get_reads() >= self.min_fragment_reads


class DelimeterReferenceGrouper(object):
    def __init__(self, delimeter):
        self.delimeter = delimeter

    def get_group(self, reference_name):
        group_list = reference_name.split(self.delimeter)
        return group_list[1]


class TrivialReferenceGrouper(object):
    def get_group(self, reference_name):
        return reference_name


class BamStatsExtractor(object):
    def __init__(self, min_reference_length, read_length, min_fragment_length, min_fragment_reads):
        self.min_reference_length = min_reference_length
        self.read_length = read_length
        self.min_fragment_length = min_fragment_length
        self.min_fragment_reads = min_fragment_reads

    def extract_secondary_stats(self, primary_storage):
        valid_barcode_storages = [ref_info.barcode_storage for ref_info in primary_storage.values() if
                                  ref_info.length > self.min_reference_length]
        logging.info('Valid references: {}'.format(len(valid_barcode_storages)))

        # Fragment stats
        fragment_filter = FragmentFilter(min_fragment_length=self.min_fragment_length,
                                         min_fragment_reads=self.min_fragment_reads)

        fragment_length = self.extract_fragment_length_distribution(valid_barcode_storages, fragment_filter=fragment_filter)
        fragment_per_barcode = self.extract_fragments_per_barcode(valid_barcode_storages)
        internal_coverage = self.extract_internal_coverage_distribution(valid_barcode_storages,
                                                                        self.read_length,
                                                                        fragment_filter)
        self.extract_read_per_fragment(valid_barcode_storages, self.read_length, self.min_fragment_length)
        reads_per_barcode = self.extract_read_per_barcode(valid_barcode_storages)

        # reference coverage
        # grouper = DelimeterReferenceGrouper("|")
        grouper = TrivialReferenceGrouper()
        reads_per_reference = self.extract_read_per_reference(primary_storage, grouper, self.read_length)
        fragments_per_reference = self.extract_fragments_per_reference(primary_storage, grouper, fragment_filter)

        return [fragment_length, fragment_per_barcode, internal_coverage,
                reads_per_barcode, reads_per_reference, fragments_per_reference]

    def extract_fragment_length_distribution(self, barcode_storages, fragment_filter):
        lengths = [fragment.len() for storage in barcode_storages
                   for barcode_entry in storage.values()
                   for fragment in barcode_entry.get_fragments() if fragment_filter(fragment)]
        logging.info("Extracting length distribution")
        length_counter = Counter(lengths)
        return FragmentLengthDistribution(dict(sorted(length_counter.items(), key=itemgetter(0))))

    def extract_fragments_per_barcode(self, barcode_storages):
        barcode_to_fragments = {}
        for storage in barcode_storages:
            for barcode in storage:
                if barcode not in barcode_to_fragments:
                    barcode_to_fragments[barcode] = storage.get_entry(barcode).get_number_of_fragments()
                else:
                    barcode_to_fragments[barcode] += storage.get_entry(barcode).get_number_of_fragments()
        counter = Counter(barcode_to_fragments.values())
        logging.info("Extracting fragments per barcode")
        # print(counter)
        return FragmentsPerBarcode(dict(sorted(counter.items(), key=itemgetter(0))))

    def extract_internal_coverage_distribution(self, barcode_storages, read_length, fragment_filter):
        coverages = [fragment.reads * read_length / float(fragment.len()) for storage in barcode_storages
                     for barcode_entry in storage.values()
                     for fragment in barcode_entry.get_fragments() if fragment_filter(fragment)]
        # print(coverages)
        logging.info("Internal coverage distribution: ")
        # print(len(coverages))
        # print("Max: ", max(coverages))
        # print("Min: ", min(coverages))
        # print("Mean: ", sum(coverages) / len(coverages))
        coverage_counter = Counter(coverages)
        return InternalCoverageDistribution(dict(sorted(coverage_counter.items(), key=itemgetter(0))))

    def extract_read_per_fragment(self, barcode_storages, read_length, min_fragment_length):
        reads = [fragment.reads for storage in barcode_storages
                 for barcode_entry in storage.values()
                 for fragment in barcode_entry.get_fragments() if fragment.len() > min_fragment_length]

        logging.info("Extracting read per fragment...")
        # print(len(reads))
        read_counter = Counter(reads)
        # print(sorted(read_counter.items(), key=itemgetter(0)))
        return ReadFragmentCoverage(dict(sorted(read_counter.items(), key=itemgetter(0))))

    def extract_read_per_barcode(self, barcode_storages):
        barcode_to_reads = {}
        for storage in barcode_storages:
            for barcode in storage:
                if barcode not in barcode_to_reads:
                    barcode_to_reads[barcode] = sum(
                        [fragment.reads for fragment in storage.get_entry(barcode).get_fragments()])
                else:
                    barcode_to_reads[barcode] += sum(
                        [fragment.reads for fragment in storage.get_entry(barcode).get_fragments()])

        logging.info("Extracting read per barcode...")
        # print(len(barcode_to_reads))
        counter = Counter(sorted(barcode_to_reads.values()))
        # print(sorted(counter.items(), key=itemgetter(0)))
        return ReadBarcodeCoverage(dict(sorted(counter.items(), key=itemgetter(0))))

    def extract_read_per_reference(self, primary_storage, grouper, read_length):
        ReadsLength = namedtuple('ReadsLength', ['reads', 'length'])
        group_to_info = {}
        for ref_info in primary_storage.values():
            # print(ref_info.name)
            group = grouper.get_group(ref_info.name)
            if group not in group_to_info:
                group_to_info[group] = ReadsLength(reads=0, length=0)
            storage = ref_info.barcode_storage
            num_reads = sum([fragment.get_reads() for barcode_entry in storage.values()
                             for fragment in barcode_entry.get_fragments()])

            group_to_info[group] = ReadsLength(reads=group_to_info[group].reads + num_reads,
                                               length=group_to_info[group].length + ref_info.length)
        # print(group_to_info)
        group_to_coverage = {ref: ((y.reads * read_length) / float(y.length))
                             for (ref, y) in group_to_info.items() if y.length > 200000}
        return ReadReferenceCoverage(dict(sorted(group_to_coverage.items(), key=itemgetter(0))))

    def extract_fragments_per_reference(self, primary_storage, grouper, fragment_filter):
        FragmentsLength = namedtuple('FragmentsLength', ['fragment_len', 'length'])
        group_to_info = {}
        for ref_info in primary_storage.values():
            # print(ref_info.name)
            group = grouper.get_group(ref_info.name)
            if group not in group_to_info:
                group_to_info[group] = FragmentsLength(fragment_len=0, length=0)
            storage = ref_info.barcode_storage
            fragment_length = sum([fragment.len() for barcode_entry in storage.values()
                                   for fragment in barcode_entry.get_fragments()])
            group_to_info[group] = FragmentsLength(fragment_len=group_to_info[group].fragment_len + fragment_length,
                                                   length=group_to_info[group].length + ref_info.length)
        # print(group_to_info)
        group_to_coverage = {ref: (y.fragment_len / float(y.length))
                             for (ref, y) in group_to_info.items() if y.length > 200000}
        return FragmentReferenceCoverage(dict(sorted(group_to_coverage.items(), key=itemgetter(0))))
