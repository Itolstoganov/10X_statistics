import logging
from collections import Counter
from collections import namedtuple
from operator import itemgetter
from multiprocessing import Pool
from itertools import chain
from itertools import combinations
from bam_statistics.bam_stats import *
import random


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
    def __init__(self, min_reference_length, read_length, min_fragment_length, min_fragment_reads, min_gap_length):
        self.min_reference_length = min_reference_length
        self.read_length = read_length
        self.min_fragment_length = min_fragment_length
        self.min_fragment_reads = min_fragment_reads
        self.min_gap_length = min_gap_length

    def extract_secondary_stats(self, primary_storage):
        valid_barcode_storages = [ref_info.barcode_storage for ref_info in primary_storage.get_entries() if
                                  ref_info.length > self.min_reference_length]
        logging.info('Valid references: {}'.format(len(valid_barcode_storages)))

        # Fragment stats
        fragment_filter = FragmentFilter(min_fragment_length=self.min_fragment_length,
                                         min_fragment_reads=self.min_fragment_reads)

        fragment_length = self.extract_fragment_length_distribution(valid_barcode_storages,
                                                                    fragment_filter=fragment_filter)
        fragment_per_barcode = self.extract_fragments_per_barcode(valid_barcode_storages,
                                                                  fragment_filter=fragment_filter)
        internal_coverage = self.extract_internal_coverage_distribution(valid_barcode_storages,
                                                                        self.read_length,
                                                                        fragment_filter)
        self.extract_read_per_fragment(valid_barcode_storages, self.read_length, self.min_fragment_length)
        barcode_to_reads = self._extract_barcode_to_reads(valid_barcode_storages)
        reads_per_barcode = self.extract_read_per_barcode(barcode_to_reads.get_data())

        grouper = TrivialReferenceGrouper()
        reads_per_reference = self.extract_read_per_reference(primary_storage, grouper,
                                                              self.read_length, self.min_reference_length)
        fragments_per_reference = self.extract_fragments_per_reference(primary_storage, grouper,
                                                                       fragment_filter, self.min_reference_length)
        references_per_barcode = self.extract_references_per_barcode(primary_storage, self.min_reference_length)
        gap_distribution = self.extract_gap_distribution(valid_barcode_storages, self.min_gap_length)

        first_len = 5000
        second_len = 20000
        samples = 20000
        score_read_threshold = 1
        # close_score_distribution = self.extract_close_score_distribution(primary_storage, first_len, second_len,
        #                                                                  samples, score_read_threshold)
        # distant_score_distribution = self.extract_distant_score_distribution(primary_storage, first_len, second_len,
        #                                                                      samples, score_read_threshold)
        # score_distributions = ScoreDistributions((close_score_distribution, distant_score_distribution))
        return [fragment_length, fragment_per_barcode, internal_coverage, reads_per_barcode, reads_per_reference,
                fragments_per_reference, references_per_barcode, gap_distribution, barcode_to_reads]

    def extract_fragment_length_distribution(self, barcode_storages, fragment_filter):
        lengths = [fragment.len() for storage in barcode_storages
                   for barcode_entry in storage.get_entries()
                   for fragment in barcode_entry.get_fragments() if fragment_filter(fragment)]
        logging.info("Extracting length distribution")
        length_counter = Counter(lengths)
        return FragmentLengthDistribution(dict(sorted(length_counter.items(), key=itemgetter(0))))

    def extract_fragments_per_barcode(self, barcode_storages, fragment_filter):
        barcode_to_fragments = {}
        for storage in barcode_storages:
            for barcode in storage:
                if barcode not in barcode_to_fragments:
                    # barcode_to_fragments[barcode] = storage.get_entry(barcode).get_number_of_fragments()
                    barcode_to_fragments[barcode] = len([fragment for fragment in
                                                         storage.get_entry(barcode).get_fragments()
                                                         if fragment_filter(fragment)])
                else:
                    barcode_to_fragments[barcode] += len([fragment for fragment in
                                                         storage.get_entry(barcode).get_fragments()
                                                         if fragment_filter(fragment)])
        counter = Counter(barcode_to_fragments.values())
        logging.info("Extracting fragments per barcode")
        # print(counter)
        return FragmentsPerContainer(dict(sorted(counter.items(), key=itemgetter(0))))

    def extract_internal_coverage_distribution(self, barcode_storages, read_length, fragment_filter):
        coverages = [fragment.reads * read_length / float(fragment.len()) for storage in barcode_storages
                     for barcode_entry in storage.get_entries()
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
                 for barcode_entry in storage.get_entries()
                 for fragment in barcode_entry.get_fragments() if fragment.len() > min_fragment_length]

        logging.info("Extracting read per fragment...")
        # print(len(reads))
        read_counter = Counter(reads)
        # print(sorted(read_counter.items(), key=itemgetter(0)))
        return ReadsPerFragment(dict(sorted(read_counter.items(), key=itemgetter(0))))

    def extract_read_per_reference(self, primary_storage, grouper, read_length, ref_length_threshold):
        ReadsLength = namedtuple('ReadsLength', ['reads', 'length'])
        group_to_info = {}
        for ref_info in primary_storage.get_entries():
            # print(ref_info.name)
            group = grouper.get_group(ref_info.name)
            if group not in group_to_info:
                group_to_info[group] = ReadsLength(reads=0, length=0)
            storage = ref_info.barcode_storage
            num_reads = sum([fragment.get_reads() for barcode_entry in storage.get_entries()
                             for fragment in barcode_entry.get_fragments()])

            group_to_info[group] = ReadsLength(reads=group_to_info[group].reads + num_reads,
                                               length=group_to_info[group].length + ref_info.length)
        # print(group_to_info)
        group_to_coverage = {ref: ((y.reads * read_length) / float(y.length))
                             for (ref, y) in group_to_info.items() if y.length > ref_length_threshold}
        return ReadReferenceCoverage(dict(sorted(group_to_coverage.items(), key=itemgetter(0))))

    def extract_fragments_per_reference(self, primary_storage, grouper, fragment_filter, ref_length_threshold):
        FragmentsLength = namedtuple('FragmentsLength', ['fragment_len', 'length'])
        group_to_info = {}
        for ref_info in primary_storage.get_entries():
            # print(ref_info.name)
            group = grouper.get_group(ref_info.name)
            if group not in group_to_info:
                group_to_info[group] = FragmentsLength(fragment_len=0, length=0)
            storage = ref_info.barcode_storage
            fragment_length = sum([fragment.len() for barcode_entry in storage.get_entries()
                                   for fragment in barcode_entry.get_fragments()])
            group_to_info[group] = FragmentsLength(fragment_len=group_to_info[group].fragment_len + fragment_length,
                                                   length=group_to_info[group].length + ref_info.length)
        # print(group_to_info)
        group_to_coverage = {ref: (y.fragment_len / float(y.length))
                             for (ref, y) in group_to_info.items() if y.length > ref_length_threshold}
        return FragmentReferenceCoverage(dict(sorted(group_to_coverage.items(), key=itemgetter(0))))

    def extract_references_per_barcode(self, primary_storage, ref_length_threshold):

        barcode_to_references = {}
        for ref_info in primary_storage.get_entries():
            if ref_info.length > ref_length_threshold:
                for barcode in ref_info.barcode_storage.keys():
                    if ref_info.barcode_storage.get_entry(barcode).get_reads() > 30:
                        if barcode not in barcode_to_references:
                            barcode_to_references[barcode] = 1
                        else:
                            barcode_to_references[barcode] += 1
        # print(barcode_to_references)
        counter = Counter(barcode_to_references.values())
        return CoveredReferencesPerBarcode(dict(sorted(counter.items(), key=itemgetter(0))))

    def extract_gap_distribution(self, valid_barcode_storages, min_gap_length):
        gaps = [gap for storage in valid_barcode_storages for entry in storage.get_entries()
                for fragment in entry.get_fragments() for gap in fragment.get_gaps()
                if gap > self.min_gap_length]
        gaps_counter = Counter(gaps)
        return GapDistribution(dict(sorted(gaps_counter.items(), key=itemgetter(0))))

    def extract_read_per_barcode(self, barcode_to_reads):
        # print(len(barcode_to_reads))
        counter = Counter(sorted(barcode_to_reads.values()))
        # print(sorted(counter.items(), key=itemgetter(0)))
        return ReadsPerContainer(dict(sorted(counter.items(), key=itemgetter(0))))

    def _extract_barcode_to_reads(self, barcode_storages):
        barcode_to_reads = {}
        for storage in barcode_storages:
            for barcode in storage:
                if barcode not in barcode_to_reads:
                    barcode_to_reads[barcode] = sum(
                        [fragment.reads for fragment in storage.get_entry(barcode).get_fragments()])
                else:
                    barcode_to_reads[barcode] += sum(
                        [fragment.reads for fragment in storage.get_entry(barcode).get_fragments()])
        return BarcodeToReads(barcode_to_reads)

    def _extract_barcode_set(self, barcode_to_reads):
        return barcode_to_reads.keys()

    def _check_pair(self, barcode_tuple):
        counter = 0
        for i in range(len(barcode_tuple[0])):
            if barcode_tuple[0][i] != barcode_tuple[1][i]:
                counter += 1
                if counter > 1:
                    logging.info(counter)
                    return None
        logging.info(counter)
        return barcode_tuple

    def pool_filter(self, pool, func, candidates):
        return [c for c, keep in zip(candidates, pool.map(func, candidates)) if keep]

    def construct_score_distribution_from_ref(self, bin_storage, first_len, second_len, distance_ranges, samples,
                                              ref_length,
                                              read_threshold):
        score_distribution = []
        first_positions = random.sample(range(first_len + 10000, ref_length - first_len - 10000), samples)
        counter = 0
        for first_pos in first_positions:
            second_pos = self.choose_from_tuples(tuples=distance_ranges, min_element=0,
                                                 max_element=ref_length - second_len, offset=second_len + 5000, pos=first_pos)
            # print(second_pos)
            first_counter = bin_storage.get_barcode_set(first_pos, first_pos + first_len)
            second_counter = bin_storage.get_barcode_set(second_pos, second_pos + second_len)
            keys = [item[0] for item in (first_counter & second_counter).items()]
            score = len(
                [key for key in keys if first_counter[key] > read_threshold and second_counter[key] > read_threshold])
            cov_second = sum(second_counter.values()) / float(second_len)
            score_distribution.append((score, cov_second))
            if counter % 1000 == 0:
                print("Generated {} samples".format(counter))
            counter += 1
        return score_distribution

    def extract_score_distribution(self, primary_storage, first_len, second_len, distance_ranges, samples,
                                   read_threshold):
        result = []
        for ref_info in primary_storage.get_entries():
            if ref_info.length > 200000:
                score_distribution = self.construct_score_distribution_from_ref(ref_info.bin_storage, first_len,
                                                                                second_len,
                                                                                distance_ranges, samples,
                                                                                ref_info.length, read_threshold)
                result += score_distribution
        print(len(result))
        return result

    def extract_close_score_distribution(self, primary_storage, first_len, second_len, samples, read_threshold):
        close_ranges = [(0, 10000)]
        score_distribution = self.extract_score_distribution(primary_storage, first_len, second_len, close_ranges,
                                                             samples, read_threshold)
        return score_distribution

    def extract_distant_score_distribution(self, primary_storage, first_len, second_len, samples, read_threshold):
        # LARGE_NUMBER = 100000000
        # distant_ranges = [(-LARGE_NUMBER, -10000), (10000, LARGE_NUMBER)]
        distant_ranges = [(-50000, -25000), (50000, 75000)]
        score_distribution = self.extract_score_distribution(primary_storage, first_len, second_len,
                                                             distant_ranges, samples, read_threshold)
        return score_distribution

    def choose_from_tuples(self, tuples, min_element, max_element, offset, pos):
        if max_element < offset:
            raise AttributeError("Len of the fragment is larger than length of the reference")
        new_tuples = [(min(max_element - offset, max(min_element, tup[0] + pos)),
                       max(min_element, min(max_element, tup[1] + pos))) for tup in tuples]
        # print(new_tuples)
        index = random.randint(0, len(new_tuples) - 1)
        # print(new_tuples[index])
        result = random.randint(new_tuples[index][0], new_tuples[index][1])
        return result

