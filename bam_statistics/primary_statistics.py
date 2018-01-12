import logging
import statistics
import numpy as np
from collections import Counter


class Fragment(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.reads = 1
        self.gaps = []

    def __str__(self):
        return "(" + str(self.start) + ", " + str(self.end) + ")"

    def __repr__(self):
        return "(" + str(self.start) + ", " + str(self.end) + ")"

    def len(self):
        return self.end - self.start

    def get_reads(self):
        return self.reads

    def update_end(self, new_end):
        self.end = new_end

    def append_gap(self, gap):
        self.gaps.append(gap)

    def increment_count(self):
        self.reads += 1

    def get_gaps(self):
        return self.gaps


class BarcodeEntry(object):
    def __init__(self):
        self.fragments = []

    def __iter__(self):
        return iter(self.fragments)

    def append_fragment(self, fragment):
        self.fragments.append(fragment)

    def get_number_of_fragments(self):
        return len(self.fragments)

    def get_last_fragment(self):
        return self.fragments[len(self.fragments) - 1]

    def update_last_fragment(self, fragment):
        self.fragments[len(self.fragments) - 1] = fragment

    def get_fragments(self):
        return self.fragments

    def get_reads(self):
        return sum([fragment.reads for fragment in self.fragments])


class BarcodeToEntryStorage(object):
    def __init__(self):
        self.barcode_to_entry = {}

    def __iter__(self):
        return iter(self.barcode_to_entry)

    def create_entry(self, barcode):
        self.barcode_to_entry[barcode] = BarcodeEntry()

    def get_entry(self, barcode):
        return self.barcode_to_entry[barcode]

    def has_entry(self, barcode):
        return barcode in self.barcode_to_entry

    def get_entries(self):
        return self.barcode_to_entry.values()

    def keys(self):
        return self.barcode_to_entry.keys()


class BarcodeIndexer(object):
    def __init__(self):
        self.storage = {}

    def get_index(self, barcode):
        if barcode in self.storage:
            return self.storage[barcode]
        size = len(self.storage)
        self.storage[barcode] = size
        return size


class ReferenceBinSetStorage(object):
    def __init__(self, bin_size, ref_length):
        list_size = ref_length // bin_size + 1
        self.pos_to_barcode_set = [Counter() for i in range(list_size)]
        self.bin_size = bin_size
        self.barcode_indexer = BarcodeIndexer()

    def get_barcode_set(self, left_pos, right_pos):
        left_bin = left_pos // self.bin_size
        right_bin = right_pos // self.bin_size + 1
        result = Counter()
        for counter in self.pos_to_barcode_set[left_bin:right_bin]:
            result |= counter
        return result

    def put_barcode(self, pos, barcode):
        barcode_index = self.barcode_indexer.get_index(barcode)
        bin = pos // self.bin_size
        self.pos_to_barcode_set[bin][barcode_index] += 1


class ReferenceInfo(object):
    def __init__(self, barcode_storage, bin_storage, length, name):
        self.name = name
        self.length = length
        self.barcode_storage = barcode_storage
        self.bin_storage = bin_storage


class PrimaryStorage(object):
    def __init__(self):
        self.ref_to_info = {}

    def __iter__(self):
        return iter(self.ref_to_info)

    def create_entry(self, barcode_storage, bin_storage, name, length):
        self.ref_to_info[name] = ReferenceInfo(barcode_storage=barcode_storage, bin_storage=bin_storage,
                                               name=name, length=length)

    def has_entry(self, ref_name):
        return ref_name in self.ref_to_info

    def get_entry(self, ref_name):
        return self.ref_to_info[ref_name]

    def get_entries(self):
        return self.ref_to_info.values()


class ProcessorStatistics(object):
    def __init__(self):
        self.unmapped = 0
        self.mapped = 0
        self.unbarcoded = 0
        self.read_length = 0
        self.insert_size = 0
        self.multiple = 0


class ProcessorParameters(object):
    def __init__(self, gap_threshold, tag, test_mode):
        self.gap_threshold = gap_threshold
        self.tag = tag
        self.test_mode = test_mode


class BamStatisticsProcessor(object):
    def __init__(self, bamfile, params):
        self.bamfile_ = bamfile
        self.params_ = params
        self.stats = ProcessorStatistics()

    def fill_processor_statistics(self):
        self.stats.mapped = self.bamfile_.mapped
        self.stats.unmapped = self.bamfile_.unmapped
        fetch_reads = 250000
        counter = 0
        read_length_sum = 0
        insert_sizes = []
        for read in self.bamfile_.fetch():
            if counter % 2 == 0:
                read_length_sum += read.query_length
                insert_sizes.append(abs(read.template_length))
            counter += 1
            if counter == fetch_reads:
                break
        self.stats.read_length = 2 * read_length_sum // counter
        self.stats.insert_size = statistics.median(insert_sizes)

    def get_primary_storage(self):
        self.fill_processor_statistics()

        bamfile = self.bamfile_
        reference_idx = {name: length for (name, length) in zip(bamfile.references, bamfile.lengths)}
        counter = 0
        barcodes = set()

        primary_storage = PrimaryStorage()
        barcode_tag = self.params_.tag
        for read in self.bamfile_.fetch():
            if read.has_tag(barcode_tag):
                barcode = read.get_tag(barcode_tag)
                if barcode not in barcodes:
                    barcodes.add(barcode)
                ref_name = read.reference_name
                if not primary_storage.has_entry(ref_name):
                    length = reference_idx[ref_name]
                    primary_storage.create_entry(barcode_storage=BarcodeToEntryStorage(),
                                                 bin_storage=ReferenceBinSetStorage(bin_size=1000, ref_length=length),
                                                 name=ref_name, length=length)
                ref_info = primary_storage.get_entry(ref_name)
                if not ref_info.barcode_storage.has_entry(barcode):
                    ref_info.barcode_storage.create_entry(barcode)
                bin_storage = ref_info.bin_storage
                self.process_read(primary_storage, read, barcode, bin_storage, ref_name)
            else:
                self.stats.unbarcoded += 1

            counter += 1
            if counter % 2000000 == 0:
                logging.info("{} reads processed.".format(counter))
            if self.params_.test_mode and counter == 500000:
                break

        logging.info('Total reads: {}'.format(self.stats.mapped + self.stats.unmapped))
        logging.info('Mapped: {}'.format(self.stats.mapped))
        logging.info('Unmapped: {}'.format(self.stats.unmapped))
        logging.info('Multiple: {}'.format(self.stats.multiple))
        logging.info('Not barcoded: {}'.format(self.stats.unbarcoded))
        logging.info('Containers: {}'.format(len(barcodes)))
        logging.info('Average read length: {}'.format(self.stats.read_length))
        logging.info('Average insert size: {}'.format(self.stats.insert_size))
        return primary_storage

    def process_read(self, primary_storage, read, barcode, bin_storage, ref_name):
        ref_info = primary_storage.get_entry(ref_name)
        barcode_entry = ref_info.barcode_storage.get_entry(barcode)
        gap_threshold = self.params_.gap_threshold
        start_pos = read.reference_start
        end_pos = read.reference_start + read.query_length
        bin_storage.put_barcode(barcode=barcode, pos=start_pos)

        if read.has_tag("XA"):
            self.stats.multiple += 1
            return

        logging.debug('Start pos: {}, end pos: {}'.format(start_pos, end_pos))

        if barcode_entry.get_number_of_fragments() == 0:
            logging.debug('No fragments')
            current_fragment = Fragment(start=start_pos, end=end_pos)
            barcode_entry.append_fragment(current_fragment)
            logging.debug('Added new fragment: {}\n'.format(current_fragment))
        else:
            current_fragment = barcode_entry.get_last_fragment()
            gap = start_pos - current_fragment.end
            logging.debug('Found fragment {} with gap {}'.format(current_fragment, gap))
            if gap < gap_threshold:
                current_fragment.update_end(end_pos)
                current_fragment.append_gap(gap)
                current_fragment.increment_count()
                barcode_entry.update_last_fragment(current_fragment)
                logging.debug('Updated current fragment: {}\n'.format(current_fragment))
            else:
                new_fragment = Fragment(start=start_pos, end=end_pos)
                logging.debug('Added new fragment: {}\n'.format(new_fragment))
                barcode_entry.append_fragment(new_fragment)
