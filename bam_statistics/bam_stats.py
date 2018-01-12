import sys


class BamStats(object):
    def __init__(self, name, data):
        self.name_ = name
        self.data_ = data

    def get_data(self):
        return self.data_

    def get_name(self):
        return self.name_


class FragmentLengthDistribution(BamStats):
    # list of lengths
    def __init__(self, lengths):
        super(FragmentLengthDistribution, self).__init__(name="fragment_length_distribution",
                                                         data=lengths)


class FragmentsPerContainer(BamStats):
    # number of fragments:number of such barcodes
    def __init__(self, fragments_per_barcode):
        super(FragmentsPerContainer, self).__init__(name="fragments_per_container",
                                                    data=fragments_per_barcode)


class InternalCoverageDistribution(BamStats):
    # list of internal coverages
    def __init__(self, internal_coverages):
        super(InternalCoverageDistribution, self).__init__(name="internal_coverage_distribution",
                                                           data=internal_coverages)


class ReadsPerFragment(BamStats):
    # number of reads: number of such fragments
    def __init__(self, reads_per_fragment):
        super(ReadsPerFragment, self).__init__(name="reads_per_fragment",
                                               data=reads_per_fragment)


class ReadsPerContainer(BamStats):
    # number of reads: number of such barcodes
    def __init__(self, reads_per_barcode):
        super(ReadsPerContainer, self).__init__(name="reads_per_container",
                                                data=reads_per_barcode)


class ReadReferenceCoverage(BamStats):
    # reference: sum length of reads / reference length
    def __init__(self, read_reference_coverage):
        super(ReadReferenceCoverage, self).__init__(name="read_reference_coverage",
                                                    data=read_reference_coverage)


class FragmentReferenceCoverage(BamStats):
    # reference: sum length of fragments / reference length
    def __init__(self, fragment_reference_coverage):
        super(FragmentReferenceCoverage, self).__init__(name="fragment_reference_coverage",
                                                        data=fragment_reference_coverage)


class CoveredReferencesPerBarcode(BamStats):
    def __init__(self, covered_references_per_barcode):
        super(CoveredReferencesPerBarcode, self).__init__(name="covered_references_per_barcode",
                                                          data=covered_references_per_barcode)


class GapDistribution(BamStats):
    def __init__(self, gap_distribution):
        super(GapDistribution, self).__init__(name="gap_distribution",
                                              data=gap_distribution)


class BarcodeToReads(BamStats):
    def __init__(self, pair_sizes):
        super(BarcodeToReads, self).__init__(name="barcode_to_reads", data=pair_sizes)


class ScoreDistributions(BamStats):
    def __init__(self, cov_to_score_distribution):
        super(ScoreDistributions, self).__init__(name="score_distributions", data=cov_to_score_distribution)



def get_class_from_stat_name(name):
    camel_case_name = ''.join(word.capitalize() for word in name.split("_"))
    return getattr(sys.modules[__name__], camel_case_name)
