
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


class FragmentsPerBarcode(BamStats):
    # number of fragments:number of such barcodes
    def __init__(self, fragments_per_barcode):
        super(FragmentsPerBarcode, self).__init__(name="fragments_per_container",
                                                  data=fragments_per_barcode)


class InternalCoverageDistribution(BamStats):
    # list of internal coverages
    def __init__(self, internal_coverages):
        super(InternalCoverageDistribution, self).__init__(name="internal_coverage_distribution",
                                                           data=internal_coverages)


class ReadFragmentCoverage(BamStats):
    # number of reads: number of such fragments
    def __init__(self, reads_per_fragment):
        super(ReadFragmentCoverage, self).__init__(name="reads_per_fragment",
                                                   data=reads_per_fragment)


class ReadBarcodeCoverage(BamStats):
    # number of reads: number of such barcodes
    def __init__(self, reads_per_barcode):
        super(ReadBarcodeCoverage, self).__init__(name="reads_per_container",
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


name_to_class = {
    "fragment_length_distribution": FragmentLengthDistribution,
    "fragments_per_container": FragmentsPerBarcode,
    "internal_coverage_distribution": InternalCoverageDistribution,
    "reads_per_fragment": ReadFragmentCoverage,
    "reads_per_container": ReadBarcodeCoverage,
    "read_reference_coverage": ReadReferenceCoverage,
    "fragment_reference_coverage": FragmentReferenceCoverage
}
