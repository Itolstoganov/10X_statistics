from collections import namedtuple
from copy import deepcopy


class CandidateInfo(object):
    def __init__(self, cand_info):
        self.shared = cand_info.shared
        self.barcodes = cand_info.barcodes
        self.cov = cand_info.cov
        self.score = cand_info.score
        self.distance = cand_info.distance


class InitialFilterExtractor(object):
    def __init__(self, stats):
        self.stats_ = stats

    def get_edge_info(self, edge):
        return self.stats_[edge].first_edge_stats

    def get_conjugate(self, edge):
        return self.stats_[edge].first_edge_stats.rc_id

    def get_next(self, edge):
        return self.stats_[edge].first_edge_stats.next_id

    def get_candidate_info(self, edge, other):
        if other not in self.stats_[edge].candidates:
            raise ValueError("{} was not found in candidates of {}".format(other, edge))
        return self.stats_[edge].candidates[other]

    def get_candidate_score(self, edge, other):
        return self.get_candidate_info(edge, other).score

    def get_candidate_distance(self, edge, other):
        return self.get_candidate_info(edge, other).distance

    def get_candidate_shared(self, edge, other):
        return self.get_candidate_info(edge, other).shared

    def check_candidate(self, edge, candidate):
        return candidate in self.stats_[edge].candidates

    def get_candidates(self, edge):
        return self.stats_[edge].candidates

    def __iter__(self):
        return iter(self.stats_.keys())


class ScoreErrorStatistics(object):
    def __init__(self):
        self.overall = 0
        self.next_less_than_max = 0
        self.next_not_reached = 0
        self.max_is_near = 0
        self.large_distance_to_next = 0

    def __repr__(self):
        return "Overall: {} \nNext not reached: {} \nNext less than max: {}, " \
               "\nMax is near: {}, \nLarge distance to next: {}".format(self.overall, self.next_not_reached,
                                                                        self.next_less_than_max, self.max_is_near,
                                                                        self.large_distance_to_next)


class TopCandidatesStatstics(object):
    def __init__(self):
        self.overall = 0
        self.top_are_not_near = 0
        self.next_not_in_top = 0
        self.top_are_not_close = 0
        self.next_not_reached = 0

    def __repr__(self):
        return "Overall: {} \nTop are not near: {}\nTop are not close: {}"\
               "\nNext not in top: {}\nNext not reached: {}".format(self.overall, self.top_are_not_near, self.top_are_not_close,
                                                                    self.next_not_in_top, self.next_not_reached)


class InitialFilterProcessor(object):
    def __init__(self, extractor):
        self.extractor = extractor

    def get_next_edge_map(self):
        path_map = {edge: self.extractor.get_next(edge) for edge in self.extractor if
                    self.extractor.get_next(edge) != 0}
        return path_map

    def are_edges_nearby(self, next_map, first, second, distance=1):
        current = first
        for _ in range(distance):
            next = next_map[current]
            if second == next:
                return True
            current = next
        current = second
        for _ in range(distance):
            next = next_map[current]
            if first == next:
                return True
            current = next
        return False

    def is_edge_near_on_either_strand(self, next_map, edge, conjugate, candidate, distance=1):
        if not self.extractor.get_next(candidate):
            return False
        is_max_near_to_edge = self.are_edges_nearby(next_map=next_map, first=edge, second=candidate, distance=distance)
        is_max_near_to_conjugate = self.are_edges_nearby(next_map=next_map, first=conjugate, second=candidate,
                                                         distance=distance)
        return is_max_near_to_edge or is_max_near_to_conjugate

    def extract_neighbours(self, next_map, edge, distance):
        conjugate = self.extractor.get_conjugate(edge)
        neighbours = [candidate for candidate in self.extractor.get_candidates(edge)
                      if self.is_edge_near_on_either_strand(next_map, edge, conjugate, candidate, distance)]
        return neighbours

    def extract_top_candidates(self, edge, number_of_candidates):
        candidates = self.extractor.get_candidates(edge).keys()
        sorted_candidates = sorted(candidates,
                                   key=lambda candidate: -self.extractor.get_candidate_score(edge, candidate))
        return sorted_candidates[:number_of_candidates]

    def get_top_candidates_statistics(self, top_threshold):
        stats = TopCandidatesStatstics()
        extractor = self.extractor
        next_map = self.get_next_edge_map()
        for edge in extractor:
            next = extractor.get_next(edge)
            if next != 0:
                stats.overall += 1
                close_neighbours = self.extract_neighbours(next_map=next_map, edge=edge, distance=1)
                near_neighbours = self.extract_neighbours(next_map=next_map, edge=edge, distance=5)
                top_candidates = self.extract_top_candidates(edge, min(top_threshold, len(near_neighbours)))
                top_not_in_near = False
                for candidate in top_candidates:
                    if candidate not in near_neighbours:
                        # print("Edge: {}".format(edge))
                        # print("Next: {}".format(next))
                        # print("Top candidates: {}".format(top_candidates))
                        # print("Close neighbours: {}".format(close_neighbours))
                        # print("Near neighbours: {}".format(near_neighbours))
                        top_not_in_near = True
                        stats.top_are_not_near += 1
                        break
                for candidate in top_candidates:
                    if candidate not in close_neighbours:
                        stats.top_are_not_close += 1
                        break
                if next not in extractor.get_candidates(edge):
                    stats.next_not_reached += 1
                    continue
                if next not in top_candidates:
                    stats.next_not_in_top += 1
                    # if top_not_in_near:
                        # self._print_next_debug_info(edge, top_candidates, near_neighbours)
        return stats

    def get_next_score_distribution(self):
        extractor = self.extractor
        next_map = self.get_next_edge_map()
        distribution = [extractor.get_candidate_score(edge, next_map[edge])
                        for edge in self.extractor if
                        (edge in next_map and next_map[edge] in extractor.get_candidates(edge))]
        return distribution

    def _print_next_debug_info(self, edge, top_candidates, near_neighbours):
        extractor = self.extractor
        print("Edge: {}".format(edge))
        print("Edge Info: {}".format(extractor.get_edge_info(edge)))
        print("Next: {}".format(next))
        print("Info: {}".format(extractor.get_candidate_info(edge, next)))
        top = top_candidates[0]
        print("Top: {}".format(top))
        print("Top info: {}".format(extractor.get_candidate_info(edge, top)))
        print("Top is near: {}".format(top_candidates[0] in near_neighbours))
        last_top = top_candidates[len(top_candidates) - 1]
        print("Last top: {}".format(last_top))
        print("Min top info: {}".format(extractor.get_candidate_info(edge, last_top)))
        print("Min top is near: {}".format(last_top in near_neighbours))
        print("\n")

    def get_score_error_statistics(self):
        result = ScoreErrorStatistics()
        extractor = self.extractor
        next_map = self.get_next_edge_map()
        for edge in extractor:
            self._process_candidate(edge=edge, next_map=next_map, result=result)
        return result

    def _process_candidate(self, edge, next_map, result):
        extractor = self.extractor
        next_id = extractor.get_next(edge)
        distance_threshold = 50000
        next_score = 0
        if next_id == 0:
            return
        result.overall += 1
        if extractor.check_candidate(edge, next_id):
            next_score = extractor.get_candidate_info(edge, next_id).score
        else:
            result.next_not_reached += 1
            return
        next_distance = extractor.get_candidate_info(edge, next_id).distance
        max_score_info = max(
            [(cand_id, cand_value.score) for (cand_id, cand_value) in extractor.get_candidates(edge).items()],
            key=lambda cand_tuple: cand_tuple[1])
        max_edge = max_score_info[0]
        max_score = max_score_info[1]
        if max_score > next_score:
            result.next_less_than_max += 1
            conjugate = extractor.get_conjugate(edge)
            if self.is_edge_near_on_either_strand(next_map=next_map, edge=edge, conjugate=conjugate,
                                                  candidate=max_edge):
                result.max_is_near += 1
                return
            if next_distance > distance_threshold:
                result.large_distance_to_next += 1
