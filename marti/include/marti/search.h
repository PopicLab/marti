/**
 * Sequence alignment routines used for target search
 */

#ifndef MARTI_INCLUDE_MARTI_SEARCH_H_
#define MARTI_INCLUDE_MARTI_SEARCH_H_

#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include "edlib.h"

#include "marti/definitions.h"
#include "marti/utils.h"

namespace marti {
    // Alignment record of a target sequence (e.g. PCR adapter or polyA) in the read sequence
    struct TargetAlignment {
        seq_t target_type; // type of target sequence
        std::string_view target_seq; // target sequence
        std::string_view read_seq;   // read sequence
        bool valid{}; // set to true if the target was found in the read
        bool score_only{}; // set to true if only the max score was computed
        int score{}; // alignment score
        int start{}; // start coordinate of the target in the read
        int end{}; // end coordinate of the target in the read
        EdlibAlignResult alignment{}; // edlib alignment result

        ~TargetAlignment() {
            edlibFreeAlignResult(alignment);
        }
        void display_alignment(std::ostream& os, bool show_full) const;
    };
    bool align_target(seq_t target_name, const std::string& target_seq, const std::string& read_seq,
                      const std::optional<int> &read_start, const std::optional<int> &read_end,
                      int max_dist, bool score_only, TargetAlignment* aln);
    bool find_polya(seq_t target_name, const std::string& read_seq, int max_dist, int min_polya_len,
                    std::vector<TargetAlignment>& alns);
}

#endif //MARTI_INCLUDE_MARTI_SEARCH_H_
