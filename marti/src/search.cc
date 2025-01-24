#include "marti/search.h"

#include <algorithm>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>
#include "edlib.h"

#include "marti/utils.h"


namespace marti {
    bool align_target(const seq_t target_name, const std::string& target_seq, const std::string& read_seq,
                      const std::optional<int> &read_start, const std::optional<int> &read_end,
                      const int max_dist, const bool score_only, TargetAlignment* aln) {
        *aln = {};
        aln->target_type = target_name;
        aln->target_seq = target_seq;
        aln->read_seq = read_seq;
        aln->score_only = score_only;
        const auto read_len = (read_end.has_value() ? read_end.value() : read_seq.size()) - read_start.value_or(0);
        aln->alignment = edlibAlign(target_seq.c_str(), (int) target_seq.size(),
                                   read_seq.c_str() + read_start.value_or(0), (int) read_len,
                                   edlibNewAlignConfig(max_dist, EDLIB_MODE_HW,
                                                       (score_only) ? EDLIB_TASK_DISTANCE : EDLIB_TASK_PATH,
                                                       nullptr, 0));
        if (aln->alignment.status != EDLIB_STATUS_OK) {
            throw std::runtime_error("[Search] edlib error during target alignment.");
        }
        aln->valid = aln->alignment.editDistance >= 0;
        if (!aln->valid) return false;
        aln->score = aln->alignment.editDistance;
        aln->end = aln->alignment.endLocations[0] + read_start.value_or(0);
        if (!score_only) {
            aln->start = aln->alignment.startLocations[0] + read_start.value_or(0);
        }
        return true;
    }

    bool find_polya(const seq_t target_name, const std::string &read_seq, const int max_dist, const int min_polya_len,
                    std::vector<TargetAlignment>& alns) {
        const char base = target_name == kPolyA ? 'A' : 'T';
        std::vector<int> polyx_mask(read_seq.length());
        int n_hits = std::count(read_seq.begin(), read_seq.begin() + min_polya_len, base);
        int n_errors = min_polya_len - n_hits;
        bool has_polya = false;
        if (n_hits + max_dist >= min_polya_len) {
            std::fill_n(polyx_mask.begin(), min_polya_len, 1);
            has_polya = true;
        }
        int prev = -1;
        for (auto i = min_polya_len; i < read_seq.length(); i++) {
            if (read_seq[i] == base) n_hits++;
            else n_errors++;
            if (read_seq[i - min_polya_len] == base) n_hits--;
            else n_errors--;

            if (n_hits + max_dist >= min_polya_len) {
                polyx_mask[i] = 1;
                has_polya = true;
                if (prev < 0 || prev != i - 1) {
                    for (auto j = 0; j < min_polya_len; j++)
                        polyx_mask[i - j] = 1;
                }
                prev = i;
            }
        }
        if (!has_polya) return false;
        for (auto i = 0; i < polyx_mask.size(); i++) {
            if (polyx_mask[i] && (i == 0 || !polyx_mask[i - 1])) {
                auto j = i;
                while (read_seq[j] != base) {
                    polyx_mask[j] = 0;
                    j++;
                }
            } else if (polyx_mask[i] && (i == polyx_mask.size() - 1 || !polyx_mask[i + 1])) {
                auto j = i;
                while (read_seq[j] != base) {
                    polyx_mask[j] = 0;
                    j--;
                }
            }
        }

        // extract the polyA runs from the mask into a list of alignments
        std::vector<std::pair<int, int>> polya_runs;
        int run_start = -1;
        for (auto i = 0; i < polyx_mask.size() + 1; i++) {
            if (i < polyx_mask.size() && polyx_mask[i] == 1) {
                if (run_start == -1) {
                    run_start = i;
                }
            } else {
                // either reached the end of the read or the end of the polyA run
                if (run_start != -1) {
                    int run_end = i - 1;
                    // refine the polyA ends (2 consecutive errors in +/- 3bp window)
                    const int window = 3;
                    if (std::count(read_seq.begin() + run_start,
                                   read_seq.begin() + run_start + window,base) == 1) {
                        run_start += window;
                    }
                    if (std::count(read_seq.begin() + run_end + 1 - window,
                                   read_seq.begin() + run_end + 1,base) == 1) {
                        run_end -= window;
                    }

                    TargetAlignment poly_run;
                    const int run_len = run_end - run_start + 1;
                    const int n_match = std::count(read_seq.begin() + run_start,
                                                   read_seq.begin() + run_end + 1,
                                                   base);
                    poly_run.target_type = target_name;
                    poly_run.read_seq = read_seq;
                    poly_run.score = run_len - n_match;
                    poly_run.start = run_start;
                    poly_run.end = run_end;
                    poly_run.valid = true;
                    alns.push_back(poly_run);
                    run_start = -1;
                }
            }
        }
        return !alns.empty();
    }

    void TargetAlignment::display_alignment(std::ostream& os, bool show_full) const {
        const int line_len = 100;
        const int name_len = 25;
        std::string read_str = show_full ? std::string(read_seq.substr(0, start)) : "";
        std::string target_str = show_full ? std::string(start, '-') : "";
        std::string align_str = show_full ? std::string(start, '-') : "";
        auto read_idx = start;
        auto target_idx = 0;
        if (!kPolyTails.contains(target_type)) {
            for (auto i = 0; i < alignment.alignmentLength; i++) {
                read_str += ((alignment.alignment[i] == EDLIB_EDOP_INSERT) ? '-' : read_seq[read_idx++]);
                target_str += ((alignment.alignment[i] == EDLIB_EDOP_DELETE) ? '-' : target_seq[target_idx++]);
                if (alignment.alignment[i] == EDLIB_EDOP_MATCH) {
                    align_str += "|";
                } else if (alignment.alignment[i] == EDLIB_EDOP_MATCH) {
                    align_str += ".";
                } else {
                    align_str += "-";
                }
            }
        } else { // special case for polyA display
            const char base = target_type == kPolyA ? 'A' : 'T';
            for (auto i = start; i < end + 1; i++) {
                read_str += read_seq[i];
                target_str += base;
                align_str += read_seq[i] == base ? '|' : '.';
            }
        }
        if (show_full) {
            read_str += read_seq.substr(end + 1, read_seq.length());
            target_str += std::string((read_seq.length() - end), '-');
            align_str += std::string((read_seq.length() - end), '-');
        }
        for (int i = 0; i < read_str.length(); i += line_len) {
            int lend = std::min(i + line_len, static_cast<int>(read_str.length()));
            os << marti::utils::ljust("READ", name_len)
               << read_str.substr(i, lend - i) << "\n"
               << marti::utils::ljust(" ", name_len)
               << align_str.substr(i, lend - i) << "\n"
               << marti::utils::ljust(std::string(target_type), name_len)
               << target_str.substr(i, lend - i) << "\n\n";
        }
    }
}
