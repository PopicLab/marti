/** 
 * @brief Read simulation utilities for various known read artifacts.
 */

#ifndef MARTI_TESTS_SIM_H_
#define MARTI_TESTS_SIM_H_

#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <ranges>
#include <unordered_map>
#include <vector>

#include "marti/definitions.h"
#include "marti/utils.h"

namespace marti {

// A basic read simulator for proper and artefactual read structures
class ReadSimulator {
public:
    struct SimulatedRead {
        std::string seq; // read sequence
        int cdna_start; // start of the cDNA sequence (when applicable)
        int cdna_end; // end of the cDNA sequence (when applicable)
        std::vector<seq_t> structure; // ordered list of read segment types
        std::vector<std::pair<int, int>> structure_coords; // the start and end of each read sub-sequence
    };
    explicit ReadSimulator();
    SimulatedRead generate_read(const std::vector<seq_t>& read_structure) const;
    inline std::vector<class_t> get_read_types() const {
        auto keys = std::views::keys(read_structures_);
        return { keys.begin(), keys.end() };
    }
    inline std::vector<std::vector<seq_t>>
    get_read_structures(class_t read_type, seq_t tso_adapter, seq_t rt_adapter) const {
        return read_structures_.at(read_type).at({tso_adapter, rt_adapter});
    }
    std::string get_sequence(seq_t seq_type) const { return target2seq_.at(seq_type); }
    void set_error_rate(float err_rate) { error_rate_ = err_rate; }
    void set_cdna_len(int cdna_len) { cdna_len_ = cdna_len; }
    void set_polya_len(int polya_len) { polya_len_ = polya_len; }
    void set_sls_gap(int sls_gap_len) { sls_gap_len_ = sls_gap_len; }
    void set_internal_targets(const std::vector<seq_t>& targets) { internal_targets_ = targets; }

public:
    void populate_structures(seq_t tso_adapter, seq_t rt_adapter);
    std::string generate_dna_seq(int seq_len) const;
    std::string generate_polya(char base) const;
    std::string generate_tso(seq_t tso) const;
    std::string insert_errors(const std::string& seq) const;
    int cdna_len_; // length of the native cDNA
    int polya_len_; // length of the polyA
    int sls_gap_len_; // distance between the PCR adapter and the SLS within the TSO
    float error_rate_; // error rate to insert when generating known sequences
    std::vector<seq_t> internal_targets_;
    std::unordered_map<seq_t, std::string> target2seq_;
    std::unordered_map<class_t, std::map<std::pair<seq_t, seq_t>, std::vector<std::vector<seq_t>>>> read_structures_;
};
}  // namespace marti

#endif // #ifndef MARTI_TESTS_SIM_H_
