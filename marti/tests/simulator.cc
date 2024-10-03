#include "simulator.h"

#include <ranges>
#include <stdexcept>

#include "marti/definitions.h"
#include "marti/config.h"
#include "marti/search.h"

namespace marti {

ReadSimulator::ReadSimulator() {
    cdna_len_ = 500;
    polya_len_ = 20;
    sls_gap_len_ = 0;
    error_rate_ = 0.0;
    target2seq_ = {
            {marti::kMars, "GCGCACACAGCTACTTAACACCGAAGCAGTGGTATCAACGCAGAGGCGC"},
            {marti::kVenus, "GCGCTACTTGTGAAGATTAACACCGCTACACGACGCTCTTCCGGCGC"},
            {marti::kSmrtbell, "TTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAA"},
            {marti::kMarsSls, "GCGCCACGCAATGAAGTCGCAGGGTCGG"},
            {marti::kMarsSplitSls, "CGCCTACACGCAATGAAGTCGCAGGGTTGGC"},
            {marti::kVenusSls, "ATAGATTAGGACGATAGCGATAGAGAG"},
            {marti::kVenusSplitSls, "CGCATAGATTAGGACGATAGCGATAGAGAGCCC"}
    };
    std::vector<seq_t> targets{std::views::keys(target2seq_).begin(),
                               std::views::keys(target2seq_).end()};
    for (auto seq_type : targets) {
        target2seq_[kSeqToRC.at(seq_type)] = utils::reverse_complement(target2seq_[seq_type]);
    }
    for (const seq_t tso_pcr_adapter: {kMars, kVenus}) {
        for (const seq_t rt_pcr_adapter: {kMars, kVenus}) {
            populate_structures(tso_pcr_adapter, rt_pcr_adapter);
        }
    }
}

ReadSimulator::SimulatedRead
ReadSimulator::generate_read(const std::vector<seq_t>& read_structure) const {
    SimulatedRead sim_read;
    sim_read.structure = read_structure;
    for (const seq_t seq_type: sim_read.structure) {
        int prev_len = sim_read.seq.size();
        if (seq_type == kCdna) {
            sim_read.seq += generate_dna_seq(cdna_len_);
            sim_read.cdna_start = prev_len;
            sim_read.cdna_end = sim_read.seq.size() - 1;
        } else if (kPolyTails.contains(seq_type)) {
            sim_read.seq += insert_errors(generate_polya(seq_type == kPolyA ? 'A' : 'T'));
        } else if (kTSOs.contains(seq_type)) {
            sim_read.seq += insert_errors(generate_tso(seq_type));
        } else {
            sim_read.seq += insert_errors(target2seq_.at(seq_type));
        }
        sim_read.structure_coords.emplace_back(prev_len, sim_read.seq.size() - 1);
    }
    return sim_read;
}

void ReadSimulator::populate_structures(seq_t tso_adapter, seq_t rt_adapter) {
    const auto tso = kAdapterToTSO.at(tso_adapter);
    const auto tso_rc = kSeqToRC.at(tso);
    const auto rt_adapter_rc = kSeqToRC.at(rt_adapter);
    const auto sls = kAdapterToSLS.at(tso_adapter);
    const auto split_sls = kAdapterToSplitSLS.at(tso_adapter);
    read_structures_[kProper][{tso_adapter, rt_adapter}] = {
          {tso, kCdna, kPolyA, rt_adapter_rc},
          {tso, kPolyA, kCdna, kPolyA, rt_adapter_rc},
          {tso, kPolyT, kCdna, kPolyA, rt_adapter_rc}};
    read_structures_[kTsoTso][{tso_adapter, rt_adapter}] = {{tso, kCdna, tso_rc}};
    read_structures_[kRtRt][{tso_adapter, rt_adapter}] = {{rt_adapter, kPolyT, kCdna, kPolyA, rt_adapter_rc}};
    read_structures_[kNoPolyA][{tso_adapter, rt_adapter}] = {{tso, kCdna, rt_adapter_rc}};
    read_structures_[kNoTSO][{tso_adapter, rt_adapter}] = {{tso_adapter, kCdna, kPolyA, rt_adapter_rc}};
    read_structures_[kOnlyPolyA][{tso_adapter, rt_adapter}] = {
          {tso,         kPolyA, rt_adapter_rc},
          {tso,         kPolyT, rt_adapter_rc},
          {tso,         kPolyA, tso_rc},
          {tso,         kPolyT, tso_rc},
          {tso_adapter, kPolyA, rt_adapter_rc},
          {tso_adapter, kPolyT, rt_adapter_rc},
          {rt_adapter,  kPolyA, rt_adapter_rc},
          {rt_adapter,  kPolyT, rt_adapter_rc}};
    read_structures_[kMissingAdapter][{tso_adapter, rt_adapter}] = { // missing one or both adapters
         {kCdna,  kPolyA, rt_adapter_rc},
         {kCdna,  rt_adapter_rc},
         {kPolyT, kCdna,  kPolyT, rt_adapter_rc},
         {kPolyT, kCdna,  kPolyA, rt_adapter_rc},
         {kPolyA, kCdna,  kPolyA, rt_adapter_rc},
         {kPolyA, kCdna,  kPolyT, rt_adapter_rc},
         {tso,    kCdna,  kPolyA},
         {tso,    kCdna},
         {tso,    kPolyA, kCdna,  kPolyA},
         {tso,    kPolyA, kCdna,  kPolyT},
         {tso,    kPolyT, kCdna,  kPolyA},
         {tso,    kPolyT, kCdna,  kPolyT},
         {kCdna,  kPolyA},
         {kCdna,  kPolyT},
         {kCdna},
         {kPolyA, kCdna,  kPolyA}};
    read_structures_[kRetainedSmrtbell][{tso_adapter, rt_adapter}] = {
        {tso, kCdna, kSmrtbell, kCdna, kPolyA, rt_adapter_rc},
        {tso_adapter, kCdna, kSmrtbell, kCdna, kPolyA, rt_adapter_rc},
        {tso_adapter, kCdna, kSmrtbell, kCdna, kPolyT, rt_adapter},
        {tso_adapter, kPolyA, kSmrtbell, kCdna, kPolyA, rt_adapter_rc},
        {tso_adapter, kPolyA, kCdna, kSmrtbell, kPolyT, tso_rc}};

    for (const auto& adapter: internal_targets_) {
        read_structures_[kInternalAdapters][{tso_adapter, rt_adapter}].insert(
                read_structures_[kInternalAdapters][{tso_adapter, rt_adapter}].end(),
                {{tso, kCdna, adapter, kCdna, kPolyA, rt_adapter_rc},
                {rt_adapter, kCdna, adapter, kCdna, kPolyA, rt_adapter_rc},
                {rt_adapter, kCdna, adapter, kCdna, kPolyT, rt_adapter_rc},
                {tso_adapter, kCdna, adapter, kCdna, rt_adapter_rc},
                {tso, kCdna, adapter, kCdna, tso_adapter},
                {tso, kPolyA, kCdna, adapter, kCdna, rt_adapter_rc},
                {tso, kPolyA, kCdna, adapter, kCdna, tso_rc},
                {tso, kPolyT, kCdna, adapter, kCdna, kPolyA, rt_adapter_rc},
                {tso, kPolyA, kCdna, adapter, kCdna, kPolyT, rt_adapter_rc},
                {tso, kPolyT, kCdna, adapter, kCdna, adapter, kCdna, kPolyA, rt_adapter_rc}});
    }
}

std::string ReadSimulator::generate_dna_seq(int seq_len) const {
    static constexpr auto dna = "ACGT";
    auto char_dist = std::uniform_int_distribution{0, 3};
    std::string read(seq_len, {});
    for (int i = 0; i < seq_len; ++i) {
        read[i] = dna[i % 4];
    }
    return read;
}

std::string ReadSimulator::generate_polya(char base) const {
    return std::string(polya_len_, base);
}

std::string ReadSimulator::generate_tso(const seq_t tso) const {
    const auto tso_adapter = kTSOToAdapter.at(tso);
    const bool split_sls = sls_gap_len_ > 0;
    const auto sls = split_sls ? kAdapterToSplitSLS.at(tso_adapter) : kAdapterToSLS.at(tso_adapter);
    if (kRCTSOs.contains(tso)) {
        return target2seq_.at(sls)
               + (split_sls ? generate_dna_seq(sls_gap_len_) : "")
               + target2seq_.at(tso_adapter);
    }
    return target2seq_.at(tso_adapter)
           + (split_sls ? generate_dna_seq(sls_gap_len_) : "")
           + target2seq_.at(sls);
}

std::string ReadSimulator::insert_errors(const std::string& seq) const {
    const std::map<char, char> substitutions = {{'A', 'C'}, {'C', 'T'}, {'G', 'A'}, {'T', 'G'}};
    const int n_errors = (std::floor(seq.size() * error_rate_));
    if (!n_errors) return seq;
    auto modified_seq = seq;
    // simple spaced substitution errors
    for (int i = 0; i < n_errors; i++) {
        const int err_idx = i * (seq.size()/n_errors);
        modified_seq[err_idx] = substitutions.at(seq[err_idx]);
    }
    return modified_seq;
}
}  // namespace marti
