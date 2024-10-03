/**
 * @brief Read classification tests using simulated reads.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <string>

#include "marti/definitions.h"
#include "marti/classifier.h"
#include "marti/config.h"
#include "marti/utils.h"
#include "simulator.h"

using namespace marti;
TEST_CASE("sim_reads_no_error") {
    ReadSimulator read_simulator;
    read_simulator.set_error_rate(0.0);
    MartiConfig::MartiParams marti_params({
            .terminal_adapter_search_buffer=100,
            .terminal_polyA_search_buffer=150,
            .min_polyA_match=20,
            .max_err_rate=0,
            .max_err_rate_polya=0,
    });
    marti_params.target2seq = {{kMars, read_simulator.get_sequence(kMars)},
                               {kVenus, read_simulator.get_sequence(kVenus)},
                               {kSmrtbell, read_simulator.get_sequence(kSmrtbell)},
                               {kMarsSls, ""},
                               {kMarsSplitSls, ""},
                               {kVenusSls, ""},
                               {kVenusSplitSls, ""}};
    SUBCASE("multi_tso_merged_sls_both") {
        read_simulator.set_sls_gap(0);
        read_simulator.set_internal_targets({kMars, kVenus, kMarsRC, kVenusRC});
        marti_params.target2seq[kMarsSls] = read_simulator.get_sequence(kMarsSls);
        marti_params.target2seq[kVenusSls] = read_simulator.get_sequence(kVenusSls);
        marti_params.tso_adapters = {std::string(kMars), std::string(kVenus)};
        marti_params.rt_adapters = {std::string(kMars), std::string(kVenus)};
    }

    SUBCASE("multi_tso_split_sls_both") {
        read_simulator.set_sls_gap(20);
        read_simulator.set_internal_targets({kMars, kVenus, kMarsRC, kVenusRC});
        marti_params.target2seq[kMarsSplitSls] = read_simulator.get_sequence(kMarsSplitSls);
        marti_params.target2seq[kVenusSplitSls] = read_simulator.get_sequence(kVenusSplitSls);
        marti_params.tso_adapters = {std::string(kMars), std::string(kVenus)};
        marti_params.rt_adapters = {std::string(kMars), std::string(kVenus)};
    }

    SUBCASE("single_tso_split_sls") {
        read_simulator.set_sls_gap(20);
        read_simulator.set_internal_targets({kMars, kVenus, kMarsRC, kVenusRC});
        marti_params.target2seq[kMarsSplitSls] = read_simulator.get_sequence(kMarsSplitSls);
        marti_params.target2seq[kVenusSplitSls] = read_simulator.get_sequence(kVenusSplitSls);
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
    }

    SUBCASE("single_adapter_split_sls") {
        read_simulator.set_sls_gap(20);
        read_simulator.set_internal_targets({kMars});
        marti_params.target2seq[kMarsSplitSls] = read_simulator.get_sequence(kMarsSplitSls);
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kMars)};
    }

    MartiConfig marti_config(marti_params, "");
    ReadClassifier read_classifier(marti_config);
    std::unordered_map<class_t, size_t> sim_reads_per_class;
    // 1. test structures with a unique class
    for (const auto& tso_adapter: marti_params.tso_adapters) {
        for (const auto& rt_adapter: marti_params.rt_adapters) {
            for (const class_t read_type: read_simulator.get_read_types()) {
                // set the number of reads with this type (accounting for fwd and rc)
                sim_reads_per_class[read_type] += 2*read_simulator.get_read_structures(read_type, tso_adapter, rt_adapter).size();
                // these settings should guarantee that the read will not be filtered out based on length
                read_simulator.set_cdna_len(2 * marti_config.params.terminal_polyA_search_buffer);
                read_simulator.set_polya_len(marti_config.params.min_polyA_match);
                if (read_type == kOnlyPolyA) { // guarantee that the resulting read is not too short
                    read_simulator.set_polya_len(2 * marti_config.params.terminal_polyA_search_buffer);
                }
                for (const auto& structure : read_simulator.get_read_structures(read_type, tso_adapter, rt_adapter)) {
                    auto sim_read = read_simulator.generate_read(structure);
                    for (bool is_rc: {false, true}) {
                        ReadClassificationRecord annot;
                        CHECK_FALSE(annot.has_type());
                        read_classifier.classify(!is_rc ? sim_read.seq :
                                                  utils::reverse_complement(sim_read.seq),
                                                          annot);
                        std::set<class_t> read_labels = {read_type};
                        auto structure_str = utils::implode(structure, "...");
                        CAPTURE(read_type);
                        CAPTURE(structure_str);
                        CAPTURE(annot.class_label);
                        CAPTURE(annot.structure_label);
                        CHECK_EQ(annot.classes.size(), 1);
                        CHECK((annot.classes == read_labels));
                        if (is_rc) continue;
                        if (read_type == kProper) {
                            if (structure.size() > 4) {
                                // proper read with an additional internal polyA
                                CHECK(annot.has_internal_polya);
                                CHECK_EQ(annot.structure_alignments.size(), 4);
                                CHECK_EQ(annot.structure.size(), 4);
                                CHECK_EQ(annot.terminal_polya_alignments.size(), 2);
                                CHECK_EQ(annot.terminal_adapter_alignments.size(), 2);
                                CHECK_EQ(annot.auxiliary_alignments.size(), 1);
                            } else {
                                // canonical proper read
                                CHECK((annot.structure_label == structure_str));
                                int allowed_deviation = 0;
                                CHECK(((annot.cdna.start >= sim_read.cdna_start - allowed_deviation) &&
                                         (annot.cdna.start <= sim_read.cdna_start + allowed_deviation)));
                                CHECK(((annot.cdna.end >= sim_read.cdna_end - allowed_deviation) &&
                                         (annot.cdna.end <= sim_read.cdna_end + allowed_deviation)));
                                CHECK_EQ(annot.auxiliary_alignments.size(), 0);
                                CHECK_EQ(annot.terminal_polya_alignments.size(), 1);
                                CHECK_EQ(annot.terminal_adapter_alignments.size(), 2);
                                CHECK_EQ(annot.structure_alignments.size(), 4);
                                CHECK_EQ(annot.structure.size(), 4);
                                CHECK(annot.is_proper());
                                CHECK(annot.has_type());
                                CHECK_FALSE(annot.is_proper_and_rc());
                            }
                        } else {
                            const auto extended_structures = {kRetainedSmrtbell, kInternalAdapters};
                            if (!utils::contains(extended_structures, read_type)) {
                                CHECK((annot.structure_label == structure_str));
                            }
                        }
                    }
                }
            }
        }
    }

    read_classifier.get_read_stats().display(std::cout);
    // check the final read counts
    for (const class_t read_type: read_simulator.get_read_types()) {
        // note: this test assumes that each read has a unique type as above
        CHECK_EQ(read_classifier.get_read_stats().class_counts.at(std::string(read_type)),
                 sim_reads_per_class[read_type]);
    }

    // 2. test structures with multiple class labels
    std::vector<std::set<class_t>> mix_labels = {{kMissingAdapter, kInternalAdapters},
                                                 {kMissingAdapter, kRetainedSmrtbell},
                                                 {kInternalAdapters, kRetainedSmrtbell}};
    std::vector<std::vector<seq_t>> mixed_structures = {{kMarsTso, kCdna, kMars, kCdna, kPolyA},
                                                        {kSmrtbell, kCdna, kPolyA, kMarsRC},
                                                        {kMarsTso, kCdna, kMarsTso, kSmrtbell, kCdna, kPolyT, kMarsRC}};
    for(int i = 0; i < mix_labels.size(); i++) {
        const auto& read_labels = mix_labels[i];
        const auto& structure = mixed_structures[i];
        auto sim_read = read_simulator.generate_read(structure);
        ReadClassificationRecord annot;
        read_classifier.classify(sim_read.seq, annot);
        CAPTURE(i);
        CHECK((annot.classes == read_labels));
    }

    // 3. test unknown structures
    std::vector<std::vector<seq_t>> unk_structures = {{kMarsTso, kPolyT, kCdna, kPolyA, kMarsTsoRC},
                                                      {kMarsTso, kCdna, kPolyA, kMarsTso},
                                                      {kMars, kCdna, kPolyA, kMars}};
    for(int i = 0; i < mix_labels.size(); i++) {
        const std::set<class_t> read_labels = {kUnk};
        const auto& structure = unk_structures[i];
        auto sim_read = read_simulator.generate_read(structure);
        ReadClassificationRecord annot;
        read_classifier.classify(sim_read.seq, annot);
        CAPTURE(i);
        CHECK((annot.classes == read_labels));
    }
}

TEST_CASE("sim_reads_exp_error") {
    ReadSimulator read_simulator;
    read_simulator.set_error_rate(0.1);
    read_simulator.set_cdna_len(500);
    read_simulator.set_polya_len(20);
    MartiConfig::MartiParams marti_params({
             .terminal_adapter_search_buffer=100,
             .terminal_polyA_search_buffer=150,
             .min_polyA_match=20,
             .max_err_rate=0.1,
             .max_err_rate_polya=0.1,
     });
    marti_params.target2seq = {{kMars,          read_simulator.get_sequence(kMars)},
                               {kVenus,         read_simulator.get_sequence(kVenus)},
                               {kSmrtbell,      read_simulator.get_sequence(kSmrtbell)},
                               {kMarsSls,       ""},
                               {kMarsSplitSls,  ""},
                               {kVenusSls,      ""},
                               {kVenusSplitSls, ""}};
    SUBCASE("multi_tso_merged_sls_both") {
        read_simulator.set_sls_gap(0);
        read_simulator.set_internal_targets({kMars, kVenus, kMarsRC, kVenusRC});
        marti_params.target2seq[kMarsSls] = read_simulator.get_sequence(kMarsSls);
        marti_params.target2seq[kVenusSls] = read_simulator.get_sequence(kVenusSls);
        marti_params.tso_adapters = {std::string(kMars), std::string(kVenus)};
        marti_params.rt_adapters = {std::string(kMars), std::string(kVenus)};
    }

    SUBCASE("multi_tso_split_sls_both") {
        read_simulator.set_sls_gap(20);
        read_simulator.set_internal_targets({kMars, kVenus, kMarsRC, kVenusRC});
        marti_params.target2seq[kMarsSplitSls] = read_simulator.get_sequence(kMarsSplitSls);
        marti_params.target2seq[kVenusSplitSls] = read_simulator.get_sequence(kVenusSplitSls);
        marti_params.tso_adapters = {std::string(kMars), std::string(kVenus)};
        marti_params.rt_adapters = {std::string(kMars), std::string(kVenus)};
    }

    MartiConfig marti_config(marti_params, "");
    ReadClassifier read_classifier(marti_config);
    for (const auto& tso_adapter: marti_params.tso_adapters) {
        for (const auto& rt_adapter: marti_params.rt_adapters) {
            auto read_type = kProper;
            for (const auto &structure : read_simulator.get_read_structures(read_type, tso_adapter, rt_adapter)) {
                auto sim_read = read_simulator.generate_read(structure);
                for (bool is_rc: {false, true}) {
                    std::set<class_t> read_labels = {read_type};
                    ReadClassificationRecord annot;
                    read_classifier.classify(!is_rc ? sim_read.seq :
                                             utils::reverse_complement(sim_read.seq),
                                             annot);
                    auto structure_str = utils::implode(structure, "...");
                    CAPTURE(structure_str);
                    CAPTURE(annot.class_label);
                    CAPTURE(annot.structure_label);
                    CHECK((annot.classes == read_labels));
                    if (!is_rc) {
                        if (structure.size() > 4) {
                            CHECK(annot.has_internal_polya);
                        } else {
                            CHECK((annot.structure_label == structure_str));
                        }
                    }
                }
            }
        }
    }
}
