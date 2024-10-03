
#ifndef MARTI_INCLUDE_MARTI_CLASSIFIER_H_
#define MARTI_INCLUDE_MARTI_CLASSIFIER_H_

#include <filesystem>
#include <iostream>
#include <fstream>
#include <iterator>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <ranges>

#include "marti/config.h"
#include "marti/search.h"
#include "marti/definitions.h"
#include "marti/utils.h"
#include "marti/io.h"

namespace marti {
    // Read classification record storing read parsing and classification results
    struct ReadClassificationRecord {
        std::set<class_t> classes; // list of assigned artifact class(es)
        std::string class_label; // label combining all the assigned read classes
        std::vector<seq_t> structure; // list of read sequence components (ordered by position in the read)
        std::vector<TargetAlignment*> structure_alignments; // unordered list of target sequence alignments found in the read
        std::string structure_label; // representation of read structure
        std::vector<std::unique_ptr<TargetAlignment>> alignments; // all processed alignments
        std::string valid_targets_label; // alphabetically ordered list of valid matches
        // adapter and polyA assignments at the terminals used for structure compilation
        std::unordered_map<ReadTerminal, TargetAlignment*> terminal_adapter_alignments;
        std::unordered_map<ReadTerminal, TargetAlignment*> terminal_split_sls_alignments;
        std::unordered_map<ReadTerminal, TargetAlignment*> terminal_polya_alignments;
        std::vector<TargetAlignment*> auxiliary_alignments; // additional (internal) adapter and (near-terminal) polyA matches
        TargetAlignment cdna; // cDNA segment info
        bool has_external_polya{}; // externally-flanking polyAs near PCR-based adapters (TODO: reporting)
        bool has_internal_polya{}; // internally-flanking polyAs near PCR-based adapters (TODO: reporting)

        std::string summarize();
        inline bool has_type() const { return !classes.empty(); }
        inline void assign_type(class_t t) { classes.insert(t); }
        inline void set_class_label() { class_label = utils::implode(classes, "_"); }
        inline void set_structure_label() { structure_label = utils::implode(structure, "..."); }
        inline void set_valid_targets_label() {
            std::ostringstream structure_ss;
            for (auto &target: get_valid_matches()) {
                structure_ss << (structure_ss.str().empty() ? "" : ",") << target;
            }
            valid_targets_label = structure_ss.str();
        }
        inline bool is_proper() const { return !classes.empty() && *classes.begin() == kProper; }
        inline bool cannot_be_classified() const { return !classes.empty() && (*classes.begin() == kTooShort ||
                                                                               *classes.begin() == kLowRq); }
        inline bool is_proper_and_rc() {
            if (!is_proper()) return false;
            const auto tso = kTSOs.contains(terminal_adapter_alignments.at(ReadTerminal::LEFT)->target_type) ?
                             terminal_adapter_alignments.at(ReadTerminal::LEFT)->target_type : terminal_adapter_alignments.at(
                            ReadTerminal::RIGHT)->target_type;
            return kRCTSOs.contains(tso);
        }
        inline std::set<seq_t> get_valid_matches() {
            std::set<seq_t> valid_matches;
            for (auto &aln: alignments) {
                if (!aln->valid) continue;
                valid_matches.insert(aln->target_type);
            }
            return valid_matches;
        }

        inline std::string get_structure_coords() {
            std::ostringstream structure_ss;
            auto prefix = "";
            for (const auto& aln : structure_alignments) {
                structure_ss << prefix << aln->target_type << ":" << aln->start << ":" << aln->end << ":"<< aln->score;
                prefix = ",";
            }
            return structure_ss.str();
        }

        inline std::string get_all_valid_coords() {
            std::ostringstream structure_ss;
            auto prefix = "";
            for (const auto& aln : alignments) {
                if (!aln->valid || aln->score_only) continue;
                structure_ss << prefix << aln->target_type << ":" << aln->start << ":" << aln->end << ":"<< aln->score;
                prefix = ",";
            }
            return structure_ss.str();
        }

        inline TargetAlignment* new_target_alignment() {
            alignments.emplace_back(std::make_unique<TargetAlignment>());
            return alignments.back().get();
        }
    };

    // Output classification statistics
    struct ReadStats {
        const int n_examples_to_record = 10;
        std::map<std::string, int> class_counts; // number of reads of a particular category (e.g. proper)
        std::map<std::pair<std::string, std::string>, int> structure_counts; // number of reads with a particular class and structure
        std::map<std::pair<std::string, std::string>, std::vector<std::string>> class_examples; // examples of read parsing results

        // Updates read counters with the given classification result
        void update(const std::string& read_seq, const std::string& read_name, ReadClassificationRecord& annotation);
        void generate_reports(const std::filesystem::path& report_dir) const;
        void generate_read_parsing_logs(const std::filesystem::path& report_dir) const;
        // Merges read and structure counters
        inline void merge(const ReadStats& other) {
            for (auto const& [class_label, n_reads]: other.class_counts) {
                class_counts[class_label] += n_reads;
            }
            for (auto const& [read_label, n_reads]: other.structure_counts) {
                structure_counts[read_label] += n_reads;
            }
            for (auto const& [read_label, examples]: other.class_examples) {
                const int n_examples = std::min(examples.size(), n_examples_to_record - class_examples[read_label].size());
                class_examples[read_label].insert(class_examples[read_label].end(), examples.begin(),
                                                  examples.begin() + n_examples);
            }
        }
        // Prints the read class counts
        inline void display(std::ostream& os) const {
            os << "*****ReadStats***** \n";
            os << marti::utils::implode_map(class_counts) << "\n";
            os << "******************* \n";
        }
    };

    class ReadClassifier {
    public:
        explicit ReadClassifier(const MartiConfig& marti_config): cfg_(marti_config) {}
        void classify_and_annotate(io::read_t* read);
        void classify(const std::string& read_seq, ReadClassificationRecord& out, const std::string& read_name = "NA");
        void add_tags(io::read_t *read, ReadClassificationRecord& out);
        [[nodiscard]] const ReadStats& get_read_stats() const { return stats_; }

    private:
        static void assemble_structure(const std::string& read_seq, ReadClassificationRecord& out);
        class_t classify_structure(std::vector<seq_t>& structure) const;
        class_t lookup_structure(std::vector<seq_t>& structure, seq_t tso_adapter, seq_t rt_adapter,
                                 bool check_proper_only) const;
        void assign_adapters_to_terminals(const std::string& read_seq, ReadClassificationRecord& out);
        void assign_polya_to_terminals(const std::string& read_seq, ReadClassificationRecord& out);
        bool find_internal_adapters(const std::string& read_seq, ReadClassificationRecord& out);
        bool find_sequencing_adapters(const std::string& read_seq, ReadClassificationRecord& out);
        bool can_reclassify_as_proper(ReadClassificationRecord& out);

        // Returns true if the target alignment is considered valid (e.g. within the expected edit distance)
        // For TSO: matches of the PCR and SLS sequences are also examined (expected in target_alignments)
        [[nodiscard]] inline bool is_valid_target_match(const seq_t& target,
                                          const std::unordered_map<seq_t, TargetAlignment*>& target_alignments) const {
            return target_alignments.at(target)->valid &&
                   (!kTSOs.contains(target) || is_valid_tso_match(target, target_alignments));
        }

        // Returns true if this is a valid TSO match (by comparison to the PCR adapter only match)
        [[nodiscard]] bool is_valid_tso_match(const seq_t& target, const std::unordered_map<seq_t,
                                              TargetAlignment*>& target_alignments) const;
        const MartiConfig& cfg_;
        ReadStats stats_;
    };
}

#endif //MARTI_INCLUDE_MARTI_CLASSIFIER_H_
