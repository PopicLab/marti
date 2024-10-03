#ifndef MARTI_INCLUDE_MARTI_CONFIG_H_
#define MARTI_INCLUDE_MARTI_CONFIG_H_

#include <iostream>
#include <optional>
#include <filesystem>
#include <fstream>
#include <set>
#include <yaml-cpp/yaml.h>

#include "marti/definitions.h"
#include "marti/io.h"

namespace marti {
    // Marti configuration
    struct MartiConfig {
        /////////// user-controlled settings
        struct MartiParams {
            std::string test_seq;  // read sequence to process directly (skips the BAM file)
            std::string input_bam; // input BAM file with cDNA reads
            /////////// read filters
            float min_rq{}; // minimum rq value required to classify a read
            std::optional<int> max_reads_to_process{}; // maximum number of reads to process
            /////////// target oligo config
            std::map<seq_t, std::string> target2seq; // maps target adapters to their corresponding sequences
            std::set<std::string> tso_adapters; // list of valid PCR adapters for the TSO
            std::set<std::string> rt_adapters; // list of valid PCR adapters for the RT primer
            /////////// read structure config
            int terminal_adapter_search_buffer{}; // length of terminal for adapter search
            int terminal_polyA_search_buffer{}; // length of terminal for polyA search
            int min_polyA_match{}; // minimum length of a polyA match
            float max_err_rate{}; // maximum error rate expected for sequence targets
            float max_err_rate_polya{}; // maximum error rate expected for polyA targets
            ////////// resource config
            int n_threads{};
            int n_max_reads_in_mem{}; // maximum number of reads loaded into memory (for batched processing)
            ////////// output config
            int num_to_record{}; // maximum number of reads to log for each artifact
            bool compute_all_internal_matches{}; // compute and show all adapter matches in the read
            bool output_proper_trimmed{};
            void load(const io::YAMLFile& config_file);
        };
        MartiParams params; // input user parameters
        /////////// derived settings
        std::vector<seq_t> pcr_targets; // PCR adapter based targets (ordered s.t. TSO matches are computed last)
        std::vector<seq_t> seq_targets; // Sequencing adapter targets (SMRTbell)
        std::vector<seq_t> polya_targets; // polyA targets
        std::map<seq_t, std::string> target2seq; // maps adapter-based targets to their corresponding sequences
        std::map<seq_t, int> targets2lev; // maps targets to their maximum allowed Levenshtein distance
        std::map<seq_t, bool> is_split_tso; // indicates whether the TSO is configured with a split SLS
        // maps a read class to its expected structure (order of target oligos)
        std::unordered_map<class_t, std::map<std::pair<seq_t, seq_t>, std::vector<seq_t>>> class2structure;
        // output folders and files
        std::filesystem::path experiment_dir;
        std::filesystem::path report_dir;
        std::filesystem::path output_bam; // output BAM file with annotated cDNA reads
        std::filesystem::path output_proper_bam; // output BAM file with trimmed proper cDNA reads
        std::filesystem::path log_file_path;
        std::ofstream log_file;

        explicit MartiConfig(MartiParams params, const std::string& output_path): params(std::move(params)),
        experiment_dir(output_path) { configure(); }

        explicit MartiConfig(const std::string& config_path) {
            io::YAMLFile config_file{config_path}; // YAML file handle
            params.load(config_file);
            experiment_dir = std::filesystem::path(config_path).parent_path();
            configure();
        }

        void configure() {
            verify();
            set_targets();
            set_structures();
            set_outputs();
            MARTI_LOG(log_file, summarize());
        }

        // Validates the user-provided parameters
        void verify() const;
        // Returns the configuration summary
        std::string summarize() const;
        // Determines the search targets based on the configured TSO and RT adapters and other expected oligos
        // (polyAs and sequencing adapters); populates the sequences and Lev distances for adapter-based targets
        void set_targets();
        // Populates all the sequence targets associated with a specific PCR adapter,
        // (e.g. the reverse complement sequences, TSO, or split SLS)
        void populate_pcr_targets(seq_t adapter, bool is_tso);
        // Populates all the predefined read structures
        inline void set_structures() {
            for (const seq_t tso_adapter: params.tso_adapters) {
                for (const seq_t rt_adapter: params.rt_adapters) {
                    class2structure[kProper][{tso_adapter, rt_adapter}] = {kAdapterToTSO.at(tso_adapter), kCdna, kPolyA, kSeqToRC.at(rt_adapter)};
                    class2structure[kTsoTso][{tso_adapter, rt_adapter}] = {kAdapterToTSO.at(tso_adapter), kCdna, kSeqToRC.at(kAdapterToTSO.at(tso_adapter))};
                    class2structure[kRtRt][{tso_adapter, rt_adapter}] = {rt_adapter, kPolyT, kCdna, kPolyA, kSeqToRC.at(rt_adapter)};
                    class2structure[kNoPolyA][{tso_adapter, rt_adapter}] = {kAdapterToTSO.at(tso_adapter), kCdna, kSeqToRC.at(rt_adapter)};
                    class2structure[kNoTSO][{tso_adapter, rt_adapter}] = {tso_adapter, kCdna, kPolyA, kSeqToRC.at(rt_adapter)};
                }
            }
        }
        // Sets up the output folders/files
        inline void set_outputs() {
            report_dir = experiment_dir / "reports";
            std::filesystem::create_directories(report_dir);
            output_bam = experiment_dir / (std::filesystem::path(params.input_bam).stem().string() + ".classified.bam");
            output_proper_bam = experiment_dir / (std::filesystem::path(params.input_bam).stem().string() + ".proper.cdna.bam");
            log_file_path = report_dir / "main.log";
            log_file.open(log_file_path);
        }
};

}
#endif // #ifndef MARTI_INCLUDE_MARTI_CONFIG_H_
