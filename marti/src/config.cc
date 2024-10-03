#include "marti/config.h"

#include <filesystem>
#include <iostream>
#include <string>
#include <ranges>

#include "marti/utils.h"

namespace marti {
    void MartiConfig::MartiParams::load(const io::YAMLFile& config_file) {
        // input/output reads
        try {
            config_file.load("test_seq", test_seq, {});
            config_file.load("input_bam", input_bam, {});
            config_file.load("min_rq", min_rq, 1.0f);
            int max_reads;
            config_file.load("max_reads_to_process", max_reads, {});
            if (max_reads) max_reads_to_process = max_reads;
            config_file.load("output_proper_trimmed", output_proper_trimmed, false);
            config_file.load("compute_all_internal_matches", compute_all_internal_matches, false);
            // search params
            config_file.load("terminal_adapter_search_buffer", terminal_adapter_search_buffer, 100);
            config_file.load("terminal_polyA_search_buffer", terminal_polyA_search_buffer, 150);
            config_file.load("min_polyA_match", min_polyA_match, 20);
            config_file.load("max_err_rate", max_err_rate, 0.1f);
            config_file.load("max_err_rate_polya", max_err_rate_polya, 0.1f);
            // target sequences
            config_file.load(kMars, target2seq[kMars], {});
            config_file.load(kMarsSls, target2seq[kMarsSls], {});
            config_file.load(kMarsSplitSls, target2seq[kMarsSplitSls], {});
            config_file.load(kVenus, target2seq[kVenus], {});
            config_file.load(kVenusSls, target2seq[kVenusSls], {});
            config_file.load(kVenusSplitSls, target2seq[kVenusSplitSls], {});
            config_file.load(kSmrtbell, target2seq[kSmrtbell],
                             std::string("TTCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAA"));
            // TSO/RT adapter assignments
            config_file.load("TSO_adapters", tso_adapters);
            config_file.load("RT_adapters", rt_adapters);
            // load resource params
            config_file.load("n_threads", n_threads, 1);
            config_file.load("n_max_reads_in_mem", n_max_reads_in_mem, 100000);
        } catch (const std::exception& e) {
            throw std::invalid_argument("[IO] YAML parsing error (check parameter types!). message=[" + std::string(e.what()) + "]");
        }
    }

    void MartiConfig::set_targets() {
        // 1. PCR adapters and derived targets
        for (const auto& tso_pcr_adapter: params.tso_adapters) {
            populate_pcr_targets(tso_pcr_adapter, true);
        }
        for (const auto& rt_pcr_adapter: params.rt_adapters) {
            if (!params.tso_adapters.contains(rt_pcr_adapter)) {
                populate_pcr_targets(rt_pcr_adapter, false);
            }
        }
        // 2. sequencing adapter targets
        for (seq_t target: kSeqAdapters) {
            if (!params.target2seq[target].empty()) {
                seq_targets.push_back(target);
                seq_targets.push_back(kSeqToRC.at(target));
                target2seq[target] = params.target2seq[target];
                target2seq[kSeqToRC.at(target)] = utils::reverse_complement(params.target2seq[target]);
                targets2lev[target] = std::ceil(target2seq[target].length() * params.max_err_rate);
                targets2lev[kSeqToRC.at(target)] = targets2lev[target];
            }
        }

        // 3. polyA targets
        for (seq_t target: kPolyTails) {
            polya_targets.push_back(target);
            targets2lev[target] = std::ceil(params.min_polyA_match * params.max_err_rate_polya);
        }
    }

    void MartiConfig::populate_pcr_targets(const seq_t adapter, const bool is_tso) {
        std::vector<seq_t> targets_fwd{adapter};
        target2seq[adapter] = params.target2seq[adapter];
        if (is_tso) {
            const auto tso = kAdapterToTSO.at(adapter);
            if (!params.target2seq.at(kAdapterToSplitSLS.at(adapter)).empty()) {
                // the split SLS is a separate target
                const auto split_sls = kAdapterToSplitSLS.at(adapter);
                target2seq[split_sls] = params.target2seq[split_sls];
                target2seq[tso] = target2seq[adapter];
                is_split_tso[tso] = true;
                is_split_tso[kSeqToRC.at(tso)] = true;
                targets_fwd.push_back(split_sls);
            } else {
                // the joined SLS is not a separate target, added to the TSO
                const auto sls = kAdapterToSLS.at(adapter);
                target2seq[tso] = target2seq[adapter] + params.target2seq[sls];
                is_split_tso[tso] = false;
                is_split_tso[kSeqToRC.at(tso)] = false;
                targets2lev[sls] = std::ceil(params.target2seq[sls].length() * params.max_err_rate);
                targets2lev[kSeqToRC.at(sls)] = targets2lev[sls];
            }
            targets_fwd.push_back(tso);
        }
        for (const auto& target : targets_fwd) {
            target2seq[kSeqToRC.at(target)] = utils::reverse_complement(target2seq.at(target));
            targets2lev[target] = std::ceil(target2seq[target].length() * params.max_err_rate);
            targets2lev[kSeqToRC.at(target)] = targets2lev[target];
            pcr_targets.push_back(target);
            pcr_targets.push_back(kSeqToRC.at(target));
        }
    }

    void MartiConfig::verify() const {
        if (!params.input_bam.empty() && !std::filesystem::exists(params.input_bam)) {
            throw std::invalid_argument("Input BAM file does not exist: " + params.input_bam);
        }
        if (params.tso_adapters.empty() || params.rt_adapters.empty()) {
            throw std::invalid_argument("At least one adapter must be specified in TSO_adapters and RT_adapters.");
        }
        // validate TSO adapters
        const std::set<seq_t> valid_adapters = {kMars, kVenus};
        for (const auto &adapter: params.tso_adapters) {
            // must be the expected adapter name
            if (!valid_adapters.contains(adapter)) {
                throw std::invalid_argument(
                        "TSO_adapters contains an unexpected value. Allowed entries: " + utils::implode(valid_adapters));
            }
            // sequence must be provided
            if (params.target2seq.at(adapter).empty()) {
                throw std::invalid_argument("Adapter sequence was not provided for: " + adapter);
            }
            // must have exactly one provided SLS sequence
            if (params.target2seq.at(kAdapterToSLS.at(adapter)).empty()
                && params.target2seq.at(kAdapterToSplitSLS.at(adapter)).empty()) {
                throw std::invalid_argument("No SLS was provided for the TSO adapter: " + adapter);
            } else if (!params.target2seq.at(kAdapterToSLS.at(adapter)).empty() &&
                       !params.target2seq.at(kAdapterToSplitSLS.at(adapter)).empty()) {
                throw std::invalid_argument("Two SLS options were provided for the TSO adapter (only one is allowed): " + adapter);
            } // TODO: should a mix of SLS types be allowed across adapters
            // must be shorter (adapter + SLS) than the search buffer
            const auto tso_len = params.target2seq.at(adapter).length() +
                                            params.target2seq.at(kAdapterToSLS.at(adapter)).length() +
                                            params.target2seq.at(kAdapterToSplitSLS.at(adapter)).length();
            if (params.terminal_adapter_search_buffer < tso_len) {
                throw std::invalid_argument("terminal_adapter_search_buffer is smaller than the longest TSO sequence");
            }
        }
        // validate RT adapters
        for (const auto &adapter: params.rt_adapters) {
            if (!valid_adapters.contains(adapter)) {
                throw std::invalid_argument("RT_adapters contains an unexpected value. Allowed entries: " +
                        utils::implode(valid_adapters));
            }
            if (params.target2seq.at(adapter).empty()) {
                throw std::invalid_argument("Adapter sequence was not provided for: " + adapter);
            }
            if (params.target2seq.at(adapter).length() > params.terminal_adapter_search_buffer) {
                throw std::invalid_argument(
                        "terminal_adapter_search_buffer is smaller than the longest RT adapter sequence");
            }
        }
        if (params.terminal_adapter_search_buffer > params.terminal_polyA_search_buffer) {
            throw std::invalid_argument(
                    "terminal_adapter_search_buffer must be smaller than terminal_polyA_search_buffer");
        }
        if (params.terminal_adapter_search_buffer > params.terminal_polyA_search_buffer) {
            throw std::invalid_argument(
                    "terminal_adapter_search_buffer must be smaller than terminal_polyA_search_buffer");
        }
        if (params.terminal_adapter_search_buffer < 0) {
            throw std::invalid_argument("terminal_adapter_search_buffer must be greater than 0.");
        }
        if (params.max_err_rate < 0 || params.max_err_rate_polya < 0) {
            throw std::invalid_argument("max error rate cannot be negative.");
        }
    }

    std::string MartiConfig::summarize() const {
        std::ostringstream s;
        s << "*********** MartiConfig ***********\n"
        << "----- Input -------\n"
        << "input_bam: " << (params.input_bam.empty() ? "(none)" : params.input_bam) << "\n"
        << "min_rq: " << params.min_rq << "\n"
        << "max_reads_to_process: " << (params.max_reads_to_process.has_value() ? std::to_string(
    params.max_reads_to_process.value()) : "(none)") << "\n"
        << "----- Target oligos -------\n"
        << "TSO adapter(s): " << utils::implode(params.tso_adapters) << "\n"
        << "RT adapter(s): " << utils::implode(params.rt_adapters) << "\n"
        << utils::implode_map(target2seq) << "\n"
        << "----- Search configuration -------\n"
        << "terminal_adapter_search_buffer: " << params.terminal_adapter_search_buffer << "\n"
        << "terminal_polyA_search_buffer: " << params.terminal_polyA_search_buffer << "\n"
        << "min_polyA_match: " << params.min_polyA_match << "\n"
        << "max_err_rate: " << params.max_err_rate << "\n"
        << "max_err_rate_polya: " << params.max_err_rate_polya << "\n"
        << "----- Maximum Levenshtein distance / oligo -------\n"
        << utils::implode_map(targets2lev) << "\n"
        << "------Runtime configuration-------\n"
        << "n_threads: " << params.n_threads << "\n"
        << "n_max_reads_in_mem: " << params.n_max_reads_in_mem << "\n"
        << "------Output-------\n"
        << "output_proper_bam: " << params.output_proper_trimmed << "\n"
        << "compute_all_internal_matches: " << params.compute_all_internal_matches << "\n"
        << "output_bam: " << output_bam << "\n"
        << "experiment_dir: " << experiment_dir << "\n"
        << "report_dir: " << report_dir << "\n"
        << "log_file_path: " << log_file_path << "\n"
        << "**********************" << std::endl;
        return s.str();
    }
}
