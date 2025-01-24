#include "marti/classifier.h"

#include <filesystem>
#include <stdexcept>

#include "marti/search.h"

namespace marti {

    //////////////////////////
    ////  Classification  ////
    //////////////////////////

    void ReadClassifier::classify_and_annotate(io::read_t* read) {
        ReadClassificationRecord result;
        const auto read_seq = io::get_seq(read);
        const auto read_name = io::get_name(read);
        // 1. apply tag-based filter(s) to determine if the read should be classified
        if (io::has_tag(read, "rq") && io::get_tag<double>(read, "rq") < cfg_.params.min_rq) {
            result.assign_type(kLowRq);
            stats_.update(read_seq, read_name, result);
        } else {
            // 2. search for targets in the read sequence and classify based on their presence/absence and order
	        classify(read_seq, result, read_name);
    	}
        // 3. generate and add new BAM tags
        add_tags(read, result);
    }

    void ReadClassifier::classify(const std::string& read_seq, ReadClassificationRecord& out,
                                  const std::string& read_name) {
        // 1. apply length-based filter(s)
        if (read_seq.length() < 2 * cfg_.params.terminal_polyA_search_buffer) {
            out.assign_type(kTooShort);
            stats_.update(read_seq, read_name, out);
            return; // cannot be parsed further
        }
        // 2. search for PCR-based adapters at the terminals
        assign_adapters_to_terminals(read_seq, out);
	    // 3. check if a terminal wasn't assigned to a target adapter
        if (out.terminal_adapter_alignments.size() < 2) {
            out.assign_type(kMissingAdapter);
        }
        // 3. search for polyAs that overlap the (extended) terminals
        assign_polya_to_terminals(read_seq, out);
        // 4. check if the entire read is spanned by a single polyA
        if (out.terminal_polya_alignments.size() == 2 &&
          ((out.terminal_polya_alignments[ReadTerminal::LEFT]->start == out.terminal_polya_alignments[ReadTerminal::RIGHT]->start))) {
	    //   (out.terminal_polya_alignments[ReadTerminal::LEFT]->end + 1 == out.terminal_polya_alignments[ReadTerminal::RIGHT]->start))) {
	        out.assign_type(kOnlyPolyA);
            out.cdna.valid = false;
        } else {
            out.cdna.valid = true;
        }
        // 5. check for sequencing adapter matches
        if (find_sequencing_adapters(read_seq, out)) {
            out.assign_type(kRetainedSmrtbell);
        }
        // 6. assemble the read structure from the terminal assignments
        assemble_structure(read_seq, out);
        // 7. check for additional PCR-based adapters internally
        if (out.cdna.valid && find_internal_adapters(read_seq, out)) {
            out.assign_type(kInternalAdapters);
        }
        // 8. match the assembled structure to a library of known artifacts
        if (!out.has_type()) {
            auto read_type = classify_structure(out.structure);
            if (read_type != kUnk) {
                out.assign_type(read_type);
            }
        }
        // 9. check if the structure can be classified as proper (ignoring additional polyAs)
        if (!out.has_type()) {
	        bool proper = can_reclassify_as_proper(out);
	        if (proper) {
                out.assign_type(kProper);
            } else {
                out.assign_type(kUnk);
            }
        }
	    // record classification for reporting
        stats_.update(read_seq, read_name, out);
    }

    void ReadClassifier::assign_adapters_to_terminals(const std::string& read_seq, ReadClassificationRecord& out) {
	    for (const auto& terminal : {ReadTerminal::LEFT, ReadTerminal::RIGHT}) {
	        std::pair<int, int> search_buffer = {0, cfg_.params.terminal_adapter_search_buffer};
            if (terminal != ReadTerminal::LEFT) {
                search_buffer = {read_seq.length() - cfg_.params.terminal_adapter_search_buffer,
                                 read_seq.length()};
            }
            std::unordered_map<seq_t, TargetAlignment*> target_alignments;
            for (const auto& target: cfg_.pcr_targets) {
		        target_alignments[target] = out.new_target_alignment();
                align_target(target, cfg_.target2seq.at(target),
                             read_seq, search_buffer.first, search_buffer.second,
                             cfg_.targets2lev.at(target), false, target_alignments[target]);
                // TODO: skip TSO alignment with split-sls (use the adapter result)
		        if (!is_valid_target_match(target, target_alignments)) continue;
		        //if (kSplitSLSSeqs.contains(target)) continue;
		        if (kTSOs.contains(target)) { // always assign the TSO if valid
                    out.terminal_adapter_alignments[terminal] = target_alignments.at(target);
                    const seq_t adapter = kTSOToAdapter.at(target);
                    if (cfg_.is_split_tso.at(target)) {
                        out.terminal_split_sls_alignments[terminal] = target_alignments.at(kAdapterToSplitSLS.at(adapter));
                    }
                    // remove the adapter alignment from valid matches
                    target_alignments[adapter]->valid = false;
                } else if (!out.terminal_adapter_alignments.contains(terminal)) {
                    out.terminal_adapter_alignments[terminal] = target_alignments.at(target);
                }
            }
        }
    }

    bool ReadClassifier::find_sequencing_adapters(const std::string& read_seq, ReadClassificationRecord& out) {
        bool found_matches = false;
        for (const auto target: cfg_.seq_targets) {
            auto* aln = out.new_target_alignment();
            align_target(target, cfg_.target2seq.at(target),
                         read_seq, std::nullopt, std::nullopt,cfg_.targets2lev.at(target),
                         !cfg_.params.compute_all_internal_matches, aln);
            if (aln->valid) {
                out.auxiliary_alignments.push_back(aln);
                if (!cfg_.params.compute_all_internal_matches) return true;
                found_matches = true;
            }
        }
        return found_matches;
    }

    bool ReadClassifier::find_internal_adapters(const std::string& read_seq, ReadClassificationRecord& out) {
        bool found_matches = false;
        for (const auto target: cfg_.pcr_targets) {
            if (kSplitSLSSeqs.contains(target)) continue;
            auto* aln = out.new_target_alignment();
            align_target(target, cfg_.target2seq.at(target),
                         read_seq, out.cdna.start, out.cdna.end + 1,
                         cfg_.targets2lev.at(target),
                         !cfg_.params.compute_all_internal_matches, aln);
            if (aln->valid) {
                out.auxiliary_alignments.push_back(aln);
		        if (!cfg_.params.compute_all_internal_matches) return true;
                found_matches = true;
            }
        }
        return found_matches;
    }

    void ReadClassifier::assign_polya_to_terminals(const std::string& read_seq, ReadClassificationRecord& out) {
        // search for all polyA and polyT runs in the read
        // allow for a small buffer overlap between the polya and the adapter matches
        const auto allowed_overlap = 5;
        std::vector<TargetAlignment> polya_runs;
        for (const auto target: cfg_.polya_targets) {
            find_polya(target, read_seq, cfg_.targets2lev.at(target),
                            cfg_.params.min_polyA_match, polya_runs);
        }
        if (polya_runs.empty()) return;
        // order runs by position
        std::sort(polya_runs.begin(), polya_runs.end(),
                  [](TargetAlignment &a, TargetAlignment &b) {return a.start < b.start;});
        // assign the left-most and right-most valid polyA runs
        // note: the same polyA run can be assigned at both terminals (for only-polyA reads)
        bool assigned_left = false;
        for (const auto& polya_run : polya_runs) {
            auto* aln = out.new_target_alignment();
            *aln = polya_run;
            if (polya_run.start < cfg_.params.terminal_polyA_search_buffer) {
                if (out.terminal_adapter_alignments.contains(ReadTerminal::LEFT) &&
                    polya_run.start < out.terminal_adapter_alignments[ReadTerminal::LEFT]->end - allowed_overlap) {
                    // polyA run starts before the adapter
                    out.auxiliary_alignments.push_back(aln);
                    out.has_external_polya = true;
                } else if (!assigned_left) {
                    out.terminal_polya_alignments[ReadTerminal::LEFT] = aln;
                    assigned_left = true;
                }
            }
            if (polya_run.end >= read_seq.length() - cfg_.params.terminal_polyA_search_buffer) {
                if (out.terminal_adapter_alignments.contains(ReadTerminal::RIGHT) &&
                    polya_run.end > out.terminal_adapter_alignments[ReadTerminal::RIGHT]->start + allowed_overlap) {
                    // polyA run starts after the adapter
                    out.auxiliary_alignments.push_back(aln);
                    out.has_external_polya = true;
                } else {
                    out.terminal_polya_alignments[ReadTerminal::RIGHT] = aln;
                }
            }
        }
    }

    // Returns true if this is a valid TSO match (by comparison to the PCR adapter only match)
    bool ReadClassifier::is_valid_tso_match(const seq_t& target, const std::unordered_map<seq_t,
                                            TargetAlignment*>& target_alignments) const {
        const seq_t adapter = kTSOToAdapter.at(target);
        if (cfg_.is_split_tso.at(target)) {
            // check if a valid match was found for the split SLS
            if(!target_alignments.at(kAdapterToSplitSLS.at(adapter))->valid) {
                target_alignments.at(target)->valid = false;
            }
            return target_alignments.at(target)->valid;
        }
        // the TSO is valid when the SLS part has at most the allowed number of errors
        // if the TSO matched but the adapter didn't match -- consider valid
        if (!target_alignments.at(adapter)->valid) return true;
        if (target_alignments.at(target)->score >
            target_alignments.at(adapter)->score + cfg_.targets2lev.at(kAdapterToSLS.at(adapter))) {
            target_alignments.at(target)->valid = false;
        }
        return target_alignments.at(target)->valid;
    }

    void ReadClassifier::assemble_structure(const std::string& read_seq, ReadClassificationRecord& out) {
        // TODO(viq): refine polyA boundaries using neighboring adapter matches (check for overlap)
        if (out.terminal_adapter_alignments.contains(ReadTerminal::LEFT)) {
            out.structure_alignments.push_back(out.terminal_adapter_alignments.at(ReadTerminal::LEFT));
        }
        if (out.terminal_polya_alignments.contains(ReadTerminal::LEFT)) {
            out.structure_alignments.push_back(out.terminal_polya_alignments.at(ReadTerminal::LEFT));
        }
        const int cdna_idx = out.structure_alignments.size();
        if (out.cdna.valid) {
            out.cdna.target_type = kCdna;
            out.cdna.start = 0;
            if (cdna_idx > 0) {
                out.cdna.start = out.structure_alignments[cdna_idx - 1]->end + 1;
                if (kTSOs.contains(out.structure_alignments[cdna_idx - 1]->target_type) &&
                    out.terminal_split_sls_alignments.contains(ReadTerminal::LEFT)) {
                    out.cdna.start = out.terminal_split_sls_alignments[ReadTerminal::LEFT]->end + 1;
                }
            }
            out.structure_alignments.push_back(&out.cdna);
        }
        if (out.cdna.valid && out.terminal_polya_alignments.contains(ReadTerminal::RIGHT)) {
            out.structure_alignments.push_back(out.terminal_polya_alignments.at(ReadTerminal::RIGHT));
        }
        if (out.terminal_adapter_alignments.contains(ReadTerminal::RIGHT)) {
            out.structure_alignments.push_back(out.terminal_adapter_alignments.at(ReadTerminal::RIGHT));
        }
        out.cdna.end = read_seq.length() - 1;
        if (cdna_idx < out.structure_alignments.size() - 1) {
            out.cdna.end = out.structure_alignments[cdna_idx + 1]->start - 1;
            if (kTSOs.contains(out.structure_alignments[cdna_idx + 1]->target_type) &&
                out.terminal_split_sls_alignments.contains(ReadTerminal::RIGHT)) {
                out.cdna.end = out.terminal_split_sls_alignments[ReadTerminal::RIGHT]->start - 1;
            }
        }

	if (out.cdna.end < out.cdna.start) {
		out.cdna.valid = false;
		out.structure_alignments.erase(out.structure_alignments.begin() + cdna_idx);
        }
        // assemble the read structure (order of targets)
        for (const auto aln : out.structure_alignments) {
            out.structure.push_back(aln->target_type);
        }
    }

    class_t ReadClassifier::classify_structure(std::vector<seq_t>& structure) const {
        class_t read_class{kUnk};
        for (const seq_t tso_adapter: cfg_.params.tso_adapters) {
            for (const seq_t rt_adapter: cfg_.params.rt_adapters) {
                read_class = lookup_structure(structure, tso_adapter, rt_adapter, false);
                if (read_class != kUnk) return read_class;
            }
        }
        return read_class;
    }

    bool ReadClassifier::can_reclassify_as_proper(ReadClassificationRecord& out) {
        if (out.terminal_adapter_alignments.size() != 2 || out.terminal_polya_alignments.size() < 2) return false;
        for (const seq_t tso_adapter: cfg_.params.tso_adapters) {
            for (const seq_t rt_adapter: cfg_.params.rt_adapters) {
                for (auto i = 0; i < out.structure.size(); i++) {
                    if (!kPolyTails.contains(out.structure[i])) continue;
                    auto alternative = out.structure;
                    alternative.erase(alternative.begin() + i);
                    class_t read_class = lookup_structure(alternative, tso_adapter, rt_adapter, true);
                    if (read_class == kProper) {
                        out.auxiliary_alignments.push_back(out.structure_alignments[i]);
                        // remove the extra polyA from reported structure and fix the cDNA coordinates
                        if (out.terminal_polya_alignments[ReadTerminal::LEFT]->start == out.structure_alignments[i]->start) {
                            out.cdna.start = out.structure_alignments[i-1]->end + 1;
                        } else {
                            out.cdna.end = out.structure_alignments[i+1]->start - 1;
                        }
                        out.structure.erase(out.structure.begin() + i);
                        out.structure_alignments.erase(out.structure_alignments.begin() + i);
                        out.has_internal_polya = true;
                        return true;
                    }
                }
            }
        }
        return false;
    }

    class_t ReadClassifier::lookup_structure(std::vector<seq_t>& structure, seq_t tso_adapter, seq_t rt_adapter,
                                             bool check_proper_only) const {
        for (const auto& [read_type, structure_map]: cfg_.class2structure) {
            if (check_proper_only && read_type != kProper) continue;
            for (const auto& structure_for_type : structure_map.at({tso_adapter, rt_adapter})) {
		if (structure.size() != structure_for_type.size()) continue;
                // forward match
                if (std::equal(structure.begin(), structure.end(),
                               structure_for_type.begin())) { return read_type; }
                // reverse complement match
                if (std::equal(structure.begin(), structure.end(),
                        /* traverse in reverse */ structure_for_type.rbegin(),
                        [&](const auto& s1, const auto& s2)
                        { return s1 == kSeqToRC.at(s2);})) { return read_type; }
            }
	}
        return kUnk;
    }

    ///////////////////////
    ////  Annotation   ////
    ///////////////////////
    void ReadClassifier::add_tags(io::read_t* read, ReadClassificationRecord& result) {
        // "pr" indicates if the read is proper
        int32_t is_proper_tag = result.is_proper();
        io::add_tag_value(read, "pr", 'i', sizeof(int32_t), reinterpret_cast<uint8_t*>(&is_proper_tag));

        // for proper reads only
        if (result.is_proper()) {
            // "cd" stores the start and end coordinates of the cDNA segment (for proper reads only)
            std::vector<int32_t> cdna_coords = {result.cdna.start, result.cdna.end};
            io::add_tag_array(read, "cd", 'I', cdna_coords.size(), reinterpret_cast<uint8_t*>(cdna_coords.data()));

            // "sf" indicates that the read strand is flipped
            int32_t should_flip = result.is_proper_and_rc();
            io::add_tag_value(read, "sf", 'i', sizeof(int32_t), reinterpret_cast<uint8_t*>(&should_flip));
        }

	    // "lb" stores the comma-separated list of all the assigned read artifact classes
        std::string labels = utils::implode(result.classes);
        io::add_tag_value(read, "lb", 'Z', labels.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(labels.c_str())));

        // "st" stores the structural representation label
        io::add_tag_value(read, "st", 'Z', result.structure_label.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(result.structure_label.c_str())));

        if (result.cannot_be_classified()) return;

        // "th" stores targets with valid alignments in the read
        io::add_tag_value(read, "th", 'Z', result.valid_targets_label.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(result.valid_targets_label.c_str())));

        // "ch" stores the structural profile (ordered sequence of sequence types)
        const auto& coords_str = result.get_structure_coords();
        io::add_tag_value(read, "ch", 'Z', coords_str.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(coords_str.c_str())));

        // "al" stores the sorted coordinates of all valid target alignments (if configured)
        if (cfg_.params.compute_all_internal_matches) {
            std::sort(result.alignments.begin(), result.alignments.end(),
                      [](auto const& a, auto const& b)
                      { return a->start < b->start;});
            const auto& all_coords_str = result.get_all_valid_coords();
            io::add_tag_value(read, "al", 'Z', all_coords_str.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(all_coords_str.c_str())));
        }
    }

    std::string ReadClassificationRecord::summarize() {
	    std::ostringstream s;
        s << "\nClass: " << utils::implode(classes) << "\n";
        s << "Structure: " << structure_label << "\n";
        s << "Matches: " << valid_targets_label << "\n\n";
        for (auto i = 0; i < structure.size(); i++) {
            const auto aln = structure_alignments[i];
            s << ">>>>>>>> " << structure[i] << ": " << aln->start << "-" << aln->end;
	        if (structure[i] != kCdna) {
                s << " (lev: " << aln->score << ")\n";
                aln->display_alignment(s, true);
            } else {
                s << "\n";
            }
            if (kTSOs.contains(structure[i])) {
                const auto terminal = i == 0 ? ReadTerminal::LEFT: ReadTerminal::RIGHT;
                if (terminal_split_sls_alignments.contains(terminal)) {
                    const auto& sls_aln = terminal_split_sls_alignments.at(terminal);
                    s << ">>>>>>>> " << sls_aln->target_type << ": " << sls_aln->start << "-" << sls_aln->end;
                    s << " (lev: " << sls_aln->score << ")\n";
		            sls_aln->display_alignment(s, true);
                }
            }
        }
        s << "Auxiliary matches: " << (auxiliary_alignments.empty() ? "(none)" : "") << "\n";
	    std::sort(auxiliary_alignments.begin(), auxiliary_alignments.end(),
                  [](TargetAlignment* a, TargetAlignment* b) {return a->start < b->start;});
	    for (const auto& aln : auxiliary_alignments) {
            const std::string target_seq = std::string(aln->target_seq);
            const std::string read_seq = std::string(aln->read_seq);
	        if (aln->score_only) {
                TargetAlignment new_aln;
		        align_target(aln->target_type, target_seq, read_seq, cdna.start, cdna.end + 1,
                             aln->score, false, &new_aln);
                s << ">>>>>>>> " << new_aln.target_type << ": " << new_aln.start << "-" << new_aln.end;
                s << " (lev: " << new_aln.score << ")\n";
                new_aln.display_alignment(s, true);
            } else {
                s << ">>>>>>>> " << aln->target_type << ": " << aln->start << "-" << aln->end;
                s << " (lev: " << aln->score << ")\n";
                aln->display_alignment(s, true);
            }
	    }
        return s.str();
    }

    ///////////////////////
    ////   Reporting   ////
    ///////////////////////

    inline void ReadStats::update(const std::string& read_seq, const std::string& read_name,
                                  ReadClassificationRecord& annotation) {
        annotation.set_class_label();
        annotation.set_structure_label();
        annotation.set_valid_targets_label();
        class_counts[annotation.class_label]++;
	    const auto &read_label = std::make_pair(annotation.class_label, annotation.structure_label);
        structure_counts[read_label]++;
	    if (structure_counts[read_label] > n_examples_to_record) return;
	    std::ostringstream s;
        s << "----------------------------\n";
        s << "Name: " << read_name << "\n";
        s << "Length: " << read_seq.length() << "\n";
        s << "Sequence: " << read_seq << "\n";
        s << annotation.summarize() << "\n";
        class_examples[read_label].push_back(s.str());
    }

    void ReadStats::generate_reports(const std::filesystem::path& report_dir) const {
        std::ofstream report_classes(report_dir / "class_counts.tsv");
        report_classes << "read_class\tcount\n";
        for (const auto &[read_type, count]: class_counts) {
            report_classes << read_type << "\t" << count << "\n";
        }
        std::ofstream report_structures(report_dir / "structure_counts.tsv");
        report_structures << "read_type\tstructure\tcount\n";
        for (const auto &[read_info, count]: structure_counts) {
            const auto &[read_type, read_structure] = read_info;
            report_structures << read_type << "\t" << read_structure << "\t" << count << "\n";
        }
        generate_read_parsing_logs(report_dir);
    }

    void ReadStats::generate_read_parsing_logs(const std::filesystem::path& report_dir) const {
        for (const auto& read_class: std::views::keys(class_counts)) {
            std::ofstream report_read_class(report_dir / (read_class + ".txt"));
            report_read_class << "==========================\n";
            report_read_class << read_class << " reads\n";
            report_read_class << "==========================\n";
            for (const auto& read_label : std::views::keys(class_examples)) {
                if (read_label.first != read_class) continue;
                report_read_class << utils::implode(class_examples.at(read_label), "\n");
            }
        }
    }
}
