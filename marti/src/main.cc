#include <filesystem>
#include <iostream>
#include <omp.h>
#include <string>

#include "marti/definitions.h"
#include "marti/classifier.h"
#include "marti/config.h"
#include "marti/io.h"

void run_marti(marti::MartiConfig& cfg) {
    using namespace marti;
    double io_time = 0;

    // (debug mode) test a specific sequence
    if (!cfg.params.test_seq.empty()) {
        auto classifier = ReadClassifier(cfg);
        ReadClassificationRecord out;
        classifier.classify(cfg.params.test_seq, out);
        std::cout << out.summarize();
        exit(0);
    }

    // setup BAM IO
    io::BamReader input_bam(cfg.params.input_bam, cfg.params.n_max_reads_in_mem, cfg.params.max_reads_to_process);
    io::BamWriter output_bam(cfg.output_bam, input_bam.get_header());
    std::optional<io::BamWriter> output_proper_bam;
    if (cfg.params.output_proper_trimmed) {
        output_proper_bam.emplace(cfg.output_proper_bam, input_bam.get_header());
    }
    // setup classifiers (one per thread)
    std::vector<ReadClassifier> classifiers(cfg.params.n_threads, ReadClassifier(cfg));
    double timer_start;
    while (true) {
        timer_start = omp_get_wtime();
        // load a batch of reads
        auto &batch = input_bam.read_next_batch();
        io_time += (omp_get_wtime() - timer_start);
        if (batch.empty()) break;
        // split the batch across multiple threads for parallel processing
        const auto reads_per_thread = static_cast<int>(std::ceil(
                static_cast<double>(batch.size()) / cfg.params.n_threads));
        #pragma omp parallel for default(none) shared(batch, reads_per_thread, classifiers)
        for (size_t worker_start = 0; worker_start < batch.size(); worker_start += reads_per_thread) {
            // each thread classifies a subset of reads from the batch
            const auto worker_end = std::min(worker_start + reads_per_thread, batch.size());
            const auto tid = omp_get_thread_num();
            for (auto read_id = worker_start; read_id < worker_end; read_id++) {
                // classify each read, add new tags with marti results
                classifiers.at(tid).classify_and_annotate(batch[read_id].get());
            }
        }
        // write the annotated read batch
        timer_start = omp_get_wtime();
        for (const auto &read: batch) {
            output_bam.write(read.get());
            // if configured, output a separate BAM file with trimmed proper reads only
            if (output_proper_bam && io::get_tag<int>(read.get(), "pr")) {
                const int cdna_start = io::get_tag_array<int>(read.get(), "cd", 0);
                const int cdna_end = io::get_tag_array<int>(read.get(), "cd", 1);
                const bool switch_strand = io::get_tag_array<int>(read.get(), "sf", 1);
                io::read_ptr_t trimmed_read(bam_init1());
                io::trim_and_orient_read(read.get(), cdna_start, cdna_end, switch_strand,trimmed_read.get());
                output_proper_bam->write(trimmed_read.get());
            }
        }
        io_time += (omp_get_wtime() - timer_start);
        MARTI_LOG(cfg.log_file, "Processed " << input_bam.get_n_processed() << " reads");
    }
    MARTI_LOG(cfg.log_file, "Total number of reads processed: " << input_bam.get_n_processed());
    MARTI_LOG(cfg.log_file, "BAM IO time: " << io_time << " seconds");

    MARTI_LOG(cfg.log_file, "Generating reports...");
    // combine per-thread classification results
    timer_start = omp_get_wtime();
    ReadStats stats;
    for (const auto &c: classifiers) {
        stats.merge(c.get_read_stats());
    }
    // generate classification reports
    stats.generate_reports(cfg.report_dir);
    MARTI_LOG(cfg.log_file, "Reporting time: " << omp_get_wtime() - timer_start << " seconds");
}

int main(int argc, char const ** argv) {
    std::cout << "**********************************************\n";
    std::cout << "* marti (" << MARTI_VERSION << "): cDNA artifact classification *\n";
    std::cout << "**********************************************\n";
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <yaml_config>\n";
        return EXIT_FAILURE;
    }
    const std::string config_fname = argv[1];
    if (!std::filesystem::exists(config_fname)) {
        std::cerr << "YAML config file " << config_fname << " does not exist\n";
        return EXIT_FAILURE;
    }

    const double start_time = omp_get_wtime();
    try {
        marti::MartiConfig cfg(config_fname);
        run_marti(cfg);
        MARTI_LOG(cfg.log_file, "Total time: " << omp_get_wtime() - start_time << " seconds");
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
