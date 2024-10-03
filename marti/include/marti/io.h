/**
 * Lightweight IO wrappers for:
 * -- HtsLib-based BAM reading and writing functionality
 * -- YAML-CPP YAML parsing functionality
 * and basic logging
 */
#ifndef MARTI_INCLUDE_MARTI_IO_H_
#define MARTI_INCLUDE_MARTI_IO_H_

#include <filesystem>
#include <htslib/sam.h>
#include <iostream>
#include <optional>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "marti/definitions.h"
#include "marti/utils.h"

namespace marti::io {
    struct ReadDeleter {
        void operator()(bam1_t *r) { bam_destroy1(r); }
    };

    struct HeaderDeleter {
        void operator()(bam_hdr_t *h) { bam_hdr_destroy(h); }
    };

    using read_ptr_t = std::unique_ptr<bam1_t, ReadDeleter>;
    using read_buffer_t = std::vector<read_ptr_t>;
    using bam_header_t = std::unique_ptr<bam_hdr_t, HeaderDeleter>;
    using read_t = bam1_t;

    class BamFile {
    public:
        explicit BamFile(const std::string &bam_file_name, const std::string &mode) {
            if (bam_file_name.empty()) {
                throw std::invalid_argument("[IO] BAM file was not provided (empty).");
            }
            fp_ = sam_open(bam_file_name.c_str(), mode.c_str());
            if(!fp_) {
                throw std::runtime_error("[IO] Failed to open BAM file: " + bam_file_name);
            }
        }

        ~BamFile() {
            sam_close(fp_);
        }

    protected:
        samFile *fp_;
    };

    class BamWriter : public BamFile {
    public:
        BamWriter(const std::string &bam_file_name, const bam_header_t &header) :
                BamFile(bam_file_name, "wb"), header_(header) {
            if (sam_hdr_write(fp_, header_.get()) < 0) {
                throw std::runtime_error("[IO] Failed to write BAM header to output file.");
            }
        }

        void write(const read_t *read) {
            if (sam_write1(fp_, header_.get(), read) < 0) {
                throw std::runtime_error("[IO] Failed to write read to output BAM file.");
            }
        }

    private:
        const bam_header_t &header_;
    };

    class BamReader : public BamFile {
    public:
        explicit BamReader(const std::string &bam_file_name,
                           const std::optional<int> batch_size = std::nullopt,
                           const std::optional<int> max_reads = std::nullopt) :
                BamFile(bam_file_name, "r"), batch_size_(batch_size), n_reads_processed_(0),
                max_reads_to_process_(max_reads) {
            header_ = bam_header_t(sam_hdr_read(fp_));
            if (!header_) {
                throw std::runtime_error("[IO] Failed to read BAM header.");
            }
            if (batch_size_.has_value()) {
                init_buffer();
            }
        }

        [[nodiscard]] inline bam_header_t const &get_header() const { return header_; }
        [[nodiscard]] inline size_t get_n_processed() const { return n_reads_processed_; }
        [[nodiscard]] inline bool limit_reached() const {
            return max_reads_to_process_.has_value() && n_reads_processed_ >= max_reads_to_process_;
        }

        read_buffer_t &read_next_batch() {
            if (!batch_size_.has_value()) {
                throw std::runtime_error("[IO] Batched-mode reading was not initialized.");
            }
            const auto processed_at_start = n_reads_processed_;
            for (auto i = 0; i < batch_size_.value(); i++) {
                if (limit_reached()) break; // reached the specified number of records
                if (sam_read1(fp_, header_.get(), buffer_[i].get()) < 0) break;
                n_reads_processed_++;
            }
            const auto n_processed = n_reads_processed_ - processed_at_start;
            if (n_processed < batch_size_.value()) {
                buffer_.resize(n_processed);
            }
            return buffer_;
        }

        std::optional<read_t *> read() {
            read_ptr_t r(bam_init1());
            const int ret = sam_read1(fp_, header_.get(), r.get());
            if (ret < 0) return std::nullopt;
            return r.get();
        }

    private:
        void init_buffer() {
            buffer_.resize(batch_size_.value());
            for (int i = 0; i < batch_size_; i++) {
                buffer_[i] = read_ptr_t(bam_init1());
            }
        }

        bam_header_t header_;
        read_buffer_t buffer_;
        size_t n_reads_processed_;
        std::optional<size_t> batch_size_;
        std::optional<size_t> max_reads_to_process_;
    };

    inline void add_tag_array(read_t *read, const std::string &tag_name, char array_data_type,
                              const int array_len, uint8_t *data) {
        if (bam_aux_update_array(read, tag_name.c_str(), array_data_type, array_len, data)) {
            throw std::runtime_error("[IO] Failed to add array tag to read.");
        }
    }

    inline void add_tag_value(read_t *read, const std::string &tag_name, char value_type,
                              const int value_size, uint8_t *data) {
        if (bam_aux_append(read, tag_name.c_str(), value_type, value_size, data)) {
            throw std::runtime_error("[IO] Failed to add tag to read.");
        }
    }

    inline bool has_tag(const read_t *read, const std::string &tag_name) {
        return bam_aux_get(read, tag_name.c_str()) != nullptr;
    }

    inline uint8_t *get_tag_data(const read_t *read, const std::string &tag_name) {
        uint8_t *tag_data = bam_aux_get(read, tag_name.c_str());
        if (!tag_data) throw std::runtime_error("Missing BAM tag: " + tag_name);
        return tag_data;
    }

    template<typename T>
    inline T get_tag(const read_t *read, const std::string &tag_name);

    template<typename T>
    inline T get_tag_array(const read_t *read, const std::string &tag_name, int array_idx);

    template<>
    inline double get_tag<double>(const read_t *read, const std::string &tag_name) {
        return bam_aux2f(get_tag_data(read, tag_name));
    }

    template<>
    inline int get_tag<int>(const read_t *read, const std::string &tag_name) {
        return bam_aux2i(get_tag_data(read, tag_name));
    }

    template<>
    inline int get_tag_array<int>(const read_t *read, const std::string &tag_name, const int array_idx) {
        return bam_auxB2i(get_tag_data(read, tag_name), array_idx);
    }

    inline std::string get_seq(const read_t *read) {
        uint8_t *seq = bam_get_seq(read);
        std::string sequence;
        for (int i = 0; i < read->core.l_qseq; ++i) {
            sequence += seq_nt16_str[bam_seqi(seq, i)];
        }
        return sequence;
    }

    inline std::string get_name(const read_t *read) {
        return bam_get_qname(read);
    }

    inline void trim_and_orient_read(read_t *orig_read, const int start, const int end, const bool strand_switch, read_t *trimmed_read) {
        std::string seq = get_seq(orig_read).substr(start, end);
        if (strand_switch) {
            seq = utils::reverse_complement(seq);
        }
        bam_set1(trimmed_read, get_name(orig_read).size(), bam_get_qname(orig_read),
                 orig_read->core.flag, orig_read->core.tid, orig_read->core.pos, orig_read->core.qual,
                 0, NULL, orig_read->core.mtid, orig_read->core.mpos, orig_read->core.isize,
                 seq.size(), seq.c_str(), reinterpret_cast<char *>(bam_get_qual(orig_read)) + start,
                 bam_get_l_aux(orig_read));
        memcpy(bam_get_aux(trimmed_read), bam_get_aux(orig_read), bam_get_l_aux(orig_read));
        trimmed_read->l_data += bam_get_l_aux(orig_read);
    }

    class YAMLFile {
    public:
        explicit YAMLFile(const std::string& config_path) {
            config_ = YAML::LoadFile(config_path);
            if (!config_.IsDefined() || !config_.IsMap()) {
                throw std::runtime_error("[IO] Error reading config file: " + config_path);
            }
        }
        template <typename T>
        inline void load(const std::string_view key, T& var, const T& default_val) const {
            if (!config_[key]) {
                var = default_val;
            } else {
                var = config_[key].as<T>();
            }
        }
        template <typename T>
        inline void load(const std::string_view& key, std::set<T>& var) const {
            if (!config_[key]) return;
            for (const auto& el : config_[key]) {
                var.insert(el.as<T>());
            }
        }
    private:
        YAML::Node config_;
    };

#define MARTI_LOG(log_file, msg) do { \
    log_file << msg << std::endl; \
    std::cout << msg << std::endl; \
} while (0)

}

#endif //MARTI_INCLUDE_MARTI_IO_H_
