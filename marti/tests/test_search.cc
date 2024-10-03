/**
 * @brief Alignment and polyA search tests.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <array>
#include <string>
#include <vector>

#include "marti/search.h"

using namespace marti;
TEST_CASE("target_search") {

    std::string target_type = "test_adapter";
    std::string read_seq = "AAAAATTTAAACCCGGGTTTTT";

    SUBCASE("exact_match") {
        std::string target_seq = "TTTAAACCCGGG";
        TargetAlignment aln;
        align_target(target_type, target_seq, read_seq, std::nullopt, std::nullopt, 0, false, &aln);
        CHECK_EQ(aln.target_type, target_type);
        CHECK_EQ(aln.target_seq, target_seq);
        CHECK_EQ(aln.read_seq, read_seq);
        CHECK(aln.valid);
        CHECK_FALSE(aln.score_only);
        CHECK_EQ(aln.score, 0);
        CHECK_EQ(aln.start, 5);
        CHECK_EQ(aln.end, 16);
    }
    SUBCASE("within_max_lev_dist") {
        std::string target_seq = "TTTAAATCCGGG";
        TargetAlignment aln;
        align_target(target_type, target_seq, read_seq, std::nullopt, std::nullopt,
                     1, false, &aln);
        CHECK(aln.valid);
        CHECK_FALSE(aln.score_only);
        CHECK_EQ(aln.score, 1);
        CHECK_EQ(aln.start, 5);
        CHECK_EQ(aln.end, 16);
    }

    SUBCASE("too_many_errors") {
        std::string target_seq = "TTTAAATTCGGG";
        TargetAlignment aln;
        align_target(target_type, target_seq, read_seq, std::nullopt, std::nullopt,
                     1, false, &aln);
        CHECK_FALSE(aln.valid);
        CHECK_FALSE(aln.score_only);
        CHECK_EQ(aln.score, 0);
        CHECK_EQ(aln.start, 0);
        CHECK_EQ(aln.end, 0);
    }

    SUBCASE("with_offset") {
        std::string target_seq = "TTTAAACCCGGG";
        TargetAlignment aln;
        align_target(target_type, target_seq, read_seq, 4, 17,0, false, &aln);
        CHECK(aln.valid);
        CHECK_FALSE(aln.score_only);
        CHECK_EQ(aln.score, 0);
        CHECK_EQ(aln.start, 5);
        CHECK_EQ(aln.end, 16);
    }

    SUBCASE("with_too_much_offset") {
        std::string target_seq = "TTTAAACCCGGG";
        TargetAlignment aln;
        align_target(target_type, target_seq, read_seq, 10, 17,0, true, &aln);
        CHECK_FALSE(aln.valid);
        CHECK_EQ(aln.score, 0);
        CHECK_EQ(aln.start, 0);
        CHECK_EQ(aln.end, 0);
    }
}

TEST_CASE("find_polya") {
    SUBCASE("polya_single_exact_match") {
        std::string read_seq = "GGGTTTTTTTTGGGAAAAAAAAGGG";
        std::vector<TargetAlignment> polya_runs;
        find_polya(kPolyA, read_seq, 0, 3, polya_runs);
        CHECK_EQ(polya_runs.size(), 1);
        CHECK(polya_runs[0].valid);
        CHECK_EQ(polya_runs[0].start, 14);
        CHECK_EQ(polya_runs[0].end, 21);

        std::vector<TargetAlignment> polyt_runs;
        find_polya(kPolyT, read_seq, 0, 3, polyt_runs);
        CHECK_EQ(polyt_runs.size(), 1);
        CHECK(polyt_runs[0].valid);
        CHECK_EQ(polyt_runs[0].start, 3);
        CHECK_EQ(polyt_runs[0].end, 10);
    }

    SUBCASE("polya_single_match_with_error") {
        std::string read_seq = "GGGTTTTTTTTGGGAGAAAAAAGGG";
        std::vector<TargetAlignment> polya_runs;
        find_polya(kPolyA, read_seq, 1, 8, polya_runs);
        CHECK_EQ(polya_runs.size(), 1);
        CHECK(polya_runs[0].valid);
        CHECK_EQ(polya_runs[0].score, 1);
        CHECK_EQ(polya_runs[0].start, 14);
        CHECK_EQ(polya_runs[0].end, 21);
    }

    SUBCASE("polya_too_much_error") {
        std::string read_seq = "GGGTTTTTTTTGGGAAATAAAAGGG";
        std::vector<TargetAlignment> polya_runs;
        find_polya(kPolyA, read_seq, 0, 8, polya_runs);
        CHECK_EQ(polya_runs.size(), 0);
    }

    SUBCASE("polya_multiple_match") {
        std::string read_seq = "GGGTTTTTTTTGGGAAAAGAAAGGG";
        std::vector<TargetAlignment> polya_runs;
        find_polya(kPolyA, read_seq, 0, 3, polya_runs);
        CHECK_EQ(polya_runs.size(), 2);
        CHECK(polya_runs[0].valid);
        CHECK_EQ(polya_runs[0].start, 14);
        CHECK_EQ(polya_runs[0].end, 17);
        CHECK(polya_runs[1].valid);
        CHECK_EQ(polya_runs[1].start, 19);
        CHECK_EQ(polya_runs[1].end, 21);
    }
}
