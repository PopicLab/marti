/**
 * @brief Marti configuration (parsing and verification) tests
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <string>

#include "marti/definitions.h"
#include "marti/config.h"

using namespace marti;
TEST_CASE("marti_config") {
    MartiConfig::MartiParams marti_params({
          .terminal_adapter_search_buffer=100,
          .terminal_polyA_search_buffer=150,
          .min_polyA_match=20,
          .max_err_rate=0,
          .max_err_rate_polya=0,
    });
    marti_params.target2seq = {{kMars, "ACTGACTGGT"},
                               {kVenus, "TCGATCGA"},
                               {kSmrtbell, "ATATATA"},
                               {kMarsSls, ""},
                               {kMarsSplitSls, ""},
                               {kVenusSls, ""},
                               {kVenusSplitSls, ""}};
    SUBCASE("valid_config") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        CHECK_NOTHROW(MartiConfig(marti_params, ""));
    }
    SUBCASE("missing_tso_adapters") {
        marti_params.tso_adapters = {};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("unexpected_adapter") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {"unexpected"};
        marti_params.target2seq[kMarsSls] = "ACTG";
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("missing_adapter_seq") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        marti_params.target2seq[kMars] = {};
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("missing_tso_sls") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = {};
        marti_params.target2seq[kMarsSplitSls] = {};
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("ambiguous_tso_sls") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        marti_params.target2seq[kMarsSplitSls] = "GTCTG";
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("adapter_buffer_too_small") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        marti_params.terminal_adapter_search_buffer = 2;
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("polya_buffer_too_small") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        marti_params.terminal_adapter_search_buffer = 20;
        marti_params.terminal_polyA_search_buffer = 10;
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("negative_buffer") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        marti_params.terminal_adapter_search_buffer = -20;
        CHECK_THROWS_AS(MartiConfig(marti_params, ""), std::invalid_argument);
    }
    SUBCASE("set_targets_seq") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        MartiConfig marti_config(marti_params, "");
        CHECK((marti_config.target2seq.contains(kMars)));
        CHECK((marti_config.target2seq.contains(kMarsTso)));
        CHECK((marti_config.target2seq.contains(kMarsTsoRC)));
        CHECK((marti_config.target2seq.contains(kVenus)));
        CHECK((marti_config.target2seq.contains(kVenusRC)));
        CHECK((!marti_config.target2seq.contains(kVenusTso)));
    }
    SUBCASE("set_targets") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSplitSls] = "ACTG";
        MartiConfig marti_config(marti_params, "");
        CHECK((utils::contains(marti_config.pcr_targets, kMars)));
        CHECK((utils::contains(marti_config.pcr_targets, kMarsRC)));
        CHECK((utils::contains(marti_config.pcr_targets, kMarsTso)));
        CHECK((utils::contains(marti_config.pcr_targets, kMarsTsoRC)));
        CHECK((utils::contains(marti_config.pcr_targets, kMarsSplitSls)));
        CHECK((utils::contains(marti_config.pcr_targets, kMarsSplitSlsRC)));
        CHECK((utils::contains(marti_config.pcr_targets, kVenus)));
        CHECK((utils::contains(marti_config.pcr_targets, kVenusRC)));
        CHECK((!utils::contains(marti_config.pcr_targets, kVenusTso)));
        CHECK((!utils::contains(marti_config.pcr_targets, kVenusTsoRC)));
        CHECK((!utils::contains(marti_config.pcr_targets, kMarsSls)));
        CHECK((!utils::contains(marti_config.pcr_targets, kMarsSlsRC)));
        CHECK((utils::contains(marti_config.seq_targets, kSmrtbell)));
        CHECK((utils::contains(marti_config.seq_targets, kSmrtbellRC)));
    }

    SUBCASE("max_lev") {
        marti_params.tso_adapters = {std::string(kMars)};
        marti_params.rt_adapters = {std::string(kVenus)};
        marti_params.target2seq[kMarsSls] = "ACTG";
        MartiConfig marti_config(marti_params, "");
        CHECK((marti_config.targets2lev[kMars] == 0));
        CHECK((marti_config.targets2lev[kMarsTso] == 0));
        CHECK((marti_config.targets2lev[kVenus] == 0));
        marti_params.max_err_rate = 0.1;
        MartiConfig marti_config2(marti_params, "");
        CHECK((marti_config2.targets2lev[kMars] == 1));
    }
}