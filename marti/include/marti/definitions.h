#ifndef MARTI_INCLUDE_MARTI_DEFINITIONS_H_
#define MARTI_INCLUDE_MARTI_DEFINITIONS_H_

#include <set>
#include <string_view>
#include <unordered_map>

namespace marti {

#define MARTI_VERSION "v0.1"

// Main Marti types
using seq_t = std::string_view;
using class_t = std::string_view;

// Oligo types definitions
#define INIT_SEQ_TYPES(seq_type, seq_name) \
    inline constexpr seq_t seq_type(seq_name); \
    inline constexpr seq_t seq_type##RC(seq_name "_RC"); \
    inline constexpr seq_t seq_type##Tso(seq_name "_TSO");\
    inline constexpr seq_t seq_type##TsoRC(seq_name "_TSO_RC"); \
    inline constexpr seq_t seq_type##Sls(seq_name "_SLS");\
    inline constexpr seq_t seq_type##SlsRC(seq_name "_SLS_RC");\
    inline constexpr seq_t seq_type##SplitSls(seq_name "_SPLIT_SLS"); \
    inline constexpr seq_t seq_type##SplitSlsRC(seq_name "_SPLIT_SLS_RC");

INIT_SEQ_TYPES(kMars, "adapterA");
INIT_SEQ_TYPES(kVenus, "adapterB");
inline constexpr seq_t kCdna = "cDNA";
inline constexpr seq_t kPolyA = "polyA";
inline constexpr seq_t kPolyT = "polyT";
inline constexpr seq_t kSmrtbell = "SMRTbell";
inline constexpr seq_t kSmrtbellRC = "SMRTbell_RC";

// Lookup

// sets of oligos grouped by type
const std::set<seq_t> kTSOs = {kMarsTso, kMarsTsoRC, kVenusTso, kVenusTsoRC};
const std::set<seq_t> kSeqAdapters = {kSmrtbell, kSmrtbellRC};
const std::set<seq_t> kPolyTails = {kPolyA, kPolyT};
const std::set<seq_t> kSplitSLSSeqs = {kMarsSplitSls, kMarsSplitSlsRC, kVenusSplitSls, kVenusSplitSlsRC};
const std::set<seq_t> kRCTSOs = {kMarsTsoRC, kVenusTsoRC};

// map from a TSO to its PCR adapter
const std::unordered_map<seq_t, seq_t> kTSOToAdapter = {
        {kMarsTso,    kMars},
        {kVenusTso,   kVenus},
        {kMarsTsoRC,  kMarsRC},
        {kVenusTsoRC, kVenusRC}
};

// map from a PCR adapter to its corresponding TSO
const std::unordered_map<seq_t, seq_t> kAdapterToTSO = {{kMars,  kMarsTso}, {kVenus, kVenusTso}};

// map from a PCR adapter to the SLS
const std::unordered_map<seq_t, seq_t> kAdapterToSLS = {
        {kMars,    kMarsSls},
        {kMarsRC,  kMarsSlsRC},
        {kVenus,   kVenusSls},
        {kVenusRC, kVenusSlsRC}
    };

// map from a PCR adapter to the split SLS
const std::unordered_map<seq_t, seq_t> kAdapterToSplitSLS = {
        {kMars,    kMarsSplitSls},
        {kMarsRC,  kMarsSplitSlsRC},
        {kVenus,   kVenusSplitSls},
        {kVenusRC, kVenusSplitSlsRC}
};

#define POPULATE_RC_ENTRIES(seq_type) \
    {seq_type,         seq_type##RC},  \
    {seq_type##RC,     seq_type},     \
    {seq_type##Tso,    seq_type##TsoRC},\
    {seq_type##TsoRC,  seq_type##Tso},\
    {seq_type##Sls,    seq_type##SlsRC},\
    {seq_type##SlsRC,       seq_type##Sls},\
    {seq_type##SplitSls,    seq_type##SplitSlsRC},\
    {seq_type##SplitSlsRC,  seq_type##SplitSls},\

// map from a sequence to its reverse complement
const std::unordered_map<seq_t, seq_t> kSeqToRC = {
        POPULATE_RC_ENTRIES(kMars)
        POPULATE_RC_ENTRIES(kVenus)
        {kPolyA,           kPolyT},
        {kPolyT,           kPolyA},
        {kSmrtbell,        kSmrtbellRC},
        {kSmrtbellRC,      kSmrtbell},
        {kCdna, kCdna},
};

// Artifact classes
inline constexpr class_t kProper = "Proper";
inline constexpr class_t kTsoTso = "TsoTso";
inline constexpr class_t kRtRt = "RtRt";
inline constexpr class_t kNoPolyA = "InternalPrimingRT";
inline constexpr class_t kNoTSO = "InternalPrimingTSO";
inline constexpr class_t kOnlyPolyA = "OnlyPolyA";
inline constexpr class_t kMissingAdapter = "MissingAdapter";
inline constexpr class_t kInternalAdapters = "InternalAdapter";
inline constexpr class_t kRetainedSmrtbell = "RetainedSMRTBell";
inline constexpr class_t kExtraSLS = "ExtraSls"; // TODO
inline constexpr class_t kExtraSplitSLS = "ExtraSplitSls"; // TODO
inline constexpr class_t kUnk = "Unk";
inline constexpr class_t kLowRq = "LowRQ";
inline constexpr class_t kTooShort = "TooShort";

// Terminal (end) of a read
enum class ReadTerminal { LEFT, RIGHT };

const std::unordered_map<char, char> kDNAComplement = {
        {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};

}
#endif //MARTI_INCLUDE_MARTI_DEFINITIONS_H_
