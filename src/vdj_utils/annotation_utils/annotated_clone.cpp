#include <verify.hpp>

#include "annotated_clone.hpp"

#include <seqan/stream.h>
#include "seqan/translation.h"

namespace annotation_utils {
    std::ostream& operator<<(std::ostream& out, const StructuralRegion &region) {
        if(region == StructuralRegion::CDR1)
            out << "CDR1";
        else if(region == StructuralRegion::CDR2)
            out << "CDR2";
        else if(region == StructuralRegion::CDR3)
            out << "CDR3";
        else
            out << "Unknown region";
        return out;
    }

    void AnnotatedClone::CheckRangeConsistencyFatal(CDRRange range) {
        VERIFY(range.Full());
        VERIFY_MSG(range.start_pos < read_.length() and range.end_pos < read_.length(), "Start pos (" <<
                range.start_pos << ") or end pos (" << range.end_pos << ") exceeds read length " << read_.length());
    }

    void AnnotatedClone::UpdateStructuralRegion(StructuralRegion region, CDRRange range) {
        //TRACE("Updating " << region << " by range " << range);
        CheckRangeConsistencyFatal(range);
        region_range_map_[region] = range;
        seqan::Dna5String cdr_seq;
        if(range.Valid())
            cdr_seq = seqan::infixWithLength(read_.seq, range.start_pos, range.length());
        region_string_map_[region] = cdr_seq;
    }

    void AnnotatedClone::Initialize(CDRLabeling cdr_labeling) {
        UpdateStructuralRegion(StructuralRegion::CDR1, cdr_labeling.cdr1);
        UpdateStructuralRegion(StructuralRegion::CDR2, cdr_labeling.cdr2);
        UpdateStructuralRegion(StructuralRegion::CDR3, cdr_labeling.cdr3);
    }

    void AnnotatedClone::InitializeAASeq() {
        using namespace seqan;
        StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > aa_seqs;
        translate(aa_seqs, read_.seq, SINGLE_FRAME);
        aa_read_seq_ = aa_seqs[0];
        // todo: compute productiveness and in-frame accurately
        productive_ = false;
        in_frame_ = false;
    }

    char AnnotatedClone::GetAminoAcidByPos(const seqan::String<seqan::AminoAcid> &aa_seq, size_t nucl_pos,
                                           unsigned orf) const {
        if(nucl_pos < orf)
            return '-';
        size_t pos = (nucl_pos - orf) / 3;
        if(pos < seqan::length(aa_seq))
            return aa_seq[pos];
        else if(pos == seqan::length(aa_seq))
            return '-';
        VERIFY_MSG(false, "Unknown position " << pos << " of amino acid sequence " << aa_seq);
        return '-';
    }

    void AnnotatedClone::InitializeSHMs(germline_utils::SegmentType segment_type) {
        const alignment_utils::ImmuneGeneReadAlignment& alignment =
                (segment_type == germline_utils::SegmentType::VariableSegment) ? v_alignment_ : j_alignment_;
        GeneSegmentSHMs& shms = (segment_type == germline_utils::SegmentType::VariableSegment) ? v_shms_ : j_shms_;
        std::cout << read_ << std::endl;
        std::cout << aa_read_seq_ << std::endl;
        std::cout << alignment.subject() << ", ORF: " << alignment.subject().ORF() << std::endl;
        auto gene_row = seqan::row(alignment.Alignment(), 0);
        auto read_row = seqan::row(alignment.Alignment(), 1);
        std::cout << gene_row << std::endl;
        std::cout << read_row << std::endl;
        std::cout << alignment.subject().aa_seq() << std::endl;
        for(size_t i = alignment.RealStartAlignmentPos(); i <= alignment.RealEndAlignmentPos(); i++) {
            if(gene_row[i] != read_row[i]) {
                size_t real_read_pos = seqan::toSourcePosition(read_row, i);
                size_t real_gene_pos = seqan::toSourcePosition(gene_row, i);
                SHM shm(real_read_pos, real_gene_pos, gene_row[i], read_row[i],
                        GetAminoAcidByPos(alignment.subject().aa_seq(), real_gene_pos, alignment.subject().ORF()),
                        GetAminoAcidByPos(aa_read_seq_, real_read_pos, 0));
                shms.AddSHM(shm);
            }
        }
        std::cout << shms;
    }

    bool AnnotatedClone::RegionIsEmpty(StructuralRegion region) const {
        if (region_range_map_.find(region) == region_range_map_.end())
            return true;
        return seqan::length(GetRegionString(region)) == 0;
    }

    seqan::Dna5String AnnotatedClone::GetRegionString(StructuralRegion region) const {
        if(region_range_map_.find(region) == region_range_map_.end())
            return seqan::Dna5String();
        //VERIFY_MSG(region_range_map_.find(region) != region_range_map_.end(),
        //           "Clone does not have information about region " << region);
        return region_string_map_.at(region);
    }

    CDRRange AnnotatedClone::GetRangeByRegion(StructuralRegion region) const {
        VERIFY_MSG(!RegionIsEmpty(region), "Clone does not have information about region " << region);
        return region_range_map_.at(region);
    }

    std::ostream& operator<<(std::ostream& out, const AnnotatedClone &obj) {
        out << obj.Read() << std::endl;
        if(!obj.RegionIsEmpty(StructuralRegion::CDR1))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR1) << std::endl;
        if(!obj.RegionIsEmpty(StructuralRegion::CDR2))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR2) << std::endl;
        if(!obj.RegionIsEmpty(StructuralRegion::CDR3))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR3) << std::endl;
        return out;
    }
}