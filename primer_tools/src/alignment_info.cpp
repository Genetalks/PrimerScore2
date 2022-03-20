#include "alignment_info.h"
#include <tbb/parallel_for.h>
#include <sstream>
#include <htslib/sam.h>
#include <cassert>

pt::alignment_info::alignment_info(const bam1_t *b, int32_t template_len):
    tid(b->core.tid),
    chr(bam_get_qname(b)),
    is_rev(bam_is_rev(b))
    
{
    this->initialize_from_bam_record(b, template_len);
}

std::string pt::alignment_info::to_string(){
    std::ostringstream oss;

    oss << "[" << this->qs << ","
        << this->qe << "]|"
        << this->chr << ":"
        << this->rs << "-"
        << this->re << "|"
        << (this->is_rev ? "-" : "+") << "|";
    for (auto &c : cigar){
        oss << bam_cigar_oplen(c) << bam_cigar_opchr(c);
    }
    oss << "\n" << this->aligned_ref_seq;
    oss << "\n" << this->aligned_str;
    oss << "\n" << this->aligned_seq; 

    return oss.str();
}


void pt::alignment_info::initialize_from_bam_record(const bam1_t *b, int32_t template_len){

    int32_t n_cigar = b->core.n_cigar;
    this->cigar.reserve(n_cigar);
    uint32_t *cgr = bam_get_cigar(b);
	
	this->qs = b->core.pos+1;
	this->qe = bam_endpos(b);
	
	if (this->qs > 1){
		this->cigar.emplace_back((this->qs - 1) << 4 | BAM_CHARD_CLIP);
	}
	for (int32_t i = 0; i < n_cigar; ++i){
		auto op = bam_cigar_op(cgr[i]);
		auto len = bam_cigar_oplen(cgr[i]);
		switch (op){
		case BAM_CMATCH:
			this->cigar.emplace_back(cgr[i]);
			break;
		case BAM_CINS:
			this->cigar.emplace_back(len << 4 | BAM_CDEL);
			break;
		case BAM_CDEL:
			this->cigar.emplace_back(len << 4 | BAM_CINS);
			break;
		case BAM_CHARD_CLIP:
			break;
		default:
			assert(0);
			break;
		}
	}
	if (this->qe < template_len){
		this->cigar.push_back((template_len - qe - 1) << 4 | BAM_CHARD_CLIP);
	}

    if (!this->is_rev){
        this->rs = bam_cigar_op(cgr[0]) == BAM_CHARD_CLIP ? bam_cigar_oplen(cgr[0])+1 : 0;
        this->re = this->rs + bam_cigar2qlen(n_cigar, cgr) - 1;
       
    }else{
        this->rs = bam_cigar_op(cgr[n_cigar - 1]) == BAM_CHARD_CLIP ? bam_cigar_oplen(cgr[n_cigar - 1])+1 : 0;
        this->re = this->rs + bam_cigar2qlen(n_cigar, cgr) - 1;
    }

    // get sequence
    auto s = bam_get_seq(b);
    this->seq.resize(b->core.l_qseq);
    for (int32_t i = 0; i < b->core.l_qseq; ++i){
        this->seq[i] = seq_nt16_str[bam_seqi(s, i)];
    }
}

void pt::alignment_info::calculate_aligned_blocks(const std::string &template_seq){

    int32_t tlen = 0;
    int32_t rlen = 0;
    for (auto &cgr : this->cigar){
        int32_t op_len = bam_cigar_oplen(cgr);
        int32_t op = bam_cigar_op(cgr);

        switch(op){
            case BAM_CMATCH:
            {
                assert(rlen + op_len <= this->seq.size());
                assert(tlen + op_len <= template_seq.size());
                std::string s1{this->seq.substr(rlen, op_len)};
                std::string s2{template_seq.substr(tlen, op_len)};
                this->aligned_ref_seq.append(s1);
                this->aligned_seq.append(s2);
                for (int32_t i = 0; i < op_len; ++i){
                    if (s1[i] == s2[i]){
                        this->aligned_str.append(1, '|');
                    }else{
                        this->aligned_str.append(1, '*');
                    }
                }
                tlen += op_len;
                rlen += op_len;
                break;
            }
            case BAM_CINS:
            {
                assert(tlen + op_len <= template_seq.size());
                this->aligned_ref_seq.append(op_len, '-');
                this->aligned_seq.append(template_seq.substr(tlen, op_len));
                this->aligned_str.append(op_len, ' ');
                tlen += op_len;
                break;
            }
            case BAM_CDEL:
            {
                assert(rlen + op_len <= this->seq.size());
                this->aligned_ref_seq.append(this->seq.substr(rlen, op_len));
                this->aligned_seq.append(op_len, '-');
                this->aligned_str.append(op_len, ' ');
                rlen += op_len;
                break;
            }
            case BAM_CHARD_CLIP:
                tlen += op_len;
                break;
            default:
                break;
        }
    }

}

pt::alignment_info_db::alignment_info_db(int32_t template_len):template_len(template_len) {}

void pt::alignment_info_db::init_alignment_info_from_bam_records(const std::vector<std::shared_ptr<bam1_t> > &alignments){
    this->alignment_infos.resize(alignments.size());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, this->alignment_infos.size()), [&](tbb::blocked_range<size_t> r){
        for (size_t i = r.begin(); i < r.end(); ++i){
            this->alignment_infos[i] = std::make_shared<alignment_info>(alignments[i].get(), this->template_len);
        }
    });

}

void pt::alignment_info_db::calculate_aligned_blocks(const std::string &template_seq){
    tbb::parallel_for(tbb::blocked_range<size_t>(0, this->alignment_infos.size()), [&](tbb::blocked_range<size_t> r){
        for (size_t i = r.begin(); i < r.end(); ++i){
            this->alignment_infos[i]->calculate_aligned_blocks(template_seq);
        }
    });
}

void pt::alignment_info_db::build_query_tree(){

    interval_vector_t aligned_intervals;
    aligned_intervals.reserve(this->alignment_infos.size());

    for(auto &a : this->alignment_infos){
        aligned_intervals.push_back(interval_t{a->get_template_start(), a->get_template_end(), a});
    }

    this->query_index.reset(new interval_tree_t{std::move(aligned_intervals)});
}
