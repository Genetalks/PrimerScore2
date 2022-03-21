#ifndef PRIMER_H
#define PRIMER_H

#include <string>
#include <cstdint>
#include <memory>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include "htslib/sam.h"
#include <tbb/parallel_for.h>
#include "alignment_info.h"

namespace pt {

class primer {
public:
    primer(const std::string &id, const std::string &tid, int32_t start, int32_t end, bool is_rev, const std::string &seq):
        id(id),
        tid(tid),
        start(start),
        end(end),
        is_rev_strand(is_rev),
		seq(seq)
    {
    }

    const std::string &get_id() const {
        return this->id;
    }
    const std::string &get_tid() const {
        return this->tid;
    }
    const std::string &get_seq() const {
        return this->seq;
    }
    int32_t get_len() const {
        return this->seq.size();
    }
    int32_t get_start() const {
        return this->start;
    }
    int32_t get_end() const {
        return this->end;
    }
    size_t size() const {
        return end >= start ? end - start + 1 : start - end + 1;
    }

    bool is_rev() const {
        return this->is_rev_strand;
    };

private:
    std::string id;
    std::string tid; // template id
    int32_t start;  // start pos of primer on template
    int32_t end;    // end pos of primer on template 
    bool is_rev_strand;  // is negative strand
    std::string seq; // primer sequence
};

typedef std::shared_ptr<primer> primer_ptr;

struct primer_interval_t {
    int32_t ps;
    int32_t pe;
    int64_t rs;
    int64_t re;
	char strand;
    std::string chr;
    std::string pseq;
    std::string rseq;
    std::string astr;

    std::string to_string() const {
        std::ostringstream oss;

        oss << this->ps << "\t"
            << this->pe << "\t"
            << this->strand << "\t"
            << this->chr << "\t"
            << this->rs << "\t"
            << this->re << "\t"
            << this->pseq << "\t"
            << this->rseq << "\t"
            << this->astr;

        return oss.str();
    }
};

typedef std::vector<primer_interval_t> primer_interval_vec_t;

class primer_template {
public:
    primer_template(const std::string &tid):tid(tid){}
    void add_primer(const primer_ptr &p){
        this->primers.emplace_back(p);
    }
    void add_alignment(const std::shared_ptr<bam1_t> b){
        this->raw_alignments.emplace_back(b);
    }

    const std::string &get_tid() const {
        return this->tid;
    }

    const std::string &get_seq() const {
        return this->seq;
    }

    void set_seq(const std::string &seq) {
        this->seq = seq;
    }

    const std::vector<primer_ptr> &get_primers() const {
        return this->primers;
    }
    
    const std::vector<std::shared_ptr<bam1_t> > &get_alignments() const {
        return this->raw_alignments;
    }

    void clear_alignments(){
        this->raw_alignments.clear();
    }

    void initialize_alignment_db(bam_hdr_t *h){
        int32_t idx = sam_hdr_name2tid(h, this->tid.c_str());

        if (idx < 0){
            std::cerr << "No alignments of template found in bam file: " << this->tid ;
            exit(EXIT_FAILURE);
        }
        
        int32_t template_len = h->target_len[idx];
        this->candidate_alignments = std::make_shared<alignment_info_db>(template_len);
    }

    alignment_info_db_ptr get_db(){
        return this->candidate_alignments;
    }

    void query_primers(int32_t min_overlap = 12, int32_t max_non_overlap_3_end = 5, bool probe_mode = false){
        this->primers_hits.resize(this->primers.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, this->primers.size()), [&](tbb::blocked_range<size_t> r){
            for (size_t i = r.begin(); i < r.end(); ++i){
				auto &p = this->primers[i];
				bool is_rev = p->is_rev();
				int32_t ps = p->get_start();
				int32_t pe = p->get_end();

				// query
                this->primers_hits[i] = std::move(this->candidate_alignments->query(ps, pe));

				auto &hits = this->primers_hits[i];
				int32_t k = 0;
				for (size_t j = 0; j < hits.size(); ++j){
					int32_t overlap = std::min(hits[j].stop - this->primers[i]->get_start(), this->primers[i]->get_end() - hits[j].start);
                    if (overlap >= min_overlap){
                        if(probe_mode || 
                           (!probe_mode 
                             && ((is_rev && ps >= hits[j].start - max_non_overlap_3_end) 
                               || (!is_rev && pe <= hits[j].stop + max_non_overlap_3_end)))){ 
				   
							hits[k++] = hits[j];
						}
					}
				}
				hits.resize(k);
            }
        });
    }

    // calcaute one hits
    void calculate_primer_interval(const primer_ptr &p, const alignment_info_ptr &template_alignment_info, primer_interval_t &primer_hits_interval){
        int32_t ps = p->get_start()+1; //1-based
        int32_t pe = p->get_end()+1;
        int32_t plen = p->get_len();
		bool is_prev = p->is_rev();
        std::string seqp = p->get_seq();
		int32_t ts = template_alignment_info->get_template_start();
        int32_t te = template_alignment_info->get_template_end();
        int64_t rs = template_alignment_info->get_ref_start();
        int64_t re = template_alignment_info->get_ref_end();
        primer_hits_interval.chr = template_alignment_info->get_chromosome();
		bool is_rev = template_alignment_info->is_rev_strand();
		std::cout<<">>"<<p->get_id()<<"\t"<<seqp<<"\t"<<ps<<"\t"<<pe<<"\t"<<plen<<"\t"<<is_prev<<"\t";
		
        // overlap section
        int32_t s = std::max(ps, ts); // start pos of alignment on template
        int32_t e = std::min(pe, te); // end pos of alignment on template
        int32_t asp = s-ps+1; // start pos of alignment on primer, 1-based
        int32_t aep = e-ps+1; // end pos of alignment on primer, 1-based
		primer_hits_interval.ps = asp;
        primer_hits_interval.pe = aep;
		primer_hits_interval.strand = is_rev==is_prev? '+': '-';
        s -= ts; // start pos of primer on aligned template region, 0-based
        e -= ts; // end pos of primer on aligned template region, 0-based
		std::cout<<asp<<"\t"<<aep<<"\t"<<s<<"\t"<<e<<std::endl;
		std::cout<<ts<<"\t"<<te<<"\t"<<rs<<"\t"<<re<<"\t"<<is_rev<<"\t"<<std::endl;
		
        assert(e >= s);
		
        // s and e offset on aligned query string
		// t: pos on template; r: pos on ref
        int32_t as = 0, ae = 0, t = 0, r = 0;
		
		// aligned template info
        auto &astring = template_alignment_info->get_aligned_str(); 
        auto &ars = template_alignment_info->get_aligned_ref_seq();
        auto &ats = template_alignment_info->get_aligned_seq();
        std::cout<<ats<<"\n"<<astring<<"\n"<<ars<<std::endl<<std::endl;
		for (size_t i = 0; i < astring.size(); ++i){
            if (astring[i] != ' '){ // MATCH and mismatch
                ++t;
                ++r;
            }else{ // ins and del
                if (ars[i] == '-'){
                    ++t;
                }else{
                    ++r;
                }
            }
            if (t == s+1){
                as = i;
				if(!is_rev){
					primer_hits_interval.rs = rs + r - 1;
				}else{
					primer_hits_interval.re = re - r + 1;
				}
            }
            if (t == e+1){
                ae = i;
				if(!is_rev){
					primer_hits_interval.re = rs + r - 1;
				}else{
					primer_hits_interval.rs = re - r + 1;
				}
			}
        }
		
		std::string prefixp = seqp.substr(0,asp-1);
		std::string prefixr (asp-1, 'N'); 
		std::string prefixa (asp-1, '#'); 
		std::string sufixp = seqp.substr(aep);
		std::string sufixr ((plen-aep), 'N'); 
		std::string sufixa ((plen-aep), '#'); 
		primer_hits_interval.pseq = prefixp+ats.substr(as, ae - as + 1)+sufixp;
        primer_hits_interval.rseq = prefixr+ars.substr(as, ae - as + 1)+sufixr;
        primer_hits_interval.astr = prefixa+astring.substr(as, ae - as + 1)+sufixa;
		std::cout<<prefixp<<","<<sufixp<<";"<<prefixr<<","<<sufixr<<std::endl;
		std::cout<<primer_hits_interval.pseq<<"\n"<<primer_hits_interval.astr<<"\n"<<primer_hits_interval.rseq<<std::endl;
	}

    void calculate_primer_intervals(){
        this->primers_intervals.resize(this->primers.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, this->primers.size()), [&](tbb::blocked_range<size_t> r){
            for (size_t i = r.begin(); i < r.end(); ++i){
                auto &hits = this->primers_hits[i];
                auto &intervals =  this->primers_intervals[i];
                intervals.resize(hits.size());
                for (size_t j = 0; j < hits.size(); ++j){
                    calculate_primer_interval(this->primers[i], hits[j].value, intervals[j]);
                }	
            }
        });
    }

    const std::vector<interval_vector_t> &get_hits() const {
        return this->primers_hits;
    }

    const std::vector<primer_interval_vec_t> &get_primer_interval() const {
        return this->primers_intervals;
    }

private:
    std::string tid;  // template id
    std::string seq;  // template sequence
    std::vector<primer_ptr> primers; // primer list on this template
    std::vector<std::shared_ptr<bam1_t> > raw_alignments;
    alignment_info_db_ptr candidate_alignments;
    std::vector<interval_vector_t> primers_hits;
    std::vector<primer_interval_vec_t> primers_intervals;
};

typedef std::shared_ptr<primer_template> primer_template_ptr;

typedef std::map<std::string, primer_template_ptr> primer_template_db_t;
typedef std::shared_ptr<primer_template_db_t> primer_template_db_ptr;
};


#endif

