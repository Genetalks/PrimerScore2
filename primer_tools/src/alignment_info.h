#ifndef ALIGNMENT_INFO_H
#define ALIGNMENT_INFO_H

#include <string>
#include <vector>
#include <memory>
#include <htslib/sam.h>
#include <cassert>
#include <utility>
#include "IntervalTree.h"

namespace pt
{

    class alignment_info
    {
    public:
        alignment_info(const bam1_t *b, int32_t template_len);

        int32_t get_template_end() const
        {
            return this->qe;
        }
        int32_t get_template_start() const
        {
            return this->qs;
        }

        int64_t get_ref_start() const
        {
            return this->rs;
        }

        int64_t get_ref_end() const
        {
            return this->re;
        }

        bool is_rev_strand() const
        {
            return this->is_rev;
        }

        const std::string &get_chromosome() const
        {
            return this->chr;
        }

        const std::string &get_ref_seq() const
        {
            return this->seq;
        }

        const std::vector<uint32_t> &get_cigar() const
        {
            return this->cigar;
        }

        const std::string &get_aligned_ref_seq() const
        {
            return this->aligned_ref_seq;
        }

        const std::string &get_aligned_str() const
        {
            return this->aligned_str;
        }

        const std::string &get_aligned_seq() const
        {
            return this->aligned_seq;
        }

        void calculate_aligned_blocks(const std::string &template_seq);

        std::string to_string();

    private:
        void initialize_from_bam_record(const bam1_t *b, int32_t template_len);

    private:
        int32_t tid; // template id
        int32_t qs;
        int32_t qe;
        std::string chr;
        int64_t rs;
        int64_t re;
        std::string seq; // reference sequence
        bool is_rev;
        std::vector<uint32_t> cigar;
        std::string aligned_ref_seq;
        std::string aligned_str;
        std::string aligned_seq;
    };

    typedef std::shared_ptr<alignment_info> alignment_info_ptr;

    typedef IntervalTree<int32_t, alignment_info *> interval_tree_t;
    typedef std::shared_ptr<interval_tree_t> interval_tree_ptr;
    typedef interval_tree_t::interval interval_t;
    typedef interval_tree_t::interval_vector interval_vector_t;
    typedef interval_tree_t::value_vector value_vector_t;

    // This class is used for one primer template
    class alignment_info_db
    {
    public:
        alignment_info_db(int32_t template_len);

        alignment_info_ptr get_alignment_info(uint32_t idx) const
        {
            assert(idx < alignment_infos.size());
            return this->alignment_infos[idx];
        }
        value_vector_t query(int32_t primer_start, int32_t primer_end) const
        {
            auto v = this->query_index->findOverlappingValues(primer_start, primer_end);
            return v;
        }
        size_t size() const
        {
            return this->alignment_infos.size();
        }

        void init_alignment_info_from_bam_records(const std::vector<std::shared_ptr<bam1_t>> &alignments);
        void calculate_aligned_blocks(const std::string &template_seq);
        void build_query_tree();

        void clear()
        {
            this->query_index.reset();
            this->alignment_infos.clear();
        }

    private:
        std::vector<alignment_info_ptr> alignment_infos; // for one template
        interval_tree_ptr query_index;
        int32_t template_len;
    };

    typedef std::shared_ptr<alignment_info_db> alignment_info_db_ptr;

}

#endif
