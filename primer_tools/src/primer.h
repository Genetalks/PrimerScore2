#ifndef PRIMER_H
#define PRIMER_H

#include <string>
#include <cstdint>
#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <array>
#include "htslib/sam.h"
#include "alignment_info.h"
#include "tm_calculator.h"
#include "linked_map.h"

namespace pt
{

    class primer
    {
    public:
        primer(const std::string &id, const std::string &tid, int32_t start, int32_t end, bool is_rev, const std::string &seq);

        const std::string &get_id() const
        {
            return this->id;
        }
        const std::string &get_tid() const
        {
            return this->tid;
        }
        const std::string &get_seq() const
        {
            return this->seq;
        }
        int32_t get_len() const
        {
            return this->seq.size();
        }
        int32_t get_start() const
        {
            return this->start;
        }
        int32_t get_end() const
        {
            return this->end;
        }
        size_t size() const
        {
            return end >= start ? end - start + 1 : start - end + 1;
        }

        bool is_rev() const
        {
            return this->is_rev_strand;
        };

    private:
        std::string id;
        std::string tid;    // template id
        int32_t start;      // start pos of primer on template
        int32_t end;        // end pos of primer on template
        bool is_rev_strand; // is negative strand
        std::string seq;    // primer sequence
    };

    typedef std::shared_ptr<primer> primer_ptr;

    class primer_interval_t
    {
    public:
        void set(int64_t pos, const std::string &chr, const std::string *pseq, const std::string &astr, double tm, char strand);
        std::string to_string() const;

        bool is_retained() const
        {
            return this->retained;
        }

    public:
        static const int16_t max_cigar_op_len = 0xfff;
        static const int16_t cgr_op_mask = 0xf;
        static const int16_t cgr_len_shift = 4;
        static const int32_t alpha = 128;
        static const char unknown_op;
        static const std::string cgr_str;
        static const std::array<uint8_t, alpha> init_cgr_map();
        static const std::array<uint8_t, alpha> cgr_map;

    private:
        int64_t pos = -1; // pos5
        std::string chr;
        const std::string *pseq = nullptr;
        std::vector<int16_t> astr;
        double tm = -1.0;
        char strand;
        bool contains_deletion = false;
        bool retained = false;
    };

    typedef std::vector<primer_interval_t> primer_interval_vec_t;
    typedef std::unordered_map<std::string, double> tm_cache_t;

    class primer_template
    {
    public:
        primer_template(const std::string &tid) : tid(tid) {}
        void add_primer(const primer_ptr &p)
        {
            this->primers.emplace_back(p);
        }
        void add_alignment(const std::shared_ptr<bam1_t> b)
        {
            this->raw_alignments.emplace_back(b);
        }

        const std::string &get_tid() const
        {
            return this->tid;
        }

        const std::string &get_seq() const
        {
            return this->seq;
        }

        void set_seq(const std::string &seq)
        {
            this->seq = seq;
        }

        const std::vector<primer_ptr> &get_primers() const
        {
            return this->primers;
        }

        const std::vector<std::shared_ptr<bam1_t>> &get_alignments() const
        {
            return this->raw_alignments;
        }

        void clear_alignments()
        {
            this->raw_alignments.clear();
        }

        alignment_info_db_ptr get_db()
        {
            return this->candidate_alignments;
        }

        const std::vector<value_vector_t> &get_hits() const
        {
            return this->primers_hits;
        }

        const std::vector<primer_interval_vec_t> &get_primer_interval() const
        {
            return this->primers_intervals;
        }

        void clear_hits()
        {
            this->primers_hits.clear();
            this->candidate_alignments->clear();
        }

        void clear_intervals()
        {
            this->primers_intervals.clear();
        }

        void initialize_alignment_db(bam_hdr_t *h);
        void query_primers(int32_t min_overlap = 12, int32_t max_non_overlap_3_end = 5, bool probe_mode = false);
        // calcaute one hits
        bool calculate_primer_interval(tm_cache_t &alignstr_filter_map, const primer &p, const alignment_info *template_alignment_info, primer_interval_t &primer_hits_interval, pt::tm_calculator *tm_calc, double min_bound_tm);
        void calculate_primer_intervals(pt::tm_calculator *tm_calc, double min_bound_tm, int32_t s = -1, int32_t e = -1);

    public:
        static void remove(std::string &str, int32_t ix, int32_t rl)
        {
            int32_t len = str.length();
            for (int32_t i = ix; i < len - rl; i++)
            {
                str[i] = str[i + rl];
            }
            str.resize(len - rl);
        }

        static void reverse(std::string &str)
        {
            int32_t len = str.length();
            for (int32_t i = 0; i < len / 2; i++)
            {
                char t = str[i];
                str[i] = str[len - 1 - i];
                str[len - 1 - i] = t;
            }
        }
        static void revcom(std::string &str)
        {
            for (char &c : str)
            {
                switch (c)
                {
                case 'A':
                    c = 'T';
                    break;
                case 'T':
                    c = 'A';
                    break;
                case 'C':
                    c = 'G';
                    break;
                case 'G':
                    c = 'C';
                    break;
                }
            }
            reverse(str);
        }

        int64_t get_total_hits()
        {
            if (this->total_hits == -1)
            {
                int64_t sum = 0;
                for (auto &hit : this->primers_hits)
                {
                    sum += hit.size();
                }
                this->total_hits = sum;
            }
            return this->total_hits;
        }

        void mark_split()
        {
            this->require_split = true;
        }

        bool is_marked_split() const
        {
            return this->require_split;
        }

    private:
        std::string tid;                 // template id
        std::string seq;                 // template sequence
        std::vector<primer_ptr> primers; // primer list on this template
        std::vector<std::shared_ptr<bam1_t>> raw_alignments;
        alignment_info_db_ptr candidate_alignments;
        std::vector<value_vector_t> primers_hits;
        std::vector<primer_interval_vec_t> primers_intervals;
        int64_t total_hits = -1;
        bool require_split = false;
    };

    typedef std::shared_ptr<primer_template> primer_template_ptr;

    typedef pt::linked_map<std::string, primer_template_ptr> primer_template_db_t;
    typedef std::shared_ptr<primer_template_db_t> primer_template_db_ptr;
};

#endif
