#include <sstream>
#include <sstream>
#include <iostream>
#include <glog/logging.h>
#include <atomic>
#include <chrono>
#include "primer.h"
#include "tbb/parallel_for.h"
#include "tbb/task_arena.h"

namespace pt
{

    // primer
    primer::primer(const std::string &id, const std::string &tid, int32_t start, int32_t end, bool is_rev, const std::string &seq) : id(id),
                                                                                                                                     tid(tid),
                                                                                                                                     start(start),
                                                                                                                                     end(end),
                                                                                                                                     is_rev_strand(is_rev),
                                                                                                                                     seq(seq)
    {
    }

    // primer interval
    const char primer_interval_t::unknown_op{static_cast<char>(255)};
    const std::string primer_interval_t::cgr_str{"|*^-#"};
    const std::array<uint8_t, primer_interval_t::alpha> primer_interval_t::cgr_map{init_cgr_map()};
    const std::array<uint8_t, primer_interval_t::alpha> primer_interval_t::init_cgr_map()
    {
        std::array<uint8_t, alpha> cgr_map;
        std::fill(cgr_map.begin(), cgr_map.end(), unknown_op);
        cgr_map['|'] = 0;
        cgr_map['*'] = 1;
        cgr_map['^'] = 2;
        cgr_map['-'] = 3;
        cgr_map['#'] = 4;

        return cgr_map;
    }

    void primer_interval_t::set(int64_t pos, const std::string &chr, const std::string *pseq, const std::string &astr, double tm, char strand)
    {
        this->pos = pos;
        this->chr = chr;
        this->pseq = pseq;
        this->tm = tm;
        this->strand = strand;

        // calculate alignment cigar
        bool has_deletion = false;
        int32_t op_len = 0;
        char last_char = unknown_op;
        for (auto &c : astr)
        {
            if (c != last_char)
            {
                if (last_char != unknown_op)
                {
                    assert(op_len <= max_cigar_op_len);
                    this->astr.emplace_back(op_len << cgr_len_shift | cgr_map[last_char]);
                    if (last_char == '-')
                    {
                        has_deletion = true;
                    }
                }
                op_len = 1;
            }
            else
            {
                ++op_len;
            }
            last_char = c;
        }
        assert(op_len <= max_cigar_op_len);
        this->astr.emplace_back(op_len << cgr_len_shift | cgr_map[last_char]);
        if (last_char == '-')
        {
            has_deletion = true;
        }
        this->astr.shrink_to_fit();
        this->contains_deletion = has_deletion;
        this->retained = true;
    }

    std::string primer_interval_t::to_string() const
    {
        std::ostringstream oss;

        oss << this->strand << "\t"
            << this->chr << "\t"
            << this->pos << "\t";
        // print seq
        if (!this->contains_deletion)
        {
            oss << *this->pseq << "\t";
        }
        else
        {
            int32_t offset = 0;
            for (auto &cgr : this->astr)
            {
                auto op_len = cgr >> cgr_len_shift;
                auto op = cgr_str[cgr & cgr_op_mask];
                if (op != '-')
                {
                    oss << this->pseq->substr(offset, op_len);
                    offset += op_len;
                }
                else
                {
                    std::string fill_chars;
                    fill_chars.append(op_len, op);
                    oss << fill_chars;
                }
            }
            oss << "\t";
        }
        oss << this->tm << "\t"
            << this->pseq->back() << "\t"; // end3base
        // print aln str
        for (auto &cgr : this->astr)
        {
            std::string block;
            block.append(cgr >> cgr_len_shift, cgr_str[cgr & cgr_op_mask]);
            oss << block;
        }

        return oss.str();
    }

    // primer template
    void primer_template::initialize_alignment_db(bam_hdr_t *h)
    {
        int32_t idx = sam_hdr_name2tid(h, this->tid.c_str());

        if (idx < 0)
        {
            std::cerr << "No alignments of template found in bam file: " << this->tid;
            exit(EXIT_FAILURE);
        }

        int32_t template_len = h->target_len[idx];
        this->candidate_alignments = std::make_shared<alignment_info_db>(template_len);
    }

    void primer_template::query_primers(int32_t min_overlap, int32_t max_non_overlap_3_end, bool probe_mode)
    {
        this->primers_hits.resize(this->primers.size());
        auto db = this->candidate_alignments.get();
        tbb::parallel_for(tbb::blocked_range<size_t>(0, this->primers.size()), [&](tbb::blocked_range<size_t> r)
                          {
        for (size_t i = r.begin(); i < r.end(); ++i){
            auto p = this->primers[i].get();
            bool is_rev = p->is_rev();
            int32_t ps = p->get_start();
            int32_t pe = p->get_end();

            // query
            this->primers_hits[i] = std::move(db->query(ps, pe));

            auto &hits = this->primers_hits[i];
            int32_t k = 0;
            for (size_t j = 0; j < hits.size(); ++j){
                int32_t hits_start = hits[j]->get_template_start();
                int32_t hits_stop = hits[j]->get_template_end();
                int32_t overlap = std::min(hits_stop - p->get_start(), p->get_end() - hits_start);
                if (overlap >= min_overlap){
                    if(probe_mode || 
                       (!probe_mode 
                         && ((is_rev && ps >= hits_start - max_non_overlap_3_end) 
                           || (!is_rev && pe <= hits_stop + max_non_overlap_3_end)))){ 
               
                        hits[k++] = hits[j];
                    }
                }
            }
            hits.resize(k);
        } });
    }

    void primer_template::calculate_primer_intervals(pt::tm_calculator *tm_calc, double min_bound_tm, int32_t s, int32_t e)
    {
        // s: start primer index, e: end primer index
        s = s == -1 ? 0 : s;
        e = e == -1 ? static_cast<int32_t>(this->primers.size()) : e;
        if (e <= s)
        {
            LOG(FATAL) << "Block index is  not properly set.";
        }

        std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
        this->primers_intervals.resize(e - s);

        std::atomic<int64_t> total_hits{0};
        std::atomic<int64_t> total_intervals{0};
        tbb::parallel_for(tbb::blocked_range<int32_t>(s, e), [&](tbb::blocked_range<int32_t> r)
                          {
        for (int32_t i = r.begin(); i < r.end(); ++i){
            auto &hits = this->primers_hits[i];
            auto &intervals =  this->primers_intervals[i-s];
            tm_cache_t alignstr_filter_map;
            intervals.clear(); 
            intervals.resize(hits.size());
            size_t n=0;
            bool is_filter;
            for (size_t j = 0; j < hits.size(); ++j){
                    is_filter = calculate_primer_interval(alignstr_filter_map, *this->primers[i], hits[j], intervals[n], tm_calc, min_bound_tm);
                if(!is_filter){
                    n++;
                }
            }
            intervals.resize(n);
            total_hits += hits.size();
            total_intervals += n;
        } });

        std::chrono::duration<double> etime = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        LOG(INFO) << "[T]: " << this->tid << ", "
                  << "[P]: " << this->primers.size() << ", "
                  << "[B]: [" << s << ", " << e << "),"
                  << "[A]: " << this->candidate_alignments->size() << ", "
                  << "[H]: " << total_hits << ", "
                  << "[I]: " << total_intervals << ", "
                  << "[E]: " << etime.count() << "s.";
    }

    bool primer_template::calculate_primer_interval(tm_cache_t &alignstr_filter_map, const primer &p, const alignment_info *template_alignment_info, primer_interval_t &primer_hits_interval, pt::tm_calculator *tm_calc, double min_bound_tm)
    {
        int32_t ps = p.get_start() + 1; // 1-based
        int32_t pe = p.get_end() + 1;
        int32_t plen = p.get_len();
        bool is_prev = p.is_rev();
        const std::string &seqp = p.get_seq();
        int32_t ts = template_alignment_info->get_template_start();
        int32_t te = template_alignment_info->get_template_end();
        const std::string &chr = template_alignment_info->get_chromosome();
        int64_t rs = template_alignment_info->get_ref_start();
        int64_t re = template_alignment_info->get_ref_end();
        bool is_rev = template_alignment_info->is_rev_strand();
        char strand = is_rev == is_prev ? '+' : '-';
#ifdef DEBUG
        std::cout << ">>" << p.get_id() << "\t" << seqp << "\t" << ps << "\t" << pe << "\t" << plen << "\t" << is_prev << "\t";
#endif

        // overlap section
        int32_t s = std::max(ps, ts); // start pos of alignment on template
        int32_t e = std::min(pe, te); // end pos of alignment on template
        int32_t asp = s - ps + 1;     // start pos of alignment on primer, 1-based
        int32_t aep = e - ps + 1;     // end pos of alignment on primer, 1-based
        s -= ts;                      // start pos of primer on aligned template region, 0-based
        e -= ts;                      // end pos of primer on aligned template region, 0-based
        assert(e >= s);

#ifdef DEBUG
        std::cout << asp << "\t" << aep << "\t" << s << "\t" << e << std::endl;
        std::cout << ts << "\t" << te << "\t" << rs << "\t" << re << "\t" << is_rev << "\t" << std::endl;
#endif

        // s and e offset on aligned query string
        // t: pos on template; r: pos on ref
        int32_t as = -1, ae = -1, t = 0, r = 0;

        int32_t rstart = -1;
        int32_t rend = -1;
        // aligned template info
        auto &astring = template_alignment_info->get_aligned_str();
        auto &ars = template_alignment_info->get_aligned_ref_seq();
        auto &ats = template_alignment_info->get_aligned_seq();
#ifdef DEBUG
        std::cout << ats << "\n"
                  << astring << "\n"
                  << ars << std::endl
                  << std::endl;
#endif

        for (size_t i = 0; i < astring.size(); ++i)
        {
            if (astring[i] != '-' && astring[i] != '^')
            { // MATCH and mismatch
                ++t;
                ++r;
            }
            else
            { // ins and del
                if (ars[i] == '-')
                {
                    ++t;
                }
                else
                {
                    ++r;
                }
            }
            if (as == -1 && t == s + 1)
            {
                as = i;
                if (!is_rev)
                {
                    rstart = rs + r - 1;
                }
                else
                {
                    rend = re - r + 1;
                }
            }
            if (ae == -1 && t == e + 1)
            {
                ae = i;
                if (!is_rev)
                {
                    rend = rs + r - 1;
                }
                else
                {
                    rstart = re - r + 1;
                }
            }
        }

        int32_t lena = ae - as + 1;
        std::string alignt = ats.substr(as, lena);
        std::string alignr = ars.substr(as, lena);
        std::string aligns = astring.substr(as, lena);
        int32_t prefixl, sufixl;
        if (!is_prev)
        {
            prefixl = asp - 1;
            sufixl = plen - aep;
        }
        else
        {
            sufixl = asp - 1;
            prefixl = plen - aep;
            revcom(alignt);
            revcom(alignr);
            reverse(aligns);
        }
        std::string prefixp = seqp.substr(0, prefixl);
        std::string sufixp = seqp.substr(plen - sufixl);

        // convert end "*/-/^" to "#"
        int32_t rmx = -1;
        int32_t rml = 0;
        for (int32_t i = 0; i < lena; i++)
        {
            if (aligns[i] == '|')
                break;
            if (aligns[i] == '^')
            {
                alignr[i] = 'N';
            }
            else if (aligns[i] == '-')
            { // end del: *---||||
                rml++;
                rmx = rmx == -1 ? i : rmx;
            }
            aligns[i] = '#';
        }
        if (rmx != -1)
        { // rm - in alignt
            remove(alignt, rmx, rml);
            remove(aligns, rmx, rml);
            remove(alignr, 0, rml);
        }

        rmx = -1;
        rml = 0;
        for (int32_t i = lena - 1; i > 0; i--)
        {
            if (aligns[i] == '|')
                break;
            if (aligns[i] == '^')
            {
                alignr[i] = 'N';
            }
            else if (aligns[i] == '-')
            {
                rml++;
                rmx = i;
            }
            aligns[i] = '#';
        }
        if (rmx != -1)
        { // rm - in alignt
            remove(alignt, rmx, rml);
            remove(aligns, rmx, rml);
            remove(alignr, lena - rml, rml);
        }

        //
        std::string prefixr(prefixl, 'N');
        std::string prefixa(prefixl, '#');
        std::string sufixr(sufixl, 'N');
        std::string sufixa(sufixl, '#');

        std::string pseq = prefixp + alignt + sufixp;
        std::string astr = prefixa + aligns + sufixa;
        // rseq need to revcom
        revcom(alignr);
        std::string rseq = sufixr + alignr + prefixr;

#ifdef DEBUG
        std::cout << prefixp << "," << sufixp << ";" << prefixr << "," << sufixr << std::endl;
        std::cout << pseq << "\n"
                  << astr << "\n"
                  << rseq << std::endl;
#endif

        double tm, dG, dS, dH;
        auto iter = alignstr_filter_map.find(astr);
        if (iter == alignstr_filter_map.end())
        {
            // caculate tm
            tm_calc->alignment_tm(astr, pseq, rseq, &tm, &dG, &dS, &dH);
#ifdef DEBUG
            std::cout << chr << "\t" << rstart << "\t" << strand << "\t" << align << "\t" << tm << "\t" << dG << "\t" << dS << "\t" << dH << std::endl;
#endif
            if(tm < min_bound_tm+5){
				alignstr_filter_map.insert(std::pair<std::string, double>{astr, tm});
            }
        }
        else
        {
            tm = iter->second;
#ifdef DEBUG
            std::cout << "iter:\t" << pseq << "\t" << astr << "\t" << tm << std::endl;
#endif
        }
#ifdef DEBUG
        std::cout << chr << "\t" << rstart << "\t" << strand << "\t" << astr << "\t" << tm << std::endl;
#endif

        if (tm < min_bound_tm)
        { // filtered out
            return true;
        }
        int64_t pos = strand == '+' ? rstart - asp + 1 : rend + asp - 1;
        primer_hits_interval.set(pos, chr, &seqp, astr, tm, strand);

        return false;
    }

}
