#include "primer.h"
#include "alignment_info.h"
#include "primer_utils.h"
#include "options.h"
#include "CLI/CLI.hpp"
#include <tbb/parallel_for.h>
#include <tbb/pipeline.h>
#include <tbb/task_arena.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <glog/logging.h>
#include "tm_calculator.h"

static void pt_cmdline(CLI::App &app, int32_t argc, char *argv[], pt::options &opt)
{
    app.add_option(
           "-i,--primer-list",
           opt.primer_list_file,
           "Input primer list file")
        ->required();

    app.add_option(
           "-f,--template-fasta",
           opt.template_fasta_file,
           "Input FASTA of template file.")
        ->required();

    app.add_option(
           "-a,--bam",
           opt.alignment_bam,
           "Blast alignment info of template sequence in BAM format.")
        ->required();

    app.add_option(
           "-p,--path",
           opt.path,
           "The path to the thermodynamic parameter files.")
        ->required();

    app.add_option(
           "-o,--output",
           opt.output_file,
           "Output file name with primer alignment info.")
        ->required();

    app.add_option(
        "-l,--min-overlap-3-end",
        opt.min_overlap_3_end,
        "Minimum overlap length between primer and aligned template blocks. Default [12].");

    app.add_option(
        "-u,--max-non-overlap-3-end",
        opt.max_non_overlap_3_end,
        "Maximum non-overlap length between primer and aligned template blocks. Default [1]");

    app.add_option(
        "-b,--min_bound_tm",
        opt.min_bound_tm,
        "Minimum tm to bound and amplify non-target template. Default [45]");

    app.add_option(
        "-m,--monovalent_conc",
        opt.mv,
        "Concentration of monovalent cations in mM. Default [50]");

    app.add_option(
        "-v,--divalent_conc",
        opt.dv,
        "Concentration of divalent cations in mM. Default [1.5]");

    app.add_option(
        "-n,--dNTP_conc",
        opt.dntp,
        "Concentration of deoxynycleotide triphosphate in mM. Default [0.6]");

    app.add_option(
        "-d,--dna_conc",
        opt.dna,
        "Concentration of DNA strands in nM. Default [50]");

    app.add_option(
        "-T,--temp",
        opt.temp,
        "Temperature at which duplex is calculated. Default [37]");

    app.add_option(
        "-t,--threads",
        opt.nthreads,
        "The number of threads.");

    app.add_option(
        "-s,--max-hits-to-split",
        opt.max_hits_to_split,
        "If the total hits of all primers exceeds this threshold, split the primer set into pieces to decrease memory consumption. Default [10000000]");

    app.add_option(
        "-M,--max-primer-bounds",
        opt.max_hits_allowed,
        "If the total hits of one primer exceeds this threshold, discard the primer. -1 means no filtering. Default [-1]");

    app.add_flag(
        "--probe-mode",
        opt.probe_mode,
        "Probe mode.");

    app.add_flag(
        "-x,--disable-split-mode",
        opt.disable_split_mode,
        "Disable spliting primer set, faster but more memory consumed.");
}

// for pipeline

// read bam for one template
struct read_bam_filter
{
    const pt::options *opt;
    std::shared_ptr<bam_hdr_t> h;
    pt::primer_template_db_ptr templates;
    pt::primer_template_db_t::iterator iter;

    read_bam_filter(const pt::options *opt, std::shared_ptr<bam_hdr_t> h, pt::primer_template_db_ptr templates);
    pt::primer_template_ptr operator()(tbb::flow_control &fc) const;
};

read_bam_filter::read_bam_filter(const pt::options *opt, std::shared_ptr<bam_hdr_t> h, pt::primer_template_db_ptr templates) : opt(opt), h(h), templates(templates), iter(templates->begin()) {}

pt::primer_template_ptr read_bam_filter::operator()(tbb::flow_control &fc) const
{
    auto current = iter;
    if (current == this->templates->end())
    {
        fc.stop();
        return nullptr;
    }
    ++const_cast<pt::primer_template_db_t::iterator &>(this->iter);

    pt::primer_template_ptr primer_template = current->second;

    // read bam file
    samFile *fp = sam_open(opt->alignment_bam.c_str(), "r");
    if (nullptr == fp)
    {
        std::cerr << "open BAM file failed.";
        exit(EXIT_FAILURE);
    }
    std::shared_ptr<samFile> auto_sam_fp{fp, [](samFile *fp)
                                         { sam_close(fp); }};

    // load bam index
    std::shared_ptr<hts_idx_t> idx{sam_index_load(fp, this->opt->alignment_bam.c_str()), [](hts_idx_t *index)
                                   { hts_idx_destroy(index); }};
    if (nullptr == idx)
    {
        std::cerr << "load BAM index failed.";
        exit(EXIT_FAILURE);
    }

    // read bam records of current primer template
    std::shared_ptr<hts_itr_t> iter{sam_itr_querys(idx.get(), this->h.get(), primer_template->get_tid().c_str()), [](hts_itr_t *itr)
                                    { hts_itr_destroy(itr); }};
    std::shared_ptr<bam1_t> b{bam_init1(), [](bam1_t *r)
                              { bam_destroy1(r); }};

    auto ret = sam_itr_next(fp, iter.get(), b.get());
    if (0 < ret)
    {
        do
        {
            std::shared_ptr<bam1_t> b1{bam_dup1(b.get()), [](bam1_t *r)
                                       { bam_destroy1(r); }};
            primer_template->add_alignment(b1);
        } while (0 < (ret = sam_itr_next(fp, iter.get(), b.get())));
    }
    return primer_template;
}

// process filter
struct process_filter
{
    const pt::options *opt;
    bam_hdr_t *h;
    pt::tm_calculator *tm_calc;

    process_filter(const pt::options *opt, bam_hdr_t *h, pt::tm_calculator *tm_calc);

    pt::primer_template_ptr operator()(pt::primer_template_ptr primer_template) const;
};

process_filter::process_filter(const pt::options *opt, bam_hdr_t *h, pt::tm_calculator *tm_calc) : opt(opt), h(h), tm_calc(tm_calc) {}

pt::primer_template_ptr process_filter::operator()(pt::primer_template_ptr primer_template) const
{
    if (nullptr == primer_template)
    {
        return nullptr;
    }
    primer_template->initialize_alignment_db(this->h);

    primer_template->get_db()->init_alignment_info_from_bam_records(primer_template->get_alignments());

    primer_template->get_db()->calculate_aligned_blocks(primer_template->get_seq());

    primer_template->clear_alignments();

    // build query index
    primer_template->get_db()->build_query_tree();

    // query primer
    primer_template->query_primers(this->opt->min_overlap_3_end, this->opt->max_non_overlap_3_end, this->opt->probe_mode);

    if (!this->opt->disable_split_mode)
    {
        int64_t total_hits = primer_template->get_total_hits();
        if (total_hits > this->opt->max_hits_to_split)
        {
            primer_template->mark_split();
            return nullptr;
        }
    }

    // calculate interval
    primer_template->calculate_primer_intervals(this->tm_calc, this->opt->min_bound_tm);

    // clear hits
    primer_template->clear_hits();
    return primer_template;
}

static void output_result(const pt::options &opt, const pt::primer_template &primer_template,  std::ostream *os, int32_t s = -1, int32_t e = -1);

struct output_filter
{
    const pt::options *opt;
    std::ostream *os;

    output_filter(const pt::options *opt, std::ostream *os);
    void operator()(pt::primer_template_ptr primer_template) const;
};


output_filter::output_filter(const pt::options *opt, std::ostream *os) : opt(opt), os(os) {}

void output_filter::operator()(pt::primer_template_ptr primer_template) const
{
    if (nullptr == primer_template)
    {
        return;
    }

    output_result(*opt, *primer_template, os);
    primer_template->clear_intervals();
}

void split_process_too_many_hits_template(const pt::options &opt, pt::primer_template_db_ptr templates, pt::tm_calculator_ptr tm_calc, std::ostream *os)
{
    for (auto it : *templates)
    {
        auto primer_template = it.second.get();
        if (!primer_template->is_marked_split())
        {
            continue;
        }
        size_t nparts = (primer_template->get_total_hits() + opt.max_hits_to_split - 1) / opt.max_hits_to_split;
        auto &primers = primer_template->get_primers();

        size_t np = primers.size();
        size_t primer_per_block = (np + nparts - 1) / nparts;

        LOG(INFO) << "[T]: " << primer_template->get_tid() << ", "
                  << "[P]: " << np << ", "
                  << "[H]: " << primer_template->get_total_hits() << ", "
                  << "[N]: " << nparts << ","
                  << "[M]: " << primer_per_block;

        auto tmc = tm_calc.get();
        for (size_t iblock = 0; iblock < nparts; ++iblock)
        {
            int32_t s = iblock * primer_per_block;
            if (static_cast<size_t>(s) >= np)
            {
                break;
            }
            int32_t e = std::min(static_cast<size_t>(s + primer_per_block), np);

            // calculate interval
            primer_template->calculate_primer_intervals(tmc, opt.min_bound_tm, s, e);

            // output
            output_result(opt, *primer_template, os, s, e);
        }
        primer_template->clear_hits();
        primer_template->clear_intervals();
    }
}

static void output_result(const pt::options &opt, const pt::primer_template &primer_template,  std::ostream *os, int32_t s, int32_t e)
{

    auto primers = primer_template.get_primers();
    // s: start primer index, e: end primer index
    s = s == -1 ? 0 : s;
    e = e == -1 ? static_cast<int32_t>(primers.size()) : e;
    if (e <= s)
    {
        LOG(FATAL) << "Block index is  not properly set.";
    }

    auto itvs = primer_template.get_primer_interval();
    assert(itvs.size() == e - s);

    for (int32_t i = s; i < e; ++i)
    {
        auto &p = *primers[i];
        auto &v = itvs[i-s];
        if (opt.max_hits_allowed != -1 && v.size() > static_cast<size_t>(opt.max_hits_allowed)){
            LOG(WARNING) << p.get_id() << " is filtered due to too many (" << v.size() << ") bounds!";
            continue;
        }
        for (auto &itv : v)
        {
            if (itv.is_retained())
            {
                *os << p.get_id() << "\t"
                    << itv.to_string() << "\n";
            }
        }
    }
}


int main(int argc, char *argv[])
{
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = 1;

    // get options from command line
    CLI::App app{"This tools is used to get primer alignment info from alignment of template."};
    pt::options opt;
    pt_cmdline(app, argc, argv, opt);
    CLI11_PARSE(app, argc, argv);

    LOG(INFO) << "Program start";
    LOG(INFO) << "Load Primer list file.";
    // read in primer list
    pt::primer_template_db_ptr templates = pt::primer_utils::get_primer_template_from_file(opt.primer_list_file);
    LOG(INFO) << templates->size() << " templates loaded.";
    LOG(INFO) << "Done load Primer list file.";

    LOG(INFO) << "Load template fasta file.";
    // read template fasta
    pt::primer_utils::get_template_sequence_from_fasta_file(opt.template_fasta_file, templates);
    LOG(INFO) << "Done load template fasta file.";

    // get bam header
    auto header = pt::primer_utils::read_bam_header(opt.alignment_bam);

    // output stream
    std::ofstream ofs(opt.output_file.c_str());

    LOG(INFO) << "Initialize TM calculator.";
    // init TM calculator
    pt::tm_calculator_ptr tm_calc = std::make_shared<pt::tm_calculator>(opt);
    LOG(INFO) << "Done initialize TM calculator.";

    // pipeline of parse primer info
    /*
        1. emit an template from templates db
        2. read alignment info from bam
        3. build alignment info from bam record list
        4. build interval tree for fast query of primer end pos
        5. query overlapping alignments for each primer in parallel
        6. output
    */
    LOG(INFO) << "Start primer process.";
    tbb::task_arena limited_arena(opt.nthreads);

    tbb::filter_t<void, pt::primer_template_ptr> read_bam(tbb::filter::serial_in_order, read_bam_filter(&opt, header, templates));
    tbb::filter_t<pt::primer_template_ptr, pt::primer_template_ptr> process(tbb::filter::parallel, process_filter(&opt, header.get(), tm_calc.get()));
    tbb::filter_t<pt::primer_template_ptr, void> output(tbb::filter::serial_out_of_order, output_filter(&opt, &ofs));
    tbb::filter_t<void, void> pipeline = read_bam & process & output;

    limited_arena.execute([&]()
                          { tbb::parallel_pipeline(opt.nthreads, pipeline); });

    // process template with large amount of hits by dividing primer into pieces
    if (!opt.disable_split_mode)
    {
        split_process_too_many_hits_template(opt, templates, tm_calc, &ofs);
    }
    LOG(INFO) << "End primer process.";

    // close ofile
    ofs.close();

    LOG(INFO) << "Program done.";
    google::ShutdownGoogleLogging();
    return 0;
}
