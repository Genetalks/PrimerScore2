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

static void pt_cmdline(CLI::App &app, int32_t argc, char *argv[], pt::options &opt){
    app.add_option(
	"-i,--primer-list", 
	opt.primer_list_file, 
	"Input primer list file"
    )->required();

    app.add_option(
	"-f,--template-fasta", 
	opt.template_fasta_file, 
	"Input FASTA of template file."
    )->required();
    
    app.add_option(
	"-a,--bam", 
	opt.alignment_bam, 
	"Blast alignment info of template sequence in BAM format."
    )->required();

    app.add_option(
	"-o,--output", 
	opt.output_file, 
	"Output file name with primer alignment info."
    )->required();

    app.add_option(
	"-l,--min-overlap-3-end", 
	opt.min_overlap_3_end, 
	"Minimum overlap length between primer and aligned template blocks. Default [12]."
    );

    app.add_option(
	"-u,--max-non-overlap-3-end", 
	opt.max_non_overlap_3_end, 
	"Maximum non-overlap length between primer and aligned template blocks. Default [1]"
    );

	app.add_option(
	"-b,--min_bound_tm", 
	opt.min_bound_tm, 
	"Minimum tm to bound and amplify non-target template. Default [45]"
    );

    app.add_option(
	"-m,--monovalent_conc",
	opt.mv, 
	"concentration of monovalent cations in mM. Default [50]"
    );

    app.add_option(
	"-v,--divalent_conc", 
	opt.dv,
	"concentration of divalent cations in mM. Default [1.5]"
    );

    app.add_option(
	"-n,--dNTP_conc", 
	opt.dntp,
	"concentration of deoxynycleotide triphosphate in mM. Default [0.6]"
    );

    app.add_option(
	"-d,--dna_conc", 
	opt.dna,
	"concentration of DNA strands in nM. Default [50]"
    );

    app.add_option(
	"-T,--temp", 
	opt.temp,
	"temperature at which duplex is calculated. Default [37]"
    );

	app.add_option(
	"-p,--path", 
	opt.path,
	"the path to the thermodynamic parameter files."
    )->required();

    app.add_option(
	"-t,--threads", 
	opt.nthreads, 
	"The number of threads."
    );

    app.add_flag(
	"--probe-mode", 
	opt.probe_mode, 
	"Probe mode."
    );
}

// for pipeline
struct emit_template_filter {
    pt::primer_template_db_ptr templates;
    pt::primer_template_db_t::iterator iter;

    emit_template_filter(pt::primer_template_db_ptr templates);
    pt::primer_template_ptr operator()(tbb::flow_control &fc) const ;
};

emit_template_filter::emit_template_filter(pt::primer_template_db_ptr templates):templates(templates), iter(templates->begin()){}
pt::primer_template_ptr emit_template_filter::operator()(tbb::flow_control &fc) const{
  auto current = iter;
  if (current == this->templates->end()){
      fc.stop();
      return nullptr;
  }
  ++const_cast<pt::primer_template_db_t::iterator &>(this->iter);

  return current->second;
}


// read bam for one template
struct read_bam_filter {
    const pt::options *opt;
    std::shared_ptr<bam_hdr_t> h;

    read_bam_filter(const pt::options *opt, std::shared_ptr<bam_hdr_t> h);
    pt::primer_template_ptr operator()(pt::primer_template_ptr primer_template) const;
};

read_bam_filter::read_bam_filter(const pt::options *opt, std::shared_ptr<bam_hdr_t> h):opt(opt), h(h) {}

pt::primer_template_ptr read_bam_filter::operator()(pt::primer_template_ptr primer_template) const {
    if (nullptr == primer_template){
        return nullptr;
    }

    // read bam file
    samFile *fp = sam_open(opt->alignment_bam.c_str(), "r");
    if (nullptr == fp){
        std::cerr << "open BAM file failed.";
        exit(EXIT_FAILURE);
    }
    std::shared_ptr<samFile> auto_sam_fp{fp, [](samFile *fp){sam_close(fp);}};
    
    // load bam index
    std::shared_ptr<hts_idx_t> idx{sam_index_load(fp, this->opt->alignment_bam.c_str()), [](hts_idx_t *index){hts_idx_destroy(index);}};
    if (nullptr == idx){
        std::cerr << "load BAM index failed.";
        exit(EXIT_FAILURE);
    }

    // read bam records of current primer template
    std::shared_ptr<hts_itr_t> iter{sam_itr_querys(idx.get(), this->h.get(), primer_template->get_tid().c_str()), [](hts_itr_t *itr){hts_itr_destroy(itr); }};
    std::shared_ptr<bam1_t> b{bam_init1(), [](bam1_t *r){bam_destroy1(r);}};

    auto ret = sam_itr_next(fp, iter.get(), b.get());
    if (0 < ret){
        do {
            std::shared_ptr<bam1_t> b1{bam_dup1(b.get()), [](bam1_t *r){bam_destroy1(r);}};
            primer_template->add_alignment(b1);
        } while (0 < (ret = sam_itr_next(fp, iter.get(), b.get())));
    }
    return primer_template;
}

// build alignment info
struct build_alignment_info_filter{
    bam_hdr_t *h;

    build_alignment_info_filter(bam_hdr_t *h);

    pt::primer_template_ptr operator()(pt::primer_template_ptr primer_template) const;
};

build_alignment_info_filter::build_alignment_info_filter(bam_hdr_t *h):h(h) {}

pt::primer_template_ptr build_alignment_info_filter::operator()(pt::primer_template_ptr primer_template) const{
    if (nullptr == primer_template){
        return nullptr;
    }
    primer_template->initialize_alignment_db(this->h);

    primer_template->get_db()->init_alignment_info_from_bam_records(primer_template->get_alignments());

    primer_template->get_db()->calculate_aligned_blocks(primer_template->get_seq());

    primer_template->clear_alignments();

    return primer_template;
}

// build query interval tree index
struct index_filter{
    const pt::options* opt;

    index_filter(const pt::options *opt);
    pt::primer_template_ptr operator()(pt::primer_template_ptr primer_template) const;
};

index_filter::index_filter(const pt::options *opt):opt(opt){}

pt::primer_template_ptr index_filter::operator()(pt::primer_template_ptr primer_template) const{
    if (nullptr == primer_template){
        return nullptr;
    }

    primer_template->get_db()->build_query_tree();
    return primer_template;
}

// query primer
struct query_filter {
    const pt::options *opt;

    query_filter(const pt::options *opt);
    pt::primer_template_ptr operator()(pt::primer_template_ptr primer_template) const;
};

query_filter::query_filter(const pt::options *opt):opt(opt) {}

pt::primer_template_ptr query_filter::operator()(pt::primer_template_ptr primer_template) const{
    if (nullptr == primer_template){
        return nullptr;
    }
    primer_template->query_primers(this->opt->min_overlap_3_end, this->opt->max_non_overlap_3_end, this->opt->probe_mode);
    return primer_template;
}



struct cal_primer_interval_filter {
    const pt::options *opt;
	cal_primer_interval_filter(const pt::options *opt);
    pt::primer_template_ptr operator()(pt::primer_template_ptr primer_template) const;
};

cal_primer_interval_filter::cal_primer_interval_filter(const pt::options *opt):opt(opt) {}

pt::primer_template_ptr cal_primer_interval_filter::operator()(pt::primer_template_ptr primer_template) const {
    if (nullptr == primer_template){
        return nullptr;
    }
    primer_template->calculate_primer_intervals(this->opt->min_bound_tm, this->opt->mv, this->opt->dv, this->opt->dntp, this->opt->dna, this->opt->temp, this->opt->path);
    return primer_template;
}

struct output_filter {
    const pt::options *opt;
    std::ostream *os;

    output_filter(const pt::options *opt, std::ostream *os);
    void operator()(pt::primer_template_ptr primer_template) const;
};

output_filter::output_filter(const pt::options *opt, std::ostream *os):opt(opt), os(os) {}

void output_filter::operator()(pt::primer_template_ptr primer_template) const{
    if (nullptr == primer_template){
        return;
    }
    
    auto primers = primer_template->get_primers();
    auto hits = primer_template->get_primer_interval();

    assert(primers.size() == hits.size());

    for (size_t i = 0; i < primers.size(); ++i){
        auto &p = primers[i];
        for (auto &itv : hits[i]){
			if(itv.is_valid){
	            *os << p->get_id() << "\t"
		            << itv.to_string() << "\n";
			}
        }
    }
}

int main(int argc, char *argv[]){
    // get options from command line
    CLI::App app{"This tools is used to get primer alignment info from alignment of template."};
    pt::options opt;
    pt_cmdline(app, argc, argv, opt);    
    CLI11_PARSE(app, argc, argv);

    // read in primer list 
    pt::primer_template_db_ptr templates = pt::primer_utils::get_primer_template_from_file(opt.primer_list_file);

    // read template fasta 
    pt::primer_utils::get_template_sequence_from_fasta_file(opt.template_fasta_file, templates);

    // get bam header
    auto header = pt::primer_utils::read_bam_header(opt.alignment_bam);

    // output stream
    std::ofstream ofs(opt.output_file.c_str());

    // pipeline of parse primer info 
    /*
        1. emit an template from templates db
        2. read alignment info from bam 
        3. build alignment info from bam record list 
        4. build interval tree for fast query of primer end pos
        5. query overlapping alignments for each primer in parallel
        6. output
    */
    tbb::task_arena limited_arena(opt.nthreads);

    tbb::filter_t<void,  pt::primer_template_ptr> emit_template(tbb::filter::serial_in_order, emit_template_filter(templates));
    tbb::filter_t<pt::primer_template_ptr, pt::primer_template_ptr> read_bam(tbb::filter::parallel, read_bam_filter(&opt, header));
    tbb::filter_t<pt::primer_template_ptr, pt::primer_template_ptr> make_alignment_info(tbb::filter::parallel, build_alignment_info_filter(header.get()));
    tbb::filter_t<pt::primer_template_ptr, pt::primer_template_ptr> make_db(tbb::filter::parallel, index_filter(&opt));
    tbb::filter_t<pt::primer_template_ptr, pt::primer_template_ptr> query_primers(tbb::filter::parallel, query_filter(&opt));
    tbb::filter_t<pt::primer_template_ptr, pt::primer_template_ptr> get_primer_intervals(tbb::filter::parallel, cal_primer_interval_filter(&opt));
    tbb::filter_t<pt::primer_template_ptr, void> output(tbb::filter::serial_out_of_order, output_filter(&opt, &ofs));
    tbb::filter_t<void, void> pipeline = emit_template & read_bam & make_alignment_info & make_db & query_primers & get_primer_intervals & output;

    limited_arena.execute([&](){
        tbb::parallel_pipeline(opt.nthreads, pipeline);
    });


    // close ofile
    ofs.close();

    return 0;
}
