#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <cstdint>
#include <thread>

namespace pt {

struct options {
    std::string template_fasta_file;
    std::string alignment_bam; // bam file contains 
    std::string primer_list_file;   // pid<tab>tid<tab>start<tab>end<tab>strand
    std::string output_file;
	int32_t max_non_overlap_3_end = 1;
	int32_t min_overlap_3_end = 12;
    int32_t nthreads = std::thread::hardware_concurrency();
    bool probe_mode = false;
	double mv = 50;
	double dv = 1.5;
	double dntp = 0.6;
	double dna = 50;
	double temp = 37;
    std::string path; //primer3 config path
};


};

#endif

