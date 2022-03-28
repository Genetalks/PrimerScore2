#ifndef PRIMER_UTILS_H
#define PRIMER_UTILS_H

#include <string>
#include <memory>
#include "primer.h"
#include "htslib/sam.h"

namespace pt
{

    class primer_utils
    {
    public:
        static primer_template_db_ptr get_primer_template_from_file(const std::string &primer_list_file);
        static std::shared_ptr<bam_hdr_t> read_bam_header(const std::string &bam);
        static void get_template_sequence_from_fasta_file(const std::string &fasta_file, primer_template_db_ptr templates);
    };

};

#endif
