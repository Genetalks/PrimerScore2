#include "primer_utils.h"
#include "htslib/hts.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <algorithm>

pt::primer_template_db_ptr pt::primer_utils::get_primer_template_from_file(const std::string &primer_list_file)
{
    // primer file is a 5-colomn tsv file
    // eg: primer_id    template_id     start   end     strand
    
    primer_template_db_ptr db = std::make_shared<primer_template_db_t>();

    std::ifstream ifs(primer_list_file.c_str());
    if (!ifs.good()){
        std::cerr << "open primer list file failed.";
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::vector<std::string> columns;
    const std::regex re("\\s+");

    while (std::getline(ifs, line)){
        columns.clear();
        std::copy(std::sregex_token_iterator(line.begin(), line.end(), re, -1), std::sregex_token_iterator(), std::back_inserter(columns));
        if (columns.size() != 5){
            std::cerr << "primer list file format error, 5 columns required, " << columns.size() << " columns detected.\n";
            std::cerr << "line: " << line << "\n";
            exit(EXIT_FAILURE);
        }

        // query tamplate in db
        const std::string &template_id = columns[1];
        const std::string &primer_id = columns[0];
        int32_t start = std::stoi(columns[2]);
        int32_t end = std::stoi(columns[3]);
        bool is_rev = columns[4][0] == '-' ? true : false;

        primer_ptr p = std::make_shared<primer>(primer_id, template_id, start, end, is_rev);

        auto itr = db->find(template_id);
        if (itr == db->end()){
            primer_template_ptr tmp = std::make_shared<primer_template>(template_id);
            tmp->add_primer(p);
            db->insert(std::pair<std::string, primer_template_ptr>{template_id, tmp});
        }else{
            itr->second->add_primer(p);
        }
    }
    ifs.close();

    return db;
}

std::shared_ptr<bam_hdr_t> pt::primer_utils::read_bam_header(const std::string &bam){
    samFile *fp = sam_open(bam.c_str(), "r");
    if (nullptr == fp){
        std::cerr << "open Bam file failed.";
        exit(EXIT_FAILURE);
    }

    std::shared_ptr<bam_hdr_t> header{sam_hdr_read(fp), [](bam_hdr_t *h){bam_hdr_destroy(h);}};
    if (nullptr == header){
        std::cerr << "Read Bam header failed.";
        exit(EXIT_FAILURE);
    }
    return header;
}

void pt::primer_utils::get_template_sequence_from_fasta_file(const std::string &fasta_file, pt::primer_template_db_ptr templates){
    std::ifstream ifs(fasta_file.c_str());
    if (!ifs.good()){
        std::cerr << "Open template fasta failed.\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::string tid;
    std::string seq;
    std::string last_tid{""};
    int64_t line_no = 0;
    while(std::getline(ifs, line)){
        ++line_no;
        if (line[0] == '>'){
            auto pos_of_ws = line.find(' ');
            if (pos_of_ws != std::string::npos){
                tid = line.substr(1, pos_of_ws);
            }else{
                tid = line.substr(1);
            }
            if (!last_tid.empty()){
                auto exist_such_tid = templates->find(last_tid);
                if (exist_such_tid != templates->end()){
                    exist_such_tid->second->set_seq(seq);
                }
            }
            seq.clear();
            last_tid = tid;
        }else{
            // check ACGT
            bool has_non_acgt_char = false;
            for (auto &c : line){
                switch (c){
                    case 'A':
                        break;
                    case 'C':
                        break;
                    case 'G':
                        break;
                    case 'T':
                        break;
                    case 'N':
                        break;
                    default:
                        has_non_acgt_char = true;
                }
            }

            if (has_non_acgt_char){
                std::cerr << "line " << line_no << " contains non-ACGTN characters.";
                exit(EXIT_FAILURE);
            }

            seq.append(line);
        }
    }
    // append last
    if (!last_tid.empty()){
        auto exist_such_tid = templates->find(last_tid);
        if (exist_such_tid != templates->end()){
            exist_such_tid->second->set_seq(seq);
        }
    }
    ifs.close();

    // validate template completeness
    for (auto &t : *templates){
        if (t.second->get_seq().empty()){
            std::cerr << "No sequence for template :" << t.second->get_tid();
            exit(EXIT_FAILURE);
        }
    }
}
