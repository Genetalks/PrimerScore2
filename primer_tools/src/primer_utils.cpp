#include "primer_utils.h"
#include "htslib/hts.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <algorithm>

using namespace std;
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
//    std::vector<std::string> columns;
//    const std::regex re("\\s+");

    while (std::getline(ifs, line)){
   //     columns.clear();
   //     std::copy(std::sregex_token_iterator(line.begin(), line.end(), re, -1), std::sregex_token_iterator(), std::back_inserter(columns));
   //     if (columns.size() != 5){
   //         std::cerr << "primer list file format error, 5 columns required, " << columns.size() << " columns detected.\n";
   //         std::cerr << "line: " << line << "\n";
   //         exit(EXIT_FAILURE);
   //     }

		smatch m;
		regex_search(line, m, regex("^(\\S+)-([FR])-(\\d+)_(\\d+)_([FR])_(\\d+)\\t(\\S+)"));

        // query tamplate in db
//        const std::string &template_id = columns[1];
//        const std::string &primer_id = columns[0];
//        int32_t start = std::stoi(columns[2]);
//        int32_t end = std::stoi(columns[3]);
//        bool is_rev = columns[4][0] == '-' ? true : false;

        const std::string &primer_id = m[0];
        const std::string &template_id = m[1];
        const std::string &seq = m[7];
		int32_t startt = std::stoi(m[3]);
		int32_t endt = std::stoi(m[4]);
		int32_t off = std::stoi(m[6]);
		int32_t len = endt-startt+1-off;
        int32_t start;
        int32_t end;
        bool is_rev;
		if(m[2] == m[5]){ // +
			is_rev = 0;
			end = endt;
			start = end-len+1;
		}else{
			is_rev = 1;
			start = startt;
			end = start+len-1;
		}
		cout<<primer_id<<","<<template_id<<","<<start<<","<<end<<","<<is_rev<<","<<seq<<endl;

        primer_ptr p = std::make_shared<primer>(primer_id, template_id, start, end, is_rev, seq);

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
            auto pos_of_ws = line.find('\t');
            if (pos_of_ws != std::string::npos){
                tid = line.substr(1, pos_of_ws-1);
            }else{
                tid = line.substr(1);
            }
			//std::cout<<tid<<std::endl;
            if (!last_tid.empty()){
                auto exist_such_tid = templates->find(last_tid);
                if (exist_such_tid != templates->end()){
                    exist_such_tid->second->set_seq(seq);
					//std::cout << last_tid << ":"<< seq<<"\n"<<std::endl;
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
                    case 'a':
						c-=32;
                        break;
                    case 'c':
						c-=32;
                        break;
                    case 'g':
						c-=32;
                        break;
                    case 't':
						c-=32;
                        break;
                    case 'N':
                        break;
                    default:
                        has_non_acgt_char = true;
                }
            }

            if (has_non_acgt_char){
                std::cerr << "line " << line_no << " contains non-ACGTN characters in file" << fasta_file << std::endl;
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
            std::cerr << "No sequence for template :" << t.second->get_tid() << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}
