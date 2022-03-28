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
    primer_template_db_ptr db = std::make_shared<primer_template_db_t>();

    std::ifstream ifs(primer_list_file.c_str());
    if (!ifs.good())
    {
        std::cerr << "open primer list file failed.";
        exit(EXIT_FAILURE);
    }

    std::string line;
    const std::regex id_pattern{"^((\\S+)-([FR])-(\\d+)_(\\d+)_([FR])_(\\d+))\\t(\\S+)"};
    const std::regex id_pattern2{"^(\\S+)\\t(\\S+)"};
    while (std::getline(ifs, line))
    {
        if (line[0] == '#')
            continue;
        smatch m;
        std::string primer_id, template_id, seq;
        int32_t start, end;
        bool is_rev;
        // eg:
        // chr10_94988855-U-R-5_32_F_0     GCTCTCTGTGTTTGCTATTTTCAGGAAA    28      62.48   0.393   34.00   NA      NA      3       -8.5    NA      0A3,7T4,15T3
        // chr10_94988855-U-R-5_32_F_2     TCTCTGTGTTTGCTATTTTCAGGAAA      26      59.45   0.346   34.00   NA      NA      3       -8.5    NA      0A3,7T4,15T3
        // chr10_94988855-U-R-5_32_F_4     TCTGTGTTTGCTATTTTCAGGAAA        24      57.29   0.333   34.00   NA      NA      3       -8.5    NA      0A3,7T4,15T3
        // chr10_94988855-U-R-23_50_F_0    TACGCATGAGGAGTAACTGCTCTCTGTG    28      65.89   0.500   46.69   NA      NA      0       -6.7    NA      NA
        bool is_match = regex_search(line, m, id_pattern);
        if (is_match)
        {
            // query tamplate in db
            primer_id = m[1];
            template_id = m[2];
            seq = m[8];
            int32_t startt = std::stoi(m[4]);
            int32_t endt = std::stoi(m[5]);
            int32_t off = std::stoi(m[7]);
            int32_t len = endt - startt + 1 - off;
            if (m[3] == m[6])
            { // +
                is_rev = 0;
                end = endt;
                start = end - len + 1;
            }
            else
            {
                is_rev = 1;
                start = startt;
                end = start + len - 1;
            }
        }
        else
        {
            // eg:
            // Sepci-1-R       GTATATTGTCATGACGGTCCAGGAG
            // Sepci-1-L       TGAAGTTCTCCAGGTGGCTGTTA
            // Sepci-2-L       CCCCCAAAGTTCGTGTGTTAGA
            // Sepci-2-R       ATTAGGAACCCACTCCCAAGATAAC
            regex_search(line, m, id_pattern2);
            primer_id = m[1];
            template_id = m[1];
            seq = m[2];
            start = 0;
            end = seq.size() - 1;
            is_rev = 0;
        }
        // cout<<primer_id<<","<<template_id<<","<<start<<","<<end<<","<<is_rev<<","<<seq<<endl;

        primer_ptr p = std::make_shared<primer>(primer_id, template_id, start, end, is_rev, seq);

        auto itr = db->find(template_id);
        if (itr == db->end())
        {
            primer_template_ptr tmp = std::make_shared<primer_template>(template_id);
            tmp->add_primer(p);
            db->insert(std::pair<std::string, primer_template_ptr>{template_id, tmp});
        }
        else
        {
            itr->second->add_primer(p);
        }
    }
    ifs.close();

    return db;
}

std::shared_ptr<bam_hdr_t> pt::primer_utils::read_bam_header(const std::string &bam)
{
    samFile *fp = sam_open(bam.c_str(), "r");
    if (nullptr == fp)
    {
        std::cerr << "open Bam file failed.";
        exit(EXIT_FAILURE);
    }

    std::shared_ptr<bam_hdr_t> header{sam_hdr_read(fp), [](bam_hdr_t *h)
                                      { bam_hdr_destroy(h); }};
    if (nullptr == header)
    {
        std::cerr << "Read Bam header failed.";
        exit(EXIT_FAILURE);
    }
    return header;
}

void pt::primer_utils::get_template_sequence_from_fasta_file(const std::string &fasta_file, pt::primer_template_db_ptr templates)
{
    std::ifstream ifs(fasta_file.c_str());
    if (!ifs.good())
    {
        std::cerr << "Open template fasta failed.\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::string tid;
    std::string seq;
    std::string last_tid{""};
    int64_t line_no = 0;
    while (std::getline(ifs, line))
    {
        ++line_no;
        if (line[0] == '>')
        {
            auto pos_of_ws = line.find('\t');
            if (pos_of_ws != std::string::npos)
            {
                tid = line.substr(1, pos_of_ws - 1);
            }
            else
            {
                tid = line.substr(1);
            }
            // std::cout<<tid<<std::endl;
            if (!last_tid.empty())
            {
                auto exist_such_tid = templates->find(last_tid);
                if (exist_such_tid != templates->end())
                {
                    exist_such_tid->second->set_seq(seq);
                    // std::cout << last_tid << ":"<< seq<<"\n"<<std::endl;
                }
            }
            seq.clear();
            last_tid = tid;
        }
        else
        {
            // check ACGT
            bool has_non_acgt_char = false;
            for (auto &c : line)
            {
                switch (c)
                {
                case 'A':
                    break;
                case 'C':
                    break;
                case 'G':
                    break;
                case 'T':
                    break;
                case 'a':
                    c -= 32;
                    break;
                case 'c':
                    c -= 32;
                    break;
                case 'g':
                    c -= 32;
                    break;
                case 't':
                    c -= 32;
                    break;
                case 'N':
                    break;
                default:
                    has_non_acgt_char = true;
                }
            }

            if (has_non_acgt_char)
            {
                std::cerr << "line " << line_no << " contains non-ACGTN characters in file" << fasta_file << std::endl;
                exit(EXIT_FAILURE);
            }

            seq.append(line);
        }
    }
    // append last
    if (!last_tid.empty())
    {
        auto exist_such_tid = templates->find(last_tid);
        if (exist_such_tid != templates->end())
        {
            exist_such_tid->second->set_seq(seq);
        }
    }
    ifs.close();

    // validate template completeness
    for (auto &t : *templates)
    {
        if (t.second->get_seq().empty())
        {
            std::cerr << "No sequence for template :" << t.second->get_tid() << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}
