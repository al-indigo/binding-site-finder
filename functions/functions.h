#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include "../ahoc/AhoCorasickPlus.h"
#include "../structs/pwm.h"
#include "../structs/chromovector.h"

void init_ahoc(AhoCorasickPlus& atm, 
               Pwm& matrix, 
               size_t patterns_allowed, 
               size_t& total_words, 
               optimization_type type = classic);

void ahoc_search(size_t sequence_id, 
                 size_t part_id,
                 AhoCorasickPlus& atm,
                 Pwm& matrix,
                 ChromoVector& sequences, 
                 const std::string& result_folder,
                 const std::string& result_filename,
                 size_t& total_words,
                 std::vector<std::vector<std::string> >& files_to_merge);

void merge_files(int id,
                 const std::string& result_folder,
                 const std::string& result_filename,
                 std::vector<std::vector<std::string> >& files_to_merge,
                 std::vector<std::string>& merged_files);

//NOTE: istream is here for purpose: we can use it with files and memory.
void recalc_scores_p_values(std::istream& fin, 
                            std::vector<size_t>& positions_to_read_again, 
                            size_t mem_allowed, 
                            Pwm& matrix, 
                            ChromoVector& sequences, 
                            size_t sequence_id, 
                            std::vector<double>& pvaluesFw, 
                            std::vector<double>& pvaluesRev, 
                            std::vector<double>& scoresFw, 
                            std::vector<double>& scoresRev);

void recalc_scores_p_values(FILE * fin, 
                            std::vector<size_t>& positions_to_read_again, 
                            size_t mem_allowed, 
                            Pwm& matrix, 
                            ChromoVector& sequences, 
                            size_t sequence_id, 
                            std::vector<double>& pvaluesFw, 
                            std::vector<double>& pvaluesRev, 
                            std::vector<double>& scoresFw, 
                            std::vector<double>& scoresRev);

void format_bed(std::ostream& fout,
                std::string& chromoname,
                std::vector<size_t>& positions_to_read_again,
                std::vector<double>& pvaluesFw,
                std::vector<double>& pvaluesRev,
                std::vector<double>& scoresFw,
                std::vector<double>& scoresRev);

void format_bed(FILE * fout,
                std::string& chromoname,
                std::vector<size_t>& positions_to_read_again,
                std::vector<double>& pvaluesFw,
                std::vector<double>& pvaluesRev,
                std::vector<double>& scoresFw,
                std::vector<double>& scoresRev);

void write_status(double percent, std::string status_folder, std::string status_filename, const char * explain, const char * result_file_path);


#endif