#include <string>
#include <fstream>

void write_status(double percent, std::string status_folder, std::string status_filename, const char * explain, const char * result_file_path) {
  std::ofstream status_out((status_folder + status_filename).c_str(), std::fstream::out);
  status_out << "{\"percent done\": " << (int)percent << ", \"result\": [\"" << result_file_path << "\"], \"explain\": \"" << explain << "\"}";
  status_out.close();
}