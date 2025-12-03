#pragma READ_WRITE
using namespace std;
#include <string>


void throw_error(const string msg);
double read_from_file (string input_name, const string param);
string read_from_file_str (string input_name, const string param, const string def_val);
void make_dir(std::string dir_name);

void read_in_out(string &in, string &out, int argc, char* argv[]);

