#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <filesystem>
#include "read_write.h"

using namespace std;
namespace fs = std::filesystem;


void throw_error(const string msg) {
  cout << msg << endl;
  exit (EXIT_FAILURE);
}
bool is_number(const string& s) {
  if (s.empty())
    return false;
  bool trigger = false;
  if (s[0] == '.') {
    return false;
  }
  for (int i = 0; i < s.size(); i ++) {
    char c = s[i];
    if (c == '.' && !trigger) {
      trigger = !trigger;
      continue;
    } else if (c == '.' && trigger) {
      return false;
    }
    if (i == 0 && c == '-') {
      continue;
    }
    if (!std::isdigit(c))
    return false;
  }
  return true;
}
double read_from_file (string input_name, const string param) {
  ifstream infile(input_name);
  string str;
  stringstream msg;
  if (infile.is_open()) {
    while (getline(infile, str)) {
      // cout << iter << ": " << str << "\n";
      istringstream iss(str);
      vector<string> words{istream_iterator<string>{iss}, istream_iterator<string>{}};
      if(words.size() < 1) continue;
      // cout << words[0] << "\n";
      if (words[0] == param) {
        if(is_number(words[1])) {
          return atof(words[1].c_str());
        } else {
          msg << "ERROR: Converting string to number for parameter '" << param << "' in the input file.";
          throw_error(msg.str());
        }
        break;
      }
    }
  } else {
    throw_error("ERROR: Cannot open the input file.");
  }
  msg << "ERROR: Cannot find the given parameter '" << param << "' in the input file.";
  throw_error(msg.str());
  return 0;
}
string read_from_file_str (string input_name, const string param, const string def_val) {
  ifstream infile(input_name);
  string str;
  stringstream msg;
  if (infile.is_open()) {
    while (getline(infile, str)) {
      istringstream iss(str);
      vector<string> words{istream_iterator<string>{iss}, istream_iterator<string>{}};
      if(words.size() < 1) continue;
      if (words[0] == param) {
        if(!is_number(words[1])) {
          return words[1].c_str();
        } else {
          msg << "ERROR: param '" << param << "' is not a string in the input file.";
          throw_error(msg.str());
        }
        break;
      }
    }
  } else {
    throw_error("ERROR: Cannot open the input file.");
  }
  return def_val;
}

void make_dir(std::string dir_name) {
  if (!fs::is_directory(dir_name) || !fs::exists(dir_name)) { // Check if folder exists
    fs::create_directory(dir_name); // create folder
  } 
}

void read_in_out(std::string &in, std::string &out, int argc, char* argv[]) {
  bool found_input = false, found_output = false;
  for (int i=0; i<argc; ++i) {
		if (typeid(argv[i]) == typeid(char*)) {
			if (strncmp(argv[i], "-i", 3) == 0) {
				if (i+1 < argc) {
					in = argv[i+1];
          found_input = true;
        }
				else {
          throw_error("ERROR: No input file given.");
				}
			}
      if (strncmp(argv[i], "-o", 3) == 0) {
				if (i+1 < argc) {
					out = argv[i+1];
          found_output = true;
        }
				else {
          throw_error("ERROR: No output path given.");
				}
			}
		}
	}
  if (!found_input) {
    throw_error("ERROR: No input file given.");
  } 
  // else {
  //   cout << "INPUT: " << in << "\n";
  // }
  if (!found_output) {
    out = "output";
  }
}
