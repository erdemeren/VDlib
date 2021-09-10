#ifndef TOPO_WRITE_H
#define TOPO_WRITE_H

#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <iterator>

#include <sys/stat.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <iostream>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <gmi_mesh.h>

bool match_regex(std::string & f, std::smatch & sm, const char* r_exp);
bool match_regex(std::string & f, std::smatch & sm, std::string &r_exp);

// Remove matched string.
std::string replace_part(std::string & f, std::string& r_exp);

// Insert suffix between matched string and the part before the matched string.
std::string insert_between(std::string& f, std::string& r_exp, std::string &ins);

void split_str_delim(std::string const &str, std::string& delim,
            std::vector<std::string> &out);

class file_list {
  public:
    std::vector<double> nbrs;
    std::vector<std::string> nbr_strs;

    bool comp_nbrs(const double a, const double b);
    void sort_nbrs(int left = -1, int right = -1);
    int partition_nbrs(int left = -1, int right = -1);

    file_list();
    ~file_list();

    void clear();
    int count_dir(const char *path);
    void list_dir(const char *path, const char *filename, const char *ext);
    // Copy constructor
    file_list(const file_list& that);
    // Copy
    file_list& operator=(const file_list& that);
};

// Write a tess file:
void gmi_write_tess(struct gmi_model* m, const char* filename);


void ReadNumbers(const std::string& filename, char sep, 
                 std::vector<std::vector<int> > & output);

void ReadNumbers(const std::string& filename, char sep, 
                 std::vector<std::vector<double> > & output);

// Read lines from filename, and add the lines separated by sep into the output.
void ReadNames(const std::string& filename, const char * sep, 
                             std::vector<std::vector<std::string> > & output);

// csvfile object:
// Taken from https://gist.github.com/rudolfovich/f250900f1a833e715260a66c87369d15
// Copyright 2021 Vladimir Shestakov

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  Vladimir Shestakov 

class csvfile;

inline static csvfile& endrow(csvfile& file);
inline static csvfile& flush(csvfile& file);

class csvfile {
    std::ofstream fs_;
    const std::string separator_;
    const std::string str_wrap_;
  public:
    csvfile(const std::string filename, const std::string separator = ";"
                                      , const std::string str_wrap = "\"")
        : fs_()
        , separator_(separator), str_wrap_(str_wrap)
    {
        fs_.exceptions(std::ios::failbit | std::ios::badbit);
        fs_.open(filename, std::ios_base::app);
    }

    ~csvfile()
    {
        flush();
        fs_.close();
    }

    void flush()
    {
        fs_.flush();
    }

    void endrow()
    {
        fs_ << std::endl;
    }

    csvfile& operator << ( csvfile& (* val)(csvfile&))
    {
        return val(*this);
    }

    csvfile& operator << (const char * val)
    {
        fs_ << str_wrap_ << val << str_wrap_ << separator_;
        return *this;
    }

    csvfile& operator << (const std::string & val)
    {
        fs_ << str_wrap_ << val << str_wrap_ << separator_;
        return *this;
    }

    template<typename T>
    csvfile& operator << (const T& val)
    {
        fs_ << val << separator_;
        return *this;
    }
};


inline static csvfile& endrow(csvfile& file)
{
    file.endrow();
    return file;
}

inline static csvfile& flush(csvfile& file)
{
    file.flush();
    return file;
}


#endif
