#ifndef TOPO_WRITE_H
#define TOPO_WRITE_H

#include <iostream>
#include <fstream>

// Regular expression http://www.cplusplus.com/reference/regex/regex_replace/
#include <string>
#include <regex>
#include <iterator>

#include <sys/stat.h>

// Directory related part adapted from https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
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

void tokenize(std::string const &str, std::string& delim,
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

// Taken from https://stackoverflow.com/questions/25201131/writing-csv-files-from-c
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


void ReadNumbers(const std::string& filename, char sep, 
                 std::vector<std::vector<int> > & output);

void ReadNumbers(const std::string& filename, char sep, 
                 std::vector<std::vector<double> > & output);

// Read lines from filename, and add the lines separated by sep into the output.
void ReadNames(const std::string& filename, const char * sep, 
                             std::vector<std::vector<std::string> > & output);

#endif
