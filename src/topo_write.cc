#include <iostream> 
#include <fstream>
#include <string>
#include <vector>

#include "topo_write.h"

//https://techoverflow.net/2013/08/21/how-to-check-if-file-exists-using-stat/
bool fileExists(const char* file) {
    struct stat buf;
    return (stat(file, &buf) == 0);
}

//struct dirent {
// ino_t d_ino; /* inode number */
// off_t d_off; /* offset to the next dirent */
// unsigned short d_reclen; /* length of this record */
// unsigned char d_type; /* type of file */
// char d_name[256]; /* filename */
//

// TODO Should I be worried about overflow in readdir? Linux doesn't seem to allow
// files names beyond a size and readdir accesses only in local dir. Still may
// unlikely cause overflow errors.

// 
bool match_regex(std::string & f, std::smatch & sm, const char* r_exp) {
  std::regex e(r_exp);
  bool b = std::regex_search(f, sm, e); 
  //std::cout << " " << sm.str(0) << std::endl;
  return b;
}

bool match_regex(std::string & f, std::smatch & sm, std::string &r_exp) {
  std::regex e(r_exp);
  bool b = std::regex_search(f, sm, e); 
  //std::cout << " " << sm.str(0) << std::endl;
  return b;
}

// Given a number string with ", _ \." as separators, return the number string 
// with "." as separator.
std::string replace_sep(std::string & f) {
  std::string rep("");
  std::regex e("([0-9]+)([_\\.,])([0-9]+)");
  rep = std::regex_replace (f, e, "$1.$3");
  return rep;
}

// Remove matched string.
std::string replace_part(std::string & f, std::string& r_exp) {
  std::string rep("");
  rep = "(.*)" + r_exp;
  std::regex e(rep.c_str());
  rep = std::regex_replace (f, e, "$1");
  return rep;
}

// Insert suffix between matched string and the part before the matched string.
std::string insert_between(std::string& f, std::string& r_exp, std::string &ins) {
  std::string rep("");
  rep = "(.*)(" + r_exp + ")";
  std::regex e(rep.c_str());
  rep = std::regex_replace (f, e, "$1") + ins + std::regex_replace (f, e, "$2");
  return rep;
}

// https://www.programmersought.com/article/51335523020/
int count_nonoverlapping_sub(const std::string& str, const std::string& sub) {
	int num = 0;
	size_t len = sub.length();
	if (len == 0)len=1;//To deal with empty substring calls
	for (size_t i=0; (i=str.find(sub,i)) != std::string::npos; num++, i+=len);
	return num;
}

// https://www.techiedelight.com/split-string-cpp-using-delimiter/
void tokenize(std::string const &str, std::string& delim,
            std::vector<std::string> &out) {
  size_t start;
  size_t end = 0;

  int sz = count_nonoverlapping_sub(str, delim);
  out.clear();
  out.reserve(sz+1);

  while ((start = str.find_first_not_of(delim.c_str(), end)) 
                                                != std::string::npos) {
    end = str.find(delim.c_str(), start);
    out.push_back(str.substr(start, end - start));
  }
}

bool file_list::comp_nbrs(const double a, const double b) { 
  return a < b; 
}

void file_list::sort_nbrs(int left, int right) {
  if(right == -1 and left == -1) {
    left = 0;
    right = nbrs.size()-1;
  }
  //std::cout << "Sorting edge set " << e_set << std::endl;
  //std::cout << left << " " << right << std::endl;

  if (left < right) {
    int part = partition_nbrs(left, right);
    sort_nbrs(left, part - 1);
    sort_nbrs(part + 1, right);
  }
}
//https://en.wikipedia.org/wiki/Quicksort
//Function to determine the partitions
// partitions the array and returns the middle subscript
int file_list::partition_nbrs(int left, int right) {
  double pivot = nbrs.at(right);
  // move the mid point value to the front.
  int i = left-1;
  int j = left;
  for (; j < right; j++) {
    if(comp_nbrs(nbrs.at(j), pivot)) {
      i++;
      std::swap(nbrs.at(i), nbrs.at(j));
      std::swap(nbr_strs.at(i), nbr_strs.at(j));
    }
  }
  std::swap(nbrs.at(i+1),nbrs.at(right));
  std::swap(nbr_strs.at(i+1),nbr_strs.at(right));
  return i + 1;
}


file_list::file_list() : nbrs(0), nbr_strs(0, std::string("")) {
}

file_list::~file_list() {
  clear();
}

void file_list::clear() {
  for(int i = 0; i < nbr_strs.size(); i++)
    nbr_strs.at(i).clear();
  nbr_strs.clear();
  nbrs.clear();
}

int file_list::count_dir(const char *path) {
  int nbr = 0;
  struct dirent *entry;
  DIR *dir = opendir(path);
  if (dir == NULL) {
    return 0;
  }

  while ((entry = readdir(dir)) != NULL)
    nbr = nbr + 1;
  return nbr;
}

void file_list::list_dir(const char *path, const char *filename, const char *ext) {
  clear();
  int nbr = count_dir(path);
  nbrs.reserve(nbr);
  nbr_strs.reserve(nbr);

  struct dirent *entry;
  DIR *dir = opendir(path);
  if (dir == NULL) {
    return;
  }

  std::string s (filename);
  s = s + "[0-9]+[_\\.,][0-9]+";
  s = s + ext;
  const char* cstr = s.c_str();
  std::string nbr_temp("");

  while ((entry = readdir(dir)) != NULL) {
    //printf("%s\n", entry->d_name);
    std::string f (entry->d_name);

    // Match the number at the end of entry. 

    std::smatch sm;

    if(match_regex(f, sm, cstr)) {
      f = sm.str(0);

      if(match_regex(f, sm, "[0-9]+[_\\.,][0-9]+")) {
        nbr_temp = sm.str(0);
        nbr_strs.push_back(nbr_temp);

        nbr_temp = replace_sep(nbr_temp);
        const char* nstr = s.c_str();
        nstr = nbr_temp.c_str();
        nbrs.push_back(atof(nstr));
      }
    }
  }
  std::cout << nbrs.size() << " files" << std::endl;
  sort_nbrs();

  //for(int i = 0; i < nbrs.size(); i++)
  //  std::cout << nbr_strs.at(i) << " " << nbrs.at(i) << " "  << std::endl;

  closedir(dir);
}


// This is to be adapted from gmi_write_dmg, if needed be. Currently it is a 
// placeholder.
void gmi_write_tess(struct gmi_model* m, const char* filename)
{
  struct gmi_iter* it;
  struct gmi_ent* e;
  struct gmi_set* s;
  FILE* f = fopen(filename, "w");
  int i;
  /* entity counts */
  fprintf(f, "%d %d %d %d\n", m->n[3], m->n[2], m->n[1], m->n[0]);
  /* bounding box */
  fprintf(f, "0 0 0\n");
  fprintf(f, "0 0 0\n");
  /* vertices */
  it = gmi_begin(m, 0);
  while ((e = gmi_next(m, it))) {
    fprintf(f, "%d 0 0 0\n", gmi_tag(m, e));
  }
  gmi_end(m, it);
  /* edges */
  it = gmi_begin(m, 1);
  while ((e = gmi_next(m, it))) {
    s = gmi_adjacent(m, e, 0);
    fprintf(f, "%d ", gmi_tag(m, e));
    if (s->n == 0)
      fprintf(f,"-42 -42\n");
    else if (s->n == 1)
      fprintf(f,"%d -42\n", gmi_tag(m, s->e[0]));
    else
      fprintf(f, "%d %d\n",
          gmi_tag(m, s->e[0]), gmi_tag(m, s->e[1]));
    gmi_free_set(s);
  }
  gmi_end(m, it);
  /* faces */
  it = gmi_begin(m, 2);
  while ((e = gmi_next(m, it))) {
    /* we're going to cheat a bit here and
       treat all edges as one giant loop */
    fprintf(f, "%d 1\n", gmi_tag(m, e));
    s = gmi_adjacent(m, e, 1);
    fprintf(f, "%d\n", s->n);
    for (i = 0; i < s->n; ++i)
      fprintf(f, " %d 0\n", gmi_tag(m, s->e[i]));
    gmi_free_set(s);
  }
  gmi_end(m, it);
  /* regions */
  it = gmi_begin(m, 3);
  while ((e = gmi_next(m, it))) {
    /* same sort of cheat, all faces are one shell */
    fprintf(f, "%d 1\n", gmi_tag(m, e));
    s = gmi_adjacent(m, e, 2);
    fprintf(f, "%d\n", s->n);
    for (i = 0; i < s->n; ++i)
      fprintf(f, " %d 0\n", gmi_tag(m, s->e[i]));
    gmi_free_set(s);
  }
  gmi_end(m, it);
  fclose(f);
}

// Taken from
//http://www.cplusplus.com/forum/beginner/182904/

void ReadNumbers(const std::string& filename, char sep, 
                             std::vector<std::vector<int> > & output)
{
  std::ifstream src(filename);

  if (!src)
  {
    std::cerr << "\aError opening file.\n\n";
    exit(EXIT_FAILURE);
  }
  std::string buffer("");
  while(std::getline(src, buffer))
  {
    std::vector<int> v_curr(0);
    size_t strpos = 0;
    size_t endpos;
    endpos = buffer.find(sep);
    while (endpos != std::string::npos) {
      v_curr.push_back(std::stoi(buffer.substr(strpos,endpos - strpos)));
      strpos = endpos + 1;
      endpos = buffer.find(sep, strpos);
    }
    //std::cout << buffer.substr(strpos,buffer.length() - strpos) << std::endl;
    v_curr.push_back(std::stoi(buffer.substr(strpos,
                                          buffer.length() - strpos)));
    output.push_back(v_curr);
  }
}

void ReadNumbers(const std::string& filename, char sep, 
                          std::vector<std::vector<double > > & output)
{
  std::ifstream src(filename);

  if (!src)
  {
    std::cerr << "\aError opening file.\n\n";
    exit(EXIT_FAILURE);
  }
  std::string buffer("");
  while(std::getline(src, buffer))
  {
    std::vector<double > v_curr(0);
    size_t strpos = 0;
    size_t endpos;
    endpos = buffer.find(sep);
    while (endpos != std::string::npos) {
      v_curr.push_back(std::stod(buffer.substr(strpos,endpos - strpos)));
      strpos = endpos + 1;
      endpos = buffer.find(sep, strpos);
    }
    //std::cout << buffer.substr(strpos,buffer.length() - strpos) << std::endl;
    v_curr.push_back(std::stod(buffer.substr(strpos,
                                          buffer.length() - strpos)));
    output.push_back(v_curr);
  }
}

// Read lines from filename, and add the lines separated by sep into the output.
void ReadNames(const std::string& filename, const char * sep, 
                             std::vector<std::vector<std::string> > & output) {
  std::ifstream src(filename);

  if (!src)
  {
    std::cout << filename << std::endl;
    std::cerr << "\aError opening file.\n\n";
    exit(EXIT_FAILURE);
  }
  for(int i = 0; i < output.size(); i++)
    output.at(i).clear();
  output.clear();

  std::string buffer("");

  int count = 0;
  while(std::getline(src, buffer)) {
    count = count + 1;
  }
  src.clear(); //< Now we can read again
  std::ifstream src2(filename);

  output.reserve(count);

  std::vector<std::string> temp(0);
  while(std::getline(src2, buffer)) {
    temp.clear();
    size_t strpos = 0;
    size_t endpos;
    endpos = buffer.find(sep);
    while (endpos != std::string::npos) {
      temp.push_back(buffer.substr(strpos,endpos - strpos));
      strpos = endpos + 1;
      endpos = buffer.find(sep, strpos);
    }
    //temp.push_back(buffer.substr(strpos, buffer.length() - strpos));
    output.push_back(temp);
  }
}
