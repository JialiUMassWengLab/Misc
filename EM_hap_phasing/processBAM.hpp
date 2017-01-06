#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include "sam.h"
#include "bam.h"

using namespace std;

typedef struct {
  string barcode;
  string seq;
  vector<int> qual;
} readGT;

typedef struct {
  int chrom,start,len;
  vector<readGT*> types;
} locusGT;

static int fetch_func(const bam1_t* b, void *data);
locusGT* get_genotype_info(samfile_t* inbam, bam_index_t* idx, const char* locus_coor);
vector<string> split(string str, char delimiter);
