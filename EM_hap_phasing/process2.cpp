#include <stdio.h>
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

static int fetch_func(const bam1_t* b, void *data)
{
  locusGT *GT = (locusGT*) data;
  uint8_t *tag = bam_aux_get(b,"BX");
  if (tag != NULL)
    {
      uint32_t readEnd = bam_calend(&(b->core),bam1_cigar(b));
      if (GT->start+GT->len > readEnd || GT->start < b->core.pos) {
	return 1;
      }//Read too short to include the whole variant position
      
      //This part finds the position in the read correspond to the variant locus
      int dist2var = 0, i = 0;
      int curr_pos = b->core.pos;
      while (curr_pos < GT->start) 
	{
	  uint32_t curr_cigar = bam1_cigar(b)[i];
	  char cigarchr = bam_cigar_opchr(curr_cigar);
	  int cigarlen = bam_cigar_oplen(curr_cigar);
	  if (cigarchr == 'S' || cigarchr == 'I')
	    {
	      dist2var += cigarlen;
	    }
	  else if (cigarchr == 'D') 
	    {
	      curr_pos += cigarlen;
	      if (curr_pos > GT->start) {return 1;} //locus deleted in this read
	    }
	  else if (cigarchr == 'M')
	    {
	      if (curr_pos + cigarlen > GT->start) {
		dist2var += GT->start - curr_pos;
	      }
	      else {
		dist2var += cigarlen;
	      }
	      curr_pos += cigarlen;
	    }
	  i++;
	}

      //Decides the length of sequence to extract from the read (take care of indels)
      int extractLen = GT->len;
      if (curr_pos == GT->start+1 && curr_pos < readEnd) 
	{
	  uint32_t curr_cigar = bam1_cigar(b)[i];
	  char cigarchr = bam_cigar_opchr(curr_cigar);
	  int cigarlen = bam_cigar_oplen(curr_cigar);
	  if (cigarchr == 'D') extractLen -= cigarlen;
	  else if (cigarchr == 'I') extractLen += cigarlen;
	}
      if (extractLen <= 0) return 2;

      readGT *cur_read = new(readGT);
      cur_read->barcode = bam_aux2Z(tag);
      cur_read->barcode = cur_read->barcode.substr(0,14);
      cur_read->seq = "";
      char *s  = (char*) bam1_seq(b);
      char *qual  = (char*) bam1_qual(b);
      for(int n=0; n<extractLen; n++) {
	char v = bam1_seqi(s,n+dist2var);
	cur_read->seq += bam_nt16_rev_table[v];
	cur_read->qual.push_back((int) qual[n+dist2var]);
      }
      //cout << cur_read->seq << endl;
      //cout << cur_read->barcode << "\t" << cur_read->seq << endl;
      //cout << GT->types.size() << endl;
      //cout << cur_read << endl;
      GT->types.push_back(cur_read);
    }

  return 0;
}

int main(int argc, char* argv[])
{
  printf("%s\n", argv[1]);
  samfile_t* inbam = NULL;

  inbam = samopen(argv[1], "rb", 0);

  int ref,start,end;
  bam_index_t *idx;
  idx = bam_index_load(argv[1]);
  if (idx == 0) {
    fprintf(stderr, "BAM indexing file is not available.\n");
    return 1;
  }
  bam_parse_region(inbam->header, argv[2], &ref,
		   &start, &end);
  if (ref < 0 || start < 0) {
    fprintf(stderr, "Invalid region %s\n", argv[2]);
    return 1;
  }

  printf("%d\t%d\t%d\n", ref, start, end);
  locusGT GT;
  GT.chrom = ref;
  GT.start = start;
  GT.len = end-start;
  GT.types = vector<readGT*>();

  bam_fetch(inbam->x.bam, idx, ref, start, end, &GT, fetch_func);

  for (vector<readGT*>::iterator it = GT.types.begin(); it != GT.types.end(); it++) {
    cout << (*it)->barcode << "\t" << (*it)->seq << endl;
    for (vector<int>::iterator it2 = (*it)->qual.begin(); it2 != (*it)->qual.end(); it2++) {
      printf(" %d",*it2);
    }
    printf("\n");
    delete(*it);
  }

  samclose(inbam);
  return 0;
}
