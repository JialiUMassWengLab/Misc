#include <stdio.h>
#include <cstring>
#include "sam.h"
#include "bam.h"

using namespace std;

typedef struct {
  char barcode[20];
  char seq[10];
  int qual[10];
} readGT;

typedef struct {
  int chrom,start,len,nRead;
  readGT types[500];
} locusGT;

static int fetch_func(const bam1_t* b, void *data)
{
  locusGT *GT = (locusGT*) data;
  uint8_t *tag = bam_aux_get(b,"BX");
  if (tag != NULL)
    {
      uint32_t cigar_1st = bam1_cigar(b)[0];
      char cigarchr = bam_cigar_opchr(cigar_1st);
      int cigarlen = bam_cigar_oplen(cigar_1st);
      int dist2var = GT->start - b->core.pos;
      if (cigarchr == 'S') 
	{
	  dist2var += cigarlen;
	}

      uint32_t readEnd = bam_calend(&(b->core),bam1_cigar(b));
      if (b->core.pos+dist2var+GT->len > readEnd+1) {
	return 1;
      }//Read too short to include the whole variant position
      
      int n=0;
      readGT *cur_read = &(GT->types[GT->nRead]);
      strncpy(cur_read->barcode,bam_aux2Z(tag),14);
      char *s  = (char*) bam1_seq(b);
      char *qual  = (char*) bam1_qual(b);
      for(n=0;n<(GT->len);n++) {
	char v = bam1_seqi(s,n+dist2var);
	strncat(cur_read->seq,&bam_nt16_rev_table[v],1);
	cur_read->qual[n] = (int) qual[n+dist2var];
      }

      GT->nRead++;
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
  locusGT *GT = (locusGT*) malloc(sizeof(locusGT));
  GT->chrom = ref;
  GT->start = start;
  GT->len = end-start;
  GT->nRead = 0;

  bam_fetch(inbam->x.bam, idx, ref, start, end, GT, fetch_func);

  int i,j;
  for (i=0; i<GT->nRead; i++) {
    readGT curr = GT->types[i];
    printf("%s\t%s\n",curr.barcode,curr.seq);
    for (j=0; j<GT->len; j++) {
      printf(" %d",curr.qual[j]);
    }
    printf("\n");
  }

  free(GT);
  samclose(inbam);
  return 0;
}
