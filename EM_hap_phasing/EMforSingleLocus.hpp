#include <iostream>
#include <armadillo>
#include <string>
#include <list>
#include <vector>
#include <tr1/unordered_map>

using namespace std;
using namespace std::tr1;
using namespace arma;

struct GenotypeNode_t {
  int chrom,start,in_degree;
  string type;
  struct GenotypeNode_t* prev;
};

typedef struct GenotypeNode_t GT_tree_node;

typedef struct {
  double frac;
  GT_tree_node *pt;
} HT_vector_node;

typedef struct {
  int parentID;
  GT_tree_node *pt;
} temp_HT_vector_node;


typedef vector<temp_HT_vector_node> tempHTvec_t;
typedef vector<HT_vector_node> HTvec_t;
typedef vector<double> frac_t;
typedef unordered_map<string,frac_t> Umap_t;

