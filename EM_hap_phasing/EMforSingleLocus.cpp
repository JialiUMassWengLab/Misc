#include "EMforSingleLocus.hpp"
#include "processBAM.hpp"
#include <fstream>
#include <deque>
#include <set>

#define PRUNE_THRE 0.05
#define PRINT_THRE 0.1
#define ITER_TERM_THRE 0.000001

int pruneTreeBranch(GT_tree_node *leaf_node)
{
  GT_tree_node *curr_node = leaf_node;
  while (curr_node->prev != NULL && curr_node->in_degree == 0)
    {
      GT_tree_node *prev_node = curr_node->prev;
      delete(curr_node);
      prev_node->in_degree--;
      curr_node = prev_node;
    }

  return 0;
}

int printHaplotypes(HTvec_t &current, GT_tree_node* root_node, bam_header_t *header)
{
  int i=0;
  for (HTvec_t::iterator it=current.begin(); it != current.end(); it++)
    {
      if (it->frac > PRINT_THRE)
	{
	  cout << it->frac << endl;
	  GT_tree_node *curr_node = it->pt;
	  while (curr_node->prev != NULL) 
	    {
	      cout << header->target_name[curr_node->chrom] << ":" << curr_node->start << ":" << curr_node->type << "\t";
	      curr_node = curr_node->prev;
	    }	  
	  cout << endl;
	}

      i++;
    }

  return 0;
}

double EMiteration(vec &frac, Umap_t &temp_mols, int unique)
{
  vec prev_f = frac;
  vec pseudo_reads = zeros<vec>(frac.n_elem);
  double total_reads = 0;
  for (Umap_t::iterator it = temp_mols.begin(); it != temp_mols.end(); it++)
    {
      uvec indices = find(vec(it->second)-max(vec(it->second))==0);
      double weight = indices.n_elem > 1 ? unique/(double)temp_mols.size() : 1.0;
  
      vec norm_read = normalise(prev_f % vec(it->second), 1);
      if (max(norm_read) > 0.99) 
	{
	  uword index = index_max(norm_read);
	  norm_read.zeros();
	  norm_read[index]=1;
	}
      else
	{
	  for (unsigned i=0; i<norm_read.n_elem; i++)
	    if (norm_read[i] < 0.01) norm_read[i] = 0;
	  norm_read = normalise(norm_read,1);
	}

      pseudo_reads += norm_read * weight;
      total_reads += weight;
    }
  
  //cout << pseudo_reads << endl;
  //vec new_f = pseudo_reads/(float)temp_mols.size();
  vec new_f = pseudo_reads/total_reads;
  frac = new_f;

  return norm(new_f - prev_f, 2);
}

vec initialize(HTvec_t &current, tempHTvec_t &tempHTvec, Umap_t &mol_prob, Umap_t &temp_mols, locusGT *GT, 
	       vector<string> &genotypes, bool &ifProceed, int &unique)
{
  tempHTvec.clear();
  int i=0;
  for (HTvec_t::iterator it = current.begin(); it != current.end(); it++)
    {
      for (vector<string>::iterator it2 = genotypes.begin(); it2 != genotypes.end(); it2++)
	{
	  GT_tree_node* new_tree_node = new(GT_tree_node);
	  new_tree_node->type = (*it2);
	  new_tree_node->chrom = GT->chrom;
	  new_tree_node->start = GT->start;
	  new_tree_node->prev = it->pt;
	  it->pt->in_degree++;

	  temp_HT_vector_node new_tempHTvec_node;
	  new_tempHTvec_node.parentID = i;
	  new_tempHTvec_node.pt = new_tree_node;
	  tempHTvec.push_back(new_tempHTvec_node);
	}
      i++;
    }

  for (Umap_t::iterator it = mol_prob.begin(); it != mol_prob.end(); it++)
    {
      string bc = it->first;
      temp_mols[bc].clear();
      for (tempHTvec_t::iterator it2 = tempHTvec.begin(); it2 != tempHTvec.end(); it2++)
	{
	  double new_frac = it->second[it2->parentID];
	  temp_mols[bc].push_back(new_frac);
	}
    }
  
  Umap_t overlap_mols;
  for (vector<readGT*>::iterator it = GT->types.begin(); it != GT->types.end(); it++) 
    {
      string bc = (*it)->barcode;
      Umap_t::iterator itt = mol_prob.find(bc);
      if (itt == mol_prob.end()) 
	{
	  temp_mols[bc].clear();
	  for (tempHTvec_t::iterator it2 = tempHTvec.begin(); it2 != tempHTvec.end(); it2++)
	    temp_mols[bc].push_back( 1.00 );
	}
      
      int i = 0;
      for (tempHTvec_t::iterator it2 = tempHTvec.begin(); it2 != tempHTvec.end(); it2++)
	{
	  GT_tree_node* tree_pt = it2->pt;
	  double prob_loc=1.0;
	  int lenDiff=abs((int)(*it)->qual.size() - (int)tree_pt->type.size());
	  int minTypeLen = tree_pt->type.size() < (*it)->qual.size() ? tree_pt->type.size() : (*it)->qual.size();
	  
	  //Penalize for mismatch (using sequencing quality)
	  for (int j=0; j<minTypeLen; j++)
	    {
	      double error_rate = pow(10.0, -(*it)->qual[j]/10);
	      if (tree_pt->type[j] != (*it)->seq[j]) prob_loc *= error_rate/3;
	      else prob_loc *= 1 - error_rate;
	    }
	  
	  //Penalize for indels
	  for (int j=0; j<lenDiff; j++)
	    prob_loc *= 0.0005;
	  
	  temp_mols[bc][i] *= prob_loc;
	  i++;
	}

      if (itt != mol_prob.end()) overlap_mols[bc] = temp_mols[bc];
    }
       
  //Remove bacodes that are unlikely to support any of the current haplotypes
  for (Umap_t::iterator it = temp_mols.begin(); it != temp_mols.end(); it++)
    {
      if (max(vec(it->second)) < 0.01)
	{
	  mol_prob.erase(it->first);
	  overlap_mols.erase(it->first);
	  temp_mols.erase(it);
	}
    }

  //Initialize frac
  vec ini_frac(tempHTvec.size());
  for (unsigned int i = 0; i < ini_frac.n_elem; i++)
    ini_frac[i] = current[tempHTvec[i].parentID].frac / (float) genotypes.size();      
  //ini_frac.ones();
  ini_frac = normalise(ini_frac,1);
  
  if (!mol_prob.empty()) 
    {
      vec covered = zeros<vec>(current.size());
      set<int> path_supported_by_unique;
      for (Umap_t::iterator it = overlap_mols.begin(); it != overlap_mols.end(); it++)
	{
	  uvec indices = find(vec(it->second)-max(vec(it->second))==0);
	  for (unsigned int j=0; j < indices.n_elem; j++)
	    {
	      covered[tempHTvec[indices[j]].parentID] += 1.0;
	    }
	  if (indices.n_elem == 1) 
	    {
	      unique++;
	      path_supported_by_unique.insert(tempHTvec[indices[0]].parentID);
	    }
	  //cout << it->first << endl << vec(it->second) << endl;
	}
      //cout << covered << endl;
      ifProceed = all(covered > 0) && unique > 0 && path_supported_by_unique.size() > 1;

      if (ifProceed) EMiteration(ini_frac,overlap_mols,unique);
    }
  else {ifProceed = true;}

  return ini_frac;
}

int update_locusInfo(HTvec_t &current, tempHTvec_t &tempHTvec, Umap_t &temp_mols, Umap_t &mol_prob, vec &frac)
{
  current.clear();
  for (unsigned int i=0; i<frac.n_elem; i++)
    if (frac[i] <= PRUNE_THRE) frac[i]=0.0;
  frac = normalise(frac,1);
      
  int i = 0;
  for (tempHTvec_t::iterator it = tempHTvec.begin(); it != tempHTvec.end(); it++)
    {
      if (frac[i] == 0) pruneTreeBranch(it->pt);
      else 
	{
	  HT_vector_node new_HTvec_node;
	  new_HTvec_node.frac = frac[i];
	  new_HTvec_node.pt = it->pt;
	  current.push_back(new_HTvec_node);
	}
      i++;
    }

  for (Umap_t::iterator it = temp_mols.begin(); it != temp_mols.end(); it++)
    {
      string bc = it->first;
      mol_prob[bc].clear();
      for (unsigned int j=0; j < frac.n_elem; j++)
	if (frac[j] > 0) mol_prob[bc].push_back(it->second[j]);

      //mol_prob[bc] = conv_to< frac_t >::from(normalise(vec(mol_prob[bc]), 1));
    }

  return 0;
}

int main(int argc, char* argv[])
{
  GT_tree_node* root_node = new(GT_tree_node);
  root_node->type = "";
  root_node->chrom = -1;
  root_node->start = -1;
  root_node->in_degree = 0;
  root_node->prev = NULL;
  
  HTvec_t curr_HTvec;
  HT_vector_node first_HTvec_node;
  first_HTvec_node.frac = 1.0;
  first_HTvec_node.pt = root_node;
  curr_HTvec.push_back(first_HTvec_node);

  tempHTvec_t tempHTvec;
  Umap_t mol_prob,temp_mols;  

  //printf("%s\n", argv[1]);
  samfile_t* inbam = NULL;
  inbam = samopen(argv[1], "rb", 0);

  bam_index_t *idx;
  idx = bam_index_load(argv[1]);
  

  //This loop iterates through each locus
  string line;
  ifstream infile (argv[2]);
  if (infile.is_open())
    {
      deque<string> lines;
      while (getline(infile, line))
	{
	  lines.push_back(line);
	}

      for (unsigned int i=0; i<lines.size()/8; i++)
	{
	  line = lines.front();
	  lines.pop_front();
	  lines.push_back(line);
	}

      string last_unprocessed = "";
      bool any_processed = true;
      while (!lines.empty())
	{	  
	  line = lines.front();
	  lines.pop_front();
	  if (line == last_unprocessed) break;

	  vector<string> fields = split(line, ',');
	  cout << fields[0].c_str() << endl;
	  locusGT* GT = get_genotype_info(inbam, idx, fields[0].c_str());
	  if (GT == NULL) return 1;
  
	  vector<string> genotypes = split(fields[1], ';');

	  bool ifProceed = false;
	  int unique = 0;
	  vec frac = initialize(curr_HTvec,tempHTvec,mol_prob,temp_mols,GT,genotypes,ifProceed,unique);
	  //cout << ifProceed << endl << unique << "\t" << temp_mols.size() << "\t" << unique/(double)temp_mols.size() << endl << frac << endl;
	  if (! ifProceed)
	    {
	      lines.push_back(line);
	      if (any_processed) 
		{
		  last_unprocessed = line;
		  any_processed = false;
		}
	      continue;
	    }
	  any_processed = true;
	  /*
	  for (Umap_t::iterator it=temp_mols.begin(); it!=temp_mols.end(); it++)
	    cout << it->first << endl << vec(it->second) << endl;
	  */
	  //EM iterations done here
	  int k=0;
	  while (EMiteration(frac,temp_mols,unique) >= ITER_TERM_THRE)
	    {
	      //printf("%d\n", k);
	      //cout << frac << endl << endl;
	      k++;
	    }

	  /*
	  cout << "After update\t" << k << " iterations" << endl;
	  cout << frac << endl;
	  for (Umap_t::iterator it=temp_mols.begin(); it!=temp_mols.end(); it++)
	    cout << it->first << endl << normalise(frac % vec(it->second), 1) << endl;
	  */
	  update_locusInfo(curr_HTvec,tempHTvec,temp_mols,mol_prob,frac);
	  /*
	  cout << "After update\t" << k << " iterations" << endl;
	  for (Umap_t::iterator it=mol_prob.begin(); it!=mol_prob.end(); it++)
	    cout << it->first << endl << vec(it->second) << endl;
	  */
	  for (HTvec_t::iterator it=curr_HTvec.begin(); it!=curr_HTvec.end(); it++)
	    cout << it->frac << "\t";
	  cout << endl;
	}
      infile.close();

      printHaplotypes(curr_HTvec, root_node, inbam->header);

      //Delete the tree
      for (HTvec_t::iterator it=curr_HTvec.begin(); it!=curr_HTvec.end(); it++) pruneTreeBranch(it->pt);
      delete(root_node);
    }
  else cout << "Can't open file!" << endl;

}
