#ifndef _SEH_GRAPH_
#define _SEH_GRAPH_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

using namespace std;

enum node_type {IN,OUT};
enum node_direc {PLUS,MINUS};

typedef struct adjacency
{
  int neighbor;
  double dist;
} adjacency;

typedef struct Edge
{
  int neighbor;
  string name;
  double dist;
  Edge* next;
} Edge;

typedef struct Node
{
  unsigned long coor;
  string chrom;
  node_type type;
  node_direc strand;
} Node;

class Graph
{
  int N,E; //Number of nodes N and number of edges E
  vector <Node> Nodelist;
  vector <Edge*> Edgelist; //linklist representation for the adjacency

 public:
  inline int get_node_number () {return N;}
  inline int get_edge_number () {return E;}

  bool add_edge_directed (int x, int y, double input_dist, string input_name);
  bool del_edge_directed (int x, int y);
  void del_all_edges_for_a_node (int x);
  double get_edge_value (int x, int y);
  string get_edge_name (int x, int y);
  bool set_edge_dist (int x, int y, double input_dist);
  void adjacency_list (int x, vector <adjacency> &ad_list);
  bool get_node_attributes (int x, Node* node);
  void Initial_graph (int NodeNum);
  Graph (string samplename);
  Graph () {E=0;}
  ~Graph ();
};


#endif // _SEH_GRAPH_
