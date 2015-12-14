#ifndef _SEH_SHORTEST_PATH_
#define _SEH_SHORTEST_PATH_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "Graph.h"
#include "PriorityQueue.h"
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

using namespace std;

class ShortestPath
{
  //Class ShortestPath stores shortest paths from the start node to all the other nodes in the graph

  typedef struct PathNode
  {
    int end_node;
    int num_bkpt;
    double path_size;
    string route;
    string bkpt;
  } PathNode;

  int start_node;
  vector <PathNode> paths;

 public:
  ShortestPath (int x);
  bool if_in_closed_set (int y);
  double get_path_size (int y);
  void print_path (int y, Graph* const G);
  void Dijkstra (Graph* const G);
};

#endif // _SEH_SHORTEST_PATH


