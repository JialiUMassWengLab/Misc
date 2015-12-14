#ifndef _SEH_PRIORITYQUEUE_
#define _SEH_PRIORITYQUEUE_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>

using namespace std;

class PriorityQueue
{
  typedef struct PQNode
  {
    int nodeID,parentID;
    double dist;
    PQNode* next;
  } PQNode;

  int queue_size;

  //Priority queue
  PQNode* queue;

  //Stores the shortest distance from the start node to a certain node i (element current_dist[i])
  //If node i is in neither open set nor closed set, current_dist[i]=-1.0
  vector <double> current_dist;

 public:
  inline int get_PQ_size() {return queue_size;}

  bool check_if_element_exists (int node);
  bool pop_PQ_top (int &node, int &parent, double &dist);
  bool insert_or_update_node (int ip_node, int ip_parent, double ip_dist);
  void insert_or_update_all_neighbors(int parent, Graph* const G);
  void release_queue();
  PriorityQueue (int NodeNum);
  ~PriorityQueue ();
};


#endif // _SEH_PRIORITYQUEUE_
