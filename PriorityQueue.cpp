#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "Graph.h"
#include "PriorityQueue.h"

bool PriorityQueue::check_if_element_exists (int node)
{
  if (current_dist[node] >= 0) return true;
  else return false;
}

//Pop the top element of the Priority Queue and returns the edge and its distance
bool PriorityQueue::pop_PQ_top (int &node, int &parent, double &dist)
{
  if (queue == NULL) return false;
  
  node   = queue->nodeID;
  parent = queue->parentID;
  dist   = queue->dist;
  
  PQNode* current = queue;
  queue = queue->next;
  delete current;
  queue_size--;
  
  return true;
}
//End of pop_PQ_top 

//Inserts an edge into the Priority Queue
//Checks if the distance >=0 and if the distance is shorter than the existing distance
bool PriorityQueue::insert_or_update_node (int ip_node, int ip_parent, double ip_dist)
{
  if (ip_dist < 0) return false;
  if ((current_dist[ip_node] >= 0) && (current_dist[ip_node] <= ip_dist)) return false;
  
  PQNode* new_node = new PQNode;
  new_node->nodeID   = ip_node;
  new_node->parentID = ip_parent;
  new_node->dist     = ip_dist;
  new_node->next     = NULL;
  
  PQNode* previous_node = queue;
  if ((previous_node == NULL) || (previous_node->dist >= ip_dist))
    {
      new_node->next=previous_node;
      queue=new_node;
      current_dist[ip_node]=ip_dist;
      queue_size++;
      return true;
    }
  
  PQNode* next_node  = previous_node->next;
  while ((next_node != NULL) && (next_node->dist < ip_dist))
    {
      previous_node = next_node;
      next_node     = next_node->next;
    }
  new_node->next      = next_node;
  previous_node->next = new_node;
  current_dist[ip_node]=ip_dist;
  queue_size++;
  
  while ((next_node != NULL) && (next_node->nodeID != ip_node))
    {
      previous_node = next_node;
      next_node     = next_node->next;
    }
  
  if (next_node != NULL)
    {
      previous_node->next = next_node->next;
      delete next_node;
      queue_size--;
    }
  
  return true;
}
//End of insert_or_update_node

//When a node is put into closed set, this function updates all its neighbors to the open set
void PriorityQueue::insert_or_update_all_neighbors(int parent, Graph* const G)
{
  vector <adjacency> ad_list;
  G->adjacency_list(parent, ad_list);
  
  for (int i=0; i < ad_list.size(); ++i)
    {	
      insert_or_update_node(ad_list[i].neighbor,parent,ad_list[i].dist);
      //cout << ad_list[i].neighbor << "\t" << ad_list[i].dist << endl;
      //cout << get_PQ_size() << endl;
    }
}
//End of insert_or_update_all_neighbors

void PriorityQueue::release_queue()
{    
  if (queue == NULL) return;
  
  PQNode* cursor=queue->next;
  while (cursor != NULL)
    {
      delete queue;
      queue=cursor;
      cursor=cursor->next;
    }
  delete queue;
  queue=NULL;
  
  return;
}
//End of release_queue

PriorityQueue::PriorityQueue (int NodeNum)
{
  queue_size=0;
  queue=NULL;
  for (int i=0; i < NodeNum; ++i)
    {
      current_dist.push_back(-1.0);
    }
}
//End of constructor

PriorityQueue::~PriorityQueue ()
{
  release_queue();
  //cout << "Priority Queue destructor invoked" << endl;
}
//End of destructor
