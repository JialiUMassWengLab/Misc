#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "ShortestPath.h"

ShortestPath::ShortestPath (int x) 
{
  start_node=x;
  PathNode new_path;
  new_path.end_node=x;
  new_path.num_bkpt=0;
  new_path.path_size=0.0;
  new_path.route=static_cast<ostringstream*>( &(ostringstream() << start_node) )->str();;
  new_path.bkpt="";
  paths.push_back(new_path);
}

//Check if node y is already in closed set
bool ShortestPath::if_in_closed_set (int y)
{
  int i=0;
  while ((paths[i].end_node != y) && (i < paths.size())) i++;
  if (i < paths.size()) return true;
  else return false;
}

//Return the distance of the shortest path between start node and node y
double ShortestPath::get_path_size (int y) 
{
  int i=0;
  while ((paths[i].end_node != y) && (i < paths.size())) i++;
  if (i < paths.size()) return paths[i].path_size;
  else return -1.0;
}

//Print the shortest path from start node to node y
void ShortestPath::print_path (int y, Graph* const G) 
{
  int i=0;
  while ((paths[i].end_node != y) && (i < paths.size())) i++;
  if ((i < paths.size()) && (paths[i].num_bkpt >= 2))
    {
      //cout << "An alternating path from  " << start_node << " to " << y << " is " << paths[i].route << " with distance " << paths[i].path_size << endl;
      vector <string> entry;
      boost::split(entry, paths[i].route, boost::is_any_of("-"));
      for (vector <string>::iterator it = entry.begin(); it < entry.end(); it++)
	{
	  Node curr_node;
	  G->get_node_attributes(boost::lexical_cast <int> (*it), &curr_node);
	  string direction="+";
	  if (curr_node.strand == MINUS) {direction="-";}
	  cout << curr_node.chrom << ":" << curr_node.coor << direction;

	  if ((curr_node.type == IN) && (it < entry.end()-1))
	    {
	      cout << "::";
	    }
	  else if ((curr_node.type == OUT) && (it < entry.end()-1))
	    {
	      if (curr_node.strand == PLUS) {cout << "-->";}
	      else {cout << "<--";}
	    }
	  else if (it == entry.end()-1)
	    {
	      cout << endl;
	    }
	}

      cout << paths[i].bkpt << endl;
    }
  /*
  else
    {
      cout << "Node " << start_node << " and " << y << " are NOT connected" << endl;
    }
  */
  return;
}

//Inplementation of Dijkstra algorithm
void ShortestPath::Dijkstra (Graph* const G)
{
  int NodeNum = G->get_node_number();
  PriorityQueue PQ(NodeNum);
  
  PQ.insert_or_update_all_neighbors(start_node, G);
  
  int node=0; int parent=0; double dist=0.0;
  while (PQ.pop_PQ_top(node,parent,dist))
    {
      if (! if_in_closed_set(node))
	{
	  //cout << node << endl;
	  string new_route;
	  int j=0;
	  while ((j < paths.size()) && (paths[j].end_node != parent)) j++;
	  new_route = paths[j].route;
	  
	  string new_node  = static_cast<ostringstream*> (&(ostringstream() << node) )->str();
	  new_route += string("-");
	  new_route += new_node;
	  
	  //cout << new_route << endl;
	  PathNode new_path_node;
	  new_path_node.end_node=node;
	  new_path_node.num_bkpt=paths[j].num_bkpt;
	  if (node % 2 == 1) {new_path_node.num_bkpt++;}
	  new_path_node.path_size=dist+paths[j].path_size;
	  new_path_node.route=new_route;
	  new_path_node.bkpt = paths[j].bkpt + G->get_edge_name(parent, node);
	  paths.push_back(new_path_node);
	  
	  PQ.insert_or_update_all_neighbors(node, G);
	}
    }

}
//End of Dijkstra
