#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include "Graph.h"
#include "ShortestPath.h"
#include "PriorityQueue.h"

int main(int argc, char** argv)
{
  //cout << argv[1] << endl;
  Graph graph(argv[1]);
  //cout << graph.get_node_number() << "\t" << graph.get_edge_number() << endl;

  for (int i=0; i < graph.get_node_number(); i+=2)
    {
      ShortestPath SP(i);
      SP.Dijkstra(&graph);
      for (int j=0; j < graph.get_node_number(); j++)
	{
	  if (i != j)
	    {
	      SP.print_path(j, &graph);
	      //cout << i << "\t" << j << endl;
	    }
	}
    }
}
