#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "Graph.h"

//This member function adds a new edge connecting x and y to the graph if it doesn't exist already
//It checks if the edge to be added already exists. If not then the edge is added at the end of the linked list and returns true
//If edge already exists the function return false
bool Graph::add_edge_directed (int x, int y, double input_dist, string input_name)
{
  //Creat new edge
  Edge* new_link = new Edge;
  new_link->neighbor=y;
  new_link->dist=input_dist;
  new_link->name=input_name;
  new_link->next=NULL;
  
  bool success=true;
  
  if (Edgelist[x] == NULL)
    {
      Edgelist[x]=new_link;
      E++;
    }
  else 
    {
      Edge* current_edge = Edgelist[x];
      while ((current_edge->next != NULL) && success)
	{
	  if (current_edge->neighbor != y)
	    {
	      current_edge = current_edge->next;
	    }
	  else 
	    {
	      success=false; //Edge already exists
	    }
	}
      if (success)
	{
	  current_edge->next = new_link;
	  E++;
	}
    }
  
  return (success);
}
//End of add_edge_directed

//This member function deletes an edge connecting x and y if it exists. If the edge does exist it return false
bool Graph::del_edge_directed (int x, int y)
{
  bool success=false;
  
  if (Edgelist[x] == NULL)
    {
      return false; //Edge doesn't exist
    }
  else
    {
      Edge* previous_edge = Edgelist[x];
      Edge* current_edge = previous_edge->next;
      while ((current_edge != NULL) && (! success))
	{
	  if (current_edge->neighbor == y)
	    {
	      previous_edge->next = current_edge->next; //Found the edge and delete it
	      delete current_edge;
	      E--;
	      success=true;
	    }
	  else
	    {
	      previous_edge = current_edge;
		current_edge  = current_edge->next;
	    }
	}
    }
  
  return success;
}
//End of del_edge_directed

//This member deletes all edges connected with node x
void Graph::del_all_edges_for_a_node (int x)
{
  if (Edgelist[x] == NULL) return;
  
  Edge* current_edge = Edgelist[x];
  Edge* next_edge    = current_edge->next;
  while (next_edge != NULL)
    {
      delete current_edge;
      current_edge=next_edge;
      next_edge=next_edge->next;
    }
  delete current_edge;
  Edgelist[x]=NULL;
  
  return;
}
//End of del_all_edges_for_a_node

//This member function returns the distance of edge (x,y) if it exists. If the edge doesn't exist it returns -1.0
double Graph::get_edge_value (int x, int y)
{
  double dist=-1.0;
  bool connect=false;
  
  if (Edgelist[x] == NULL)
    {
      return dist; //Edge doesn't exist
    }
  else
    {
      Edge* current_edge = Edgelist[x];
      while ((current_edge != NULL) && (! connect))
	{
	  if (current_edge->neighbor == y)
	    {
	      connect=true;
	      dist=current_edge->dist;
	    }
	  else
	    {
	      current_edge  = current_edge->next;
	    }
	}
    }
  
  return dist;
}
//End of get_edge_value

//This member function returns the name of edge (x,y) if it exists. If the edge doesn't exist it returns "".
string Graph::get_edge_name (int x, int y)
{
  string edge_name="";
  bool connect=false;

  if (Edgelist[x] == NULL)
    {
      return edge_name; //Edge doesn't exist
    }
  else
    {
      Edge* current_edge = Edgelist[x];
      while ((current_edge != NULL) && (! connect))
        {
          if (current_edge->neighbor == y)
            {
              connect=true;
              edge_name=current_edge->name;
            }
          else
            {
              current_edge  = current_edge->next;
            }
        }
    }

  return edge_name;
}
//End of get_edge_name

//This function sets the distance for edge (x,y) and returns true if the edge exists
//It returns false if the edge doesn't exist
bool Graph::set_edge_dist (int x, int y, double input_dist)
{
  bool success=false;
  
  if (Edgelist[x] == NULL)
    {
      return false; //Edge doesn't exist
    }
  else
    {
      Edge* current_edge = Edgelist[x];
      while ((current_edge != NULL) && (! success))
	{
	  if (current_edge->neighbor == y)
	    {
	      current_edge->dist = input_dist; //Found the edge and delete it
	      success=true;
	    }
	  else
	    {
	      current_edge  = current_edge->next;
	    }
	}
    }
  
  return success;
}
//End of set_edge_dist

//This function returns all the neighbors of node x (and the distance of the edges)
void Graph::adjacency_list (int x, vector <adjacency> &ad_list)
{
  if (Edgelist[x] == NULL) {return;}
  else 
    {
      Edge* current_edge = Edgelist[x];
      while (current_edge != NULL) 
	{
	  adjacency current;
	  current.neighbor = current_edge->neighbor;
	  current.dist     = current_edge->dist;
	  ad_list.push_back(current);
	  current_edge = current_edge->next;
	}
    }
  return;
}
//End of adjacency_list

//Get Graph node attributes
bool Graph::get_node_attributes (int x, Node* node)
{
  if ((x >= 0) && (x < N))
    {
      node->chrom  = Nodelist[x].chrom;
      node->coor   = Nodelist[x].coor;
      node->type   = Nodelist[x].type;
      node->strand = Nodelist[x].strand;
      return true;
    }
  else {return false;}
}
//End of get_node_attributes

//Constructor. Read from files
Graph::Graph (string samplename)
{
  N=0; E=0;
  string nodefile=samplename+'.'+"nodes";
  string edgefile=samplename+'.'+"edges";

  //Get nodes
  ifstream ifs;
  ifs.open(nodefile.c_str(), ifstream::in);
  string line;
  while (getline(ifs,line))
    {
      vector <string> entry;
      boost::split(entry, line, boost::is_any_of(":"));
      Node tmp_node;
      tmp_node.chrom = entry[0];
      tmp_node.coor  = boost::lexical_cast <unsigned long> (entry[1]);
      if (entry[2] == "+") {tmp_node.strand=PLUS;}
      else {tmp_node.strand=MINUS;}

      tmp_node.type=IN;
      Nodelist.push_back(tmp_node);
      Edgelist.push_back(NULL);
      N++;
      //cout << N << "\t" << tmp_node.chrom << "\t" << tmp_node.coor << "\t" << tmp_node.type << "\t" << tmp_node.strand << endl;

      tmp_node.type = OUT;
      Nodelist.push_back(tmp_node);
      Edgelist.push_back(NULL);
      N++;
      //cout << N << "\t" << tmp_node.chrom << "\t" << tmp_node.coor << "\t" << tmp_node.type << "\t" << tmp_node.strand << endl;
    }
  ifs.close();
  
  //Add breakpoint edges
  ifs.open(edgefile.c_str(), ifstream::in);
  while (getline(ifs,line))
    {
      vector <string> entry;
      boost::split(entry, line, boost::is_any_of("\t"));
      int source = boost::lexical_cast <int> (entry[1])*2;
      int sink   = boost::lexical_cast <int> (entry[2])*2;
      add_edge_directed(source, sink+1, 1.0, entry[0]);
      add_edge_directed(sink, source+1, 1.0, entry[0]);
    }
  ifs.close();

  //Add adjancency edges
  for (int i=0; i < N; i++) 
    {
      //cout << i << "\t" << Nodelist[i].chrom << "\t" << Nodelist[i].coor << "\t" << Nodelist[i].type << "\t" << Nodelist[i].strand << endl;
      if (Nodelist[i].type == OUT)
        {
          for (int j=0; j < N; j++)
            {
              if ((Nodelist[j].type == IN) && (Nodelist[i].chrom == Nodelist[j].chrom))
                {
                  if ((Nodelist[i].strand == PLUS) && (Nodelist[j].strand == MINUS) && (Nodelist[i].coor <= Nodelist[j].coor)
                      && (Nodelist[j].coor-Nodelist[i].coor <= 100000))
                    {
                      add_edge_directed(i, j, static_cast<double> (Nodelist[j].coor-Nodelist[i].coor), "AD");
                    }
                  if ((Nodelist[i].strand == MINUS) && (Nodelist[j].strand == PLUS) && (Nodelist[i].coor >= Nodelist[j].coor)
                      && (Nodelist[i].coor-Nodelist[j].coor <= 100000))
                    {
                      add_edge_directed(i, j, static_cast<double> (Nodelist[j].coor-Nodelist[i].coor), "AD");
                    }
                }
            }
        }
    }

}
//End of Graph

void Graph::Initial_graph (int NodesNum)
{
  N = NodesNum;
  E = 0;
  for (int i=0; i < NodesNum; ++i)
    {
      Node tmp_node;
      tmp_node.chrom="";
      tmp_node.coor=0;
      tmp_node.type=IN;
      tmp_node.strand=PLUS;
      Nodelist.push_back(tmp_node);
      Edgelist.push_back(NULL);
    }
}
//End of Graph(int)

Graph::~Graph ()
{
  for (int i=0; i < N; i++)
    {
      del_all_edges_for_a_node(i);
    }
  //cout << "Graph destructor invoked" << endl;
}
//End of destructor

