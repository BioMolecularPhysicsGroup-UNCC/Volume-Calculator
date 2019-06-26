/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_algorithms_Fast_Graph_cpp
#define bmpg_uncc_edu_algorithms_Fast_Graph_cpp

#include <iostream>
#include <iterator>
#include <bmpg_uncc_edu/algorithms/FastGraph.hpp>

namespace bmpg_uncc_edu {
namespace algorithms {
using namespace std;




FastGraph::FastGraph()
{
	
}

/**
	Add a vertex to the graph. The vertex is duplicated if it does not exist in the graph.
	@param v vertex
	@return pointer to the duplicated vertex in the graph
*/
FastGraph::vertex_descriptor_t FastGraph::add_vertex(const vertex_descriptor_t& v)
{
	if(look_up(v))
		return v;							//found

	for(size_t i = nbr_list_vec_.size(); i <= v; i++)
		nbr_list_vec_.push_back(nbr_list_descriptor_t(new nbr_list_t));
	
	return v;
}

/**
	Add an edge to the graph. The edge and the vertices corresponding to the edge are duplicated.
	@param e edge
*/

void FastGraph::add_edge(const edge_t& e)
{
	if(look_up(e)){								//edge found, return
		return;
	}
	
	vertex_descriptor_t v1 = add_vertex(e.first());				//get the handle of the vertex
	vertex_descriptor_t v2 = add_vertex(e.second());
	
	nbr_list_vec_[v1]->push_back(v2);
	nbr_list_vec_[v2]->push_back(v1);
}


FastGraph::edge_t FastGraph::add_edge(const vertex_descriptor_t& v1,
			 	      const vertex_descriptor_t& v2)
{
	edge_t e(v1,v2);
	if(!look_up(e)){
		add_vertex(v1);
		add_vertex(v2);
	}
	
	if(!is_neighbor(v1,v2)){
		nbr_list_vec_[v1]->push_back(v2);
		nbr_list_vec_[v2]->push_back(v1);
	}
	return e;
}

/**
	Remove an edge.
	@param e edge
*/

void FastGraph::remove_edge(const edge_t& e)
{
	vertex_descriptor_t v1 = e.first();
	vertex_descriptor_t v2 = e.second();
		
	nbr_list_descriptor_t list1 = nbr_list_vec_[v1];
	nbr_list_iterator_t it1(find(list1->begin(),
				     list1->end(),
				     v2)
			       );
	
	if(it1 != list1->end()){
		list1->erase(it1);
	}
	
	nbr_list_descriptor_t list2 = nbr_list_vec_[v2];
	
	nbr_list_iterator_t it2(find(list2->begin(),
				     list2->end(),
				     v1)
			       );
	
	if(it2 != list2->end()){
		list2->erase(it2);
	}
}


/**
	Remove an edge.
	@param e edge
*/

void FastGraph::remove_edge(const vertex_descriptor_t& v1,
			    const vertex_descriptor_t& v2)
{
	if(!look_up(v1))
		return;
	
	if(!look_up(v2))
		return;
	
	edge_t e(v1,v2);
	remove_edge(e);
}

bool FastGraph::is_neighbor(const vertex_descriptor_t& v1,
			    const vertex_descriptor_t& v2) const
{
	if(!look_up(v1) || !look_up(v2))
		return false;
	
	nbr_list_descriptor_t list = nbr_list_vec_[v1];
	nbr_list_const_iterator_t it(find(list->begin(),
					  list->end(),
					  v2)
				     );
	return it != list->end();
}

/**
	Checks if a vertex exists in the graph.
	@param v vertex
	@return pointer to vertex if found; NULL if not
*/

bool FastGraph::look_up(const vertex_descriptor_t& v) const
{
	return v < nbr_list_vec_.size();	
}

/**
	Checks if an edge exists in the graph.
	@param e edge
	@return pointer to edge if found; NULL if not
*/
bool FastGraph::look_up(const vertex_descriptor_t& v1,
			const vertex_descriptor_t& v2) const
{
	if(!look_up(v1) || !look_up(v2))
		return false;
	
	nbr_list_descriptor_t list = nbr_list_vec_[v1];
	nbr_list_const_iterator_t it(find(list->begin(),
					  list->end(),
					  v2)
				     );
	return it != list->end();
}


/**
	Checks if an edge exists in the graph.
	@param e edge
	@return pointer to edge if found; NULL if not
*/
bool FastGraph::look_up(const edge_t& e) const
{
	return look_up(e.first(),e.second());
}
FastGraph::vertex_list_t FastGraph::neighbors(const vertex_descriptor_t& v) const
{
	vertex_list_t list;
	if(!look_up(v))
		return list;
	return *(nbr_list_vec_[v]);
}

inline size_t FastGraph::number_of_vertices() const
{
	return nbr_list_vec_.size();
}


FastGraph::vertex_list_t& FastGraph::neighbors(const vertex_descriptor_t& v)
{
	return *(nbr_list_vec_[v]);
}

FastGraph::vertex_list_t FastGraph::neighbors(const vertex_descriptor_t& v,
					      const size_t& degree) const
{
	std::map<vertex_descriptor_t, size_t> dist;
	std::map<vertex_descriptor_t, size_t>::const_iterator iter;

	vertex_list_t result;
	
	if(!look_up(v) || is_isolated(v))
		return result;
	
	vertex_list_t Q;
	for(size_t i = 0; i < nbr_list_vec_.size(); i++){
		if(!is_isolated(i))
			Q.push_back(i);
	}
	
	dist[v] = 0;
	
	vertex_list_t::const_iterator itv;
	
	vertex_descriptor_t u;
	
	while(Q.size() > 0){
		u = min_distance(Q,dist);
		if( !(dist[u] <= degree) )
			break;

		Q.remove(u);
		result.push_back(u);	//all removed vertices are within range
		const vertex_list_t& nbr = neighbors(u);
		for(itv = nbr.begin(); itv != nbr.end(); itv++){
			edge_t e(u,*itv);
			const size_t w = 1;
			size_t alt = dist[u] + w;
			if(dist.find(*itv) == dist.end() || alt < dist[*itv])
				dist[*itv] = alt;
		}
		
	}

	result.remove(v);
	
	return result;
}


size_t FastGraph::number_of_neighbors(const vertex_descriptor_t& v) const
{
	if(!look_up(v))
		return 0;
	return nbr_list_vec_[v]->size();
}


/**
	Find the vertex that has the minimum distance which is stored in the edge weight map.
	@param Q vertex list
	@param dist edge weight map
	@return vertex with the smallest distance
*/
FastGraph::vertex_descriptor_t FastGraph::min_distance(const vertex_list_t& Q,
						       const std::map<vertex_descriptor_t,
						       size_t>& dist) const
{
	std::map<vertex_descriptor_t, size_t>::const_iterator itd;
	vertex_list_t::const_iterator itq;
	size_t d = numeric_limits<size_t>::max();
	vertex_descriptor_t u = Q.front();
	
	
	for(itq = Q.begin(); itq != Q.end(); itq++){
		itd = dist.find(*itq);
		if(itd == dist.end())
			continue;

		if(itd->second < d){
			d = itd->second;
			u = itd->first;
		}
	}
	return u;
}

ostream& operator<< (ostream& os, FastGraph& rhs)
{
	for(size_t i = 0; i < rhs.nbr_list_vec_.size(); i++){
		os << i << " : ";
		FastGraph::nbr_list_descriptor_t list = rhs.nbr_list_vec_[i];
		FastGraph::nbr_list_const_iterator_t ite;
		copy(list->begin(),
		     list->end(),
		     ostream_iterator<FastGraph::vertex_descriptor_t>(os, " ")
		     );
		os << endl;
	}
	
	return os;
}


}	//namespace bmpg_uncc_edu::algorithms
}	//namespace bmpg_uncc_edu

#endif

