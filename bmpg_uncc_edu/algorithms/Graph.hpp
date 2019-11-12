/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_algorithms_Graph_hpp
#define bmpg_uncc_edu_algorithms_Graph_hpp

#include <string>
#include <list>
#include <map>
#include <boost/shared_ptr.hpp>
#include <bmpg_uncc_edu/algorithms/GraphTraits.hpp>

namespace bmpg_uncc_edu {
namespace algorithms {
using namespace std;

template<typename Vertex, typename EdgeWeight> class Graph;
template<typename Vertex, typename EdgeWeight> std::ostream& operator<<(std::ostream& os, const Graph<Vertex,EdgeWeight>& rhs);

template<typename EdgeWeight> std::ostream& operator<<(std::ostream& os, const Graph<int,EdgeWeight>& rhs);
template<typename EdgeWeight> std::ostream& operator<<(std::ostream& os, const Graph<string,EdgeWeight>& rhs);

/**
	The Graph class is an undirectional graph that defines a generic graph containing vertices and edges. The vertices in the graph are labeled by
	vertex descriptors, where the data type of the descriptors depends on the data type of the vertex itself. The edges in the
	graph can be associated to a weight. Certain properties are required for both the vertex type and the weight type.
	
	The only requirement for vertex descriptor is to implement "compare equal", or the operator==(). The vertex of a 
	graph can be primitive,
	or non-primative type. A primative type, as specified in GraphTraits class, is an integer, char, std::string, or
	a pointer to an object. A data type that is not primative is non-primative. The descriptor of vertex is determined by
	the type of the vertex. If the vertex type is primative, then the descriptor type is the same as the vertex type; e.g., if
	the vertex type is std::string, then the descriptor of the vertex is also string. If the vertex type is non-primative, then
	the vertex descriptor is a boost::shared_ptr to the vertex type. Hence the object associated with the vertex needs to be pointed
	to by a boost::shared_ptr in ordered to be stored in the graph. Primative vertex descriptors are supposedly more efficient, because
	copying of the descriptors frequently occur in graph operations.
	
	As an example, the following lines are valid:
	
	Graph<string,int> g0;
	g0.add_vertex("Charlotte");
	g0.add_edge("Charlotte","Atlanta",3);
	
	or,
	
	class A;
	Graph<A,int> g1;
	boost::shared_ptr<A> a1(new A);
	g1.add_vertex(A);
	
	or,
	
	class A;
	Graph<A*,double> g2;
	A* a1 = new A;
	g2.add_vertex(a1);
	
	while it is not valid to say
	
	class A;
	Graph<A,int> g3;
	A* a1 = new A;
	g3.add_vertex(a1);
	
	because the descriptor of g3 is boost::shared_ptr, not A*. Note that sometimes a regular pointer is more desirable because
	of efficiency concerns, the user
	is responsible for releasing the memory of the object, or making sure that the pointers that are stored in the graph are not
	defunct (dangling pointers).
	
	The weight of an edge is required to be have the following operations: operator+, operator<, infinity function, and zero function.
	@see GraphWeight.hpp

*/


template<typename vertex_descriptor_type>
class edge
{
public:
	edge(vertex_descriptor_type v1,
	     vertex_descriptor_type v2) :
	     first_(v1),
	     second_(v2)
	     {
	     }
	     
	const vertex_descriptor_type first() const {return first_;}
	const vertex_descriptor_type second() const {return second_;}
	
	vertex_descriptor_type target(const vertex_descriptor_type& v) const
	{
		if(v == first_){
			return second_;
		} else if(v == second_) {
			return first_;
		} else {
			throw "edge::target: vertex does not belong to the edge.";
		}
	}
	
	bool has(const vertex_descriptor_type& v) const {return (first_ == v) || (second_ == v);}
	
	bool operator==(const edge& rhs) {return (first_ == rhs.first_ && second_ == rhs.second_) || (second_ == rhs.first_ && first_ == rhs.second_);}
	
	vertex_descriptor_type first_;
	vertex_descriptor_type second_;
};

template<class Vertex, class EdgeWeight = int>					//default weight is integer
class Graph
{
public:
	typedef EdgeWeight edge_weight_t;
	typedef Vertex vertex_t;
	typedef typename if_<vertex_t, boost::shared_ptr<vertex_t> >::type vertex_descriptor_t;
	
	typedef edge<vertex_descriptor_t> edge_t;
	typedef boost::shared_ptr<edge_t> edge_descriptor_t;
	
	typedef typename std::map<edge_descriptor_t,edge_weight_t> edge_weight_map_t;
	typedef typename std::map<vertex_descriptor_t,edge_weight_t> vertex_weight_map_t;
	
	typedef typename std::list<edge_descriptor_t> edge_list_t;
	typedef typename std::list<vertex_descriptor_t> vertex_list_t;
	
	typedef typename std::map<vertex_descriptor_t, edge_list_t> edge_list_map_t;
	typedef typename edge_list_map_t::iterator edge_list_map_iterator_t;
	typedef typename edge_list_map_t::const_iterator edge_list_map_const_iterator_t;
	
	typedef typename edge_list_t::iterator edge_list_iterator_t;
	typedef typename vertex_list_t::iterator vertex_iterator_t;
	typedef typename edge_list_t::const_iterator edge_const_iterator_t;
	typedef typename vertex_list_t::const_iterator vertex_const_iterator_t;
	
	typedef typename std::pair<edge_list_iterator_t,edge_list_iterator_t> edge_list_iterator_pair_t;	
	
	Graph();
	explicit Graph(const edge_list_t& edges);
	
	edge_weight_t Dijkstra(const vertex_descriptor_t& v1,
			       const vertex_descriptor_t& v2);
	
	vertex_weight_map_t Dijkstra(const vertex_descriptor_t& v);
	
	void add_edge(const edge_descriptor_t& e,
		      const edge_weight_t& w);
	
	void add_edge(const vertex_descriptor_t& v1,
		      const vertex_descriptor_t& v2, 
		      const edge_weight_t& w);
	
	void add_edge(const vertex_descriptor_t& v1,
		      const vertex_descriptor_t& v2);
	
	void remove_edge(const edge_descriptor_t& e);
	
	void remove_edge(const vertex_descriptor_t& v1,
			 const vertex_descriptor_t& v2);
	
	bool has_edge(const vertex_descriptor_t& v1,
		      const vertex_descriptor_t& v2) const;
	
	vertex_descriptor_t add_vertex(const vertex_descriptor_t& v);
	void remove_vertex(const vertex_descriptor_t& v);
	bool has(const vertex_descriptor_t& v) const {return look_up(v);}
	
	void add_weight(const edge_descriptor_t& e,
			const edge_weight_t& w);
	
	void add_weight(const vertex_descriptor_t& v1,
			const vertex_descriptor_t& v2, 
			const edge_weight_t& w);
	
	size_t num_of_nonisolated_vertices() const;
	size_t num_of_vertices() const;

	inline vertex_list_t vertices() const;
	inline edge_list_t edges(const vertex_descriptor_t& v) const;
	inline edge_list_t edges() const;
	
	inline vertex_list_t neighbors(const vertex_descriptor_t& v) const;
	inline vertex_list_t neighbors(const vertex_descriptor_t& v,
				       const edge_weight_t& degree) const;
	
	inline bool is_neighbor(const vertex_descriptor_t& v1,
				const vertex_descriptor_t& v2) const;
	
	inline size_t number_of_neighbors(const vertex_descriptor_t& v) const;
	inline double weight(const edge_weight_t& w);	
	
	Graph<Vertex, EdgeWeight>& operator=(const Graph<Vertex,EdgeWeight>& rhs);

	template<class EdgeWeight2>
	Graph<Vertex, EdgeWeight>& operator=(const Graph<Vertex,EdgeWeight2>& rhs);

	friend ostream& operator<< <>(ostream& os, const Graph<Vertex,EdgeWeight>& rhs);
	friend ostream& operator<< <>(ostream& os, const Graph<int,EdgeWeight>& rhs);
	friend ostream& operator<< <>(ostream& os, const Graph<string,EdgeWeight>& rhs);

	edge_list_map_t edge_list_map() const {return edge_list_map_;}
	
protected:
	edge_list_map_t edge_list_map_;
	edge_weight_map_t edge_weight_map_;
	
private:
	edge_descriptor_t remove(edge_list_t& list,
				 const vertex_descriptor_t& v);
	
	bool look_up(const vertex_descriptor_t& v) const;
	bool look_up(const edge_descriptor_t& e) const;
	edge_descriptor_t look_up(const vertex_descriptor_t& v1,
				  const vertex_descriptor_t& v2) const;
	
	vertex_descriptor_t min_distance(const vertex_list_t& Q,
					 const vertex_weight_map_t& dist) const;
};


template<class Vertex, class EdgeWeight>
Graph<Vertex,EdgeWeight>::Graph()
{
	
}


template<class Vertex, class EdgeWeight>
Graph<Vertex,EdgeWeight>::Graph(const edge_list_t& edge_list)
{
	edge_const_iterator_t it;
	for(it = edge_list.begin(); it != edge_list.end(); it++){
		add_edge(**it,edge_weight_t());
	}
}

template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::edge_weight_t Graph<Vertex,EdgeWeight>::Dijkstra(const vertex_descriptor_t& v1,
										    const vertex_descriptor_t& v2)
{
	vertex_weight_map_t dist;
	
	Limit<edge_weight_t> lim;
	edge_weight_t inf = lim.inf();
	edge_list_map_iterator_t it;
	vertex_list_t Q;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		Q.push_back(it->first);
	}

	//The initial mapping does not contain all vertices, only vertices are visited will be stored in the map to save time.
	dist[v1] = lim.zero();
	
	vertex_const_iterator_t itv;
	
	vertex_descriptor_t u = v1;
	edge_weight_t d = inf;
	
	while(Q.size() > 0){
		u = min_distance(Q,dist);
		Q.remove(u);
		if(u == v2)
			return dist[v2];
		
		const vertex_list_t& nbr = neighbors(u);
		for(itv = nbr.begin(); itv != nbr.end(); itv++){
			edge_weight_t& w = edge_weight_map_[look_up(u,*itv)];
			edge_weight_t alt = dist[u] + w;

			if(dist.find(*itv) == dist.end() || alt < dist[*itv])
				dist[*itv] = alt;
		}
		
	}
	return inf;
}

template<typename Vertex, typename EdgeWeight>
typename Graph<Vertex,EdgeWeight>::vertex_weight_map_t Graph<Vertex,EdgeWeight>::Dijkstra(const vertex_descriptor_t& v)
{
	vertex_weight_map_t dist;
	
	Limit<edge_weight_t> lim;
	
	edge_list_map_iterator_t it;
	vertex_list_t Q;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		Q.push_back(it->first);
	}
	dist[v] = lim.zero();
	
	vertex_const_iterator_t itv;
	
	vertex_descriptor_t u;
	
	while(Q.size() > 0){
		u = min_distance(Q,dist);
		Q.remove(u);
		const vertex_list_t& nbr = neighbors(u);
		for(itv = nbr.begin(); itv != nbr.end(); itv++){
			edge_weight_t& w = edge_weight_map_[look_up(u,*itv)];
			edge_weight_t alt = dist[u] + w;
			if(dist.find(*itv) == dist.end() || alt < dist[*itv])
				dist[*itv] = alt;
		}
		
	}
	return dist;
}

/**
	Find the vertex that has the minimum distance which is stored in the edge weight map.
	@param Q vertex list
	@param dist edge weight map
	@return vertex with the smallest distance
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::vertex_descriptor_t Graph<Vertex,EdgeWeight>::min_distance(const vertex_list_t& Q,
											      const vertex_weight_map_t& dist) const
{
	typename vertex_weight_map_t::const_iterator itd;
	vertex_const_iterator_t itq;
	Limit<edge_weight_t> lim;
	edge_weight_t d = lim.inf();
	vertex_descriptor_t u;
	
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

/**
	Get the weight.
	@param w weight object
	@return weight
*/
template<class Vertex, class EdgeWeight>
double Graph<Vertex,EdgeWeight>::weight(const edge_weight_t& w)
{
	return w();
}



/**
	Add a vertex to the graph. The vertex is duplicated if it does not exist in the graph.
	@param v vertex
	@return pointer to the duplicated vertex in the graph
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::vertex_descriptor_t Graph<Vertex,EdgeWeight>::add_vertex(const vertex_descriptor_t& v)
{
	if(look_up(v))
		return v;
	
	edge_list_t list;							//create an empty list
	edge_list_map_.insert(std::make_pair(v,list));
	return v;
}

/**
	Add an edge and the corresponding weight to the graph.
	The edge and the vertices corresponding to the edge are created.
	@param e edge
	@param w weight
*/
template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::add_edge(const edge_descriptor_t& e,
					const edge_weight_t& w)
{
	if(look_up(e)){								//already has edge
		edge_weight_map_[e] = w;					//store the weight
		return;
	}
		
	vertex_descriptor_t v1 = add_vertex(e->first());			//get the handle of the vertex
	vertex_descriptor_t v2 = add_vertex(e->second());
	
	edge_list_map_[v1].push_back(e);
	edge_list_map_[v2].push_back(e);
	
	edge_weight_map_.insert(std::make_pair(e,w));				//store the weight
}


/**
	Add an edge and the corresponding weight to the graph.
	The edge and the vertices corresponding to the edge are created.
	@param e edge
	@param w weight
*/
template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::add_weight(const edge_descriptor_t& e,
					  const edge_weight_t& w)
{
	if(look_up(e)){								//already has edge
		edge_weight_t& w0 = edge_weight_map_[e];
		w0 = w0 + w;							//add the weight
		return;
	}
	
	vertex_descriptor_t v1 = add_vertex(e->first());			//get the handle of the vertex
	vertex_descriptor_t v2 = add_vertex(e->second());
	
	edge_list_map_[v1].push_back(e);
	edge_list_map_[v2].push_back(e);
	
	edge_weight_map_.insert(std::make_pair(e,w));				//store the weight
}


/**
	Add a certain weight to an existing edge. A new edge with t the specified weight 
	will be created if no edge exists between the two vertices.
*/
template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::add_edge(const vertex_descriptor_t& v1,
					const vertex_descriptor_t& v2,
					const edge_weight_t& w)
{
	edge_descriptor_t e = look_up(v1,v2);
	edge_descriptor_t empty;
	if(e != empty){
		edge_weight_map_[e] = w;
		return;
	}
	
	vertex_descriptor_t pv1 = add_vertex(v1);				//get the handle of the vertex
	vertex_descriptor_t pv2 = add_vertex(v2);
	edge_descriptor_t edge_ptr(new edge_t(v1,v2));				//make a new copy of the edge
	edge_list_map_[pv1].push_back(edge_ptr);
	edge_list_map_[pv2].push_back(edge_ptr);
	
	edge_weight_map_.insert(std::make_pair(edge_ptr,w));			//store the weight	
}

/**
	Add a certain weight to an existing edge. A new edge with t the specified weight 
	will be created if no edge exists between the two vertices.
*/
template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::add_weight(const vertex_descriptor_t& v1,
					  const vertex_descriptor_t& v2, 
					  const edge_weight_t& w)
{
	edge_descriptor_t e = look_up(v1,v2);
	edge_descriptor_t empty;
	if(e != empty){
		edge_weight_t& w0 = edge_weight_map_[e];
		w0 = w0 + w;							//add the weight
		return;
	}
	
	vertex_descriptor_t pv1 = add_vertex(v1);				//get the handle of the vertex
	vertex_descriptor_t pv2 = add_vertex(v2);
	edge_descriptor_t edge_ptr(new edge_t(v1,v2));				//make a new copy of the edge
	edge_list_map_[pv1].push_back(edge_ptr);
	edge_list_map_[pv2].push_back(edge_ptr);
	
	edge_weight_map_.insert(std::make_pair(edge_ptr,w));			//store the weight	
}

template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::add_edge(const vertex_descriptor_t& v1,
					const vertex_descriptor_t& v2)
{
	edge_weight_t w;
	add_edge(v1,v2,w);
}

/**
	Remove an edge.
	@param e edge
*/
template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::remove_edge(const edge_descriptor_t& e)
{
	remove_edge(e->first(),e->second());
}


/**
	Remove an edge.
	@param e edge
*/
template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::remove_edge(const vertex_descriptor_t& v1,
					   const vertex_descriptor_t& v2)
{
	if(!look_up(v1))
		return;
	
	if(!look_up(v2))
		return;
	
	edge_list_t& list1 = edge_list_map_[v1];
	remove(list1,v2);
	
	edge_list_t& list2 = edge_list_map_[v2];
	edge_descriptor_t e = remove(list2,v1);
	edge_weight_map_.erase(e);
}


/**
	Remove a vertex. All the edges connected to the vertex are removed as well.
	@param v vertex
*/
template<class Vertex, class EdgeWeight>
void Graph<Vertex,EdgeWeight>::remove_vertex(const vertex_descriptor_t& v)
{
	if(!look_up(v))
		return;
	
	edge_list_t list = edge_list_map_[v];
	
	edge_list_iterator_t it;
	for(it = list.begin(); it != list.end(); it++){
		edge_descriptor_t e = *it;
		remove_edge(e);
	}
}
	
/**
	Remove all edges that are connected to a vertex from the list.
	@param list edge list
	@param v vertex
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::edge_descriptor_t Graph<Vertex,EdgeWeight>::remove(edge_list_t& list,
										      const vertex_descriptor_t& v)
{
	edge_descriptor_t result;
	
	edge_list_iterator_t it;
	for(it = list.begin(); it != list.end(); it++){
		if((*it)->has(v)){
			result = *it;
			list.erase(it);
			return result;
		}
	}
	edge_descriptor_t empty;
	return empty;
}


/**
	Checks if a vertex exists in the graph.
	@param v vertex
	@return pointer to vertex if found; NULL if not
*/
template<class Vertex, class EdgeWeight>
bool Graph<Vertex,EdgeWeight>::look_up(const vertex_descriptor_t& v) const
{
	
	edge_list_map_const_iterator_t it = edge_list_map_.find(v);
	return (it != edge_list_map_.end());
}


/**
	Checks if an edge exists in the graph.
	@param e edge
	@return pointer to edge if found; NULL if not
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::edge_descriptor_t Graph<Vertex,EdgeWeight>::look_up(const vertex_descriptor_t& v1,
										       const vertex_descriptor_t& v2) const
{
	edge_list_map_const_iterator_t it;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		const vertex_descriptor_t v = it->first;
		if( (v != v1) && (v != v2) ){
			continue;
		}
		edge_const_iterator_t ite;
		
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			if( (*ite)->has(v1) && (*ite)->has(v2) )
				return *ite;
		}
	}
	edge_descriptor_t empty;
	return empty;
	
}



/**
	Checks if an edge exists in the graph.
	@param e edge
	@return pointer to edge if found; NULL if not
*/
template<class Vertex, class EdgeWeight>
bool Graph<Vertex,EdgeWeight>::look_up(const edge_descriptor_t& e) const
{
	edge_list_map_const_iterator_t it;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		const vertex_descriptor_t v = it->first;
		if(!e->has(v)){
			continue;
		}
		edge_const_iterator_t ite;
		
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			if(*ite == e)
				return true;
		}
	}
	return false;
	
}


/**
	Checks if an edge exists between two vertices in the graph.
	@param e edge
	@return pointer to edge if found; NULL if not
*/
template<class Vertex, class EdgeWeight>
bool Graph<Vertex,EdgeWeight>::has_edge(const vertex_descriptor_t& v1,
					const vertex_descriptor_t& v2) const
{
	edge_list_map_const_iterator_t it;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		const vertex_descriptor_t v = it->first;
		if( (v != v1) && (v != v2) ){
			continue;
		}
		edge_const_iterator_t ite;
		
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			if( (*ite)->has(v1) && (*ite)->has(v2) )
				return true;
		}
	}
	return false;
}

/**
	Get the list of edges associated with a vertex.
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::edge_list_t Graph<Vertex,EdgeWeight>::edges(const vertex_descriptor_t& v) const
{
	edge_list_map_const_iterator_t it = edge_list_map_.find(v);
	if(it == edge_list_map_.end())
		return edge_list_t();						//return an empty list if vertex is not found
	return it->second;
}

/**
	Get the list of edges of the graph.
*/
template<class Vertex, class EdgeWeight>
inline typename Graph<Vertex,EdgeWeight>::edge_list_t Graph<Vertex,EdgeWeight>::edges() const
{
	edge_list_t result;
	edge_list_map_const_iterator_t it;
	typename edge_list_t::const_iterator ite;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		const vertex_descriptor_t v = it->first;
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			const edge_descriptor_t& e = *ite;
			if(e->first() != v)continue;				//avoid double counting
			result.push_front(e);
		}
	}
	return result;
}

/**
	Count the number of non-isolated vertices in the graph.
*/
template<class Vertex, class EdgeWeight>
size_t Graph<Vertex,EdgeWeight>::num_of_nonisolated_vertices() const
{
	size_t count = 0;
	edge_list_map_const_iterator_t it;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		if(!it->second.empty())
			count++;
	}
	return count;
}


/**
	Count the number of vertices in the graph, including isolated vertices.
*/
template<class Vertex, class EdgeWeight>
size_t Graph<Vertex,EdgeWeight>::num_of_vertices() const
{
	size_t count = edge_list_map_.size();
	return count;
}


/**
	Get the list of vertices (including isolated ones).
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::vertex_list_t Graph<Vertex,EdgeWeight>::vertices() const
{
	vertex_list_t result;
	edge_list_map_const_iterator_t it;
	typename edge_list_t::const_iterator ite;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		const vertex_descriptor_t v = it->first;
		result.push_front(v);
	}
	return result;
}

/**
	Get the list of neighobrs of a vertex.
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::vertex_list_t Graph<Vertex,EdgeWeight>::neighbors(const vertex_descriptor_t& v) const
{
	vertex_list_t list;
	if(!look_up(v))
		return list;
	
	edge_list_map_const_iterator_t it_map = edge_list_map_.find(v);
	
	edge_list_t elist = it_map->second;
	edge_list_iterator_t it;
	for(it = elist.begin(); it != elist.end(); it++){
		list.push_back((*it)->target(v));
	}
	return list;
}

/**
	Get the list of neighbors within a certain distance (weight).
*/
template<class Vertex, class EdgeWeight>
typename Graph<Vertex,EdgeWeight>::vertex_list_t Graph<Vertex,EdgeWeight>::neighbors(const vertex_descriptor_t& v,
										     const edge_weight_t& degree) const
{vertex_weight_map_t dist;
	
	Limit<edge_weight_t> lim;
	edge_weight_t inf = lim.inf();
	
	edge_list_map_const_iterator_t it_map;
	vertex_list_t result;
	
	vertex_list_t Q;
	for(it_map = edge_list_map_.begin(); it_map != edge_list_map_.end(); it_map++){
		Q.push_back(it_map->first);
	}
	dist[v] = lim.zero();
	
	vertex_const_iterator_t itv;
	
	vertex_descriptor_t u;
	edge_list_map_const_iterator_t it;
	
	while(Q.size() > 0){
		u = min_distance(Q,dist);
		if( !(dist[u] < degree) && !(dist[u] == degree) )		//only < and == are defined
			break;
		
		Q.remove(u);
		result.push_back(u);						//all removed vertices are within range
		const vertex_list_t& nbr = neighbors(u);
		for(itv = nbr.begin(); itv != nbr.end(); itv++){
			edge_descriptor_t e = look_up(u,*itv);
			typename edge_weight_map_t::const_iterator iter = edge_weight_map_.find(e);
			const edge_weight_t& w = iter->second;
			edge_weight_t alt = dist[u] + w;
			if(dist.find(*itv) == dist.end() || alt < dist[*itv])
				dist[*itv] = alt;
		}
		
	}
	
	result.remove(v);
	
	return result;
}

/**
	Check if two vertices are neighbors.
*/
template<class Vertex, class EdgeWeight>
bool Graph<Vertex,EdgeWeight>::is_neighbor(const vertex_descriptor_t& v1,
					   const vertex_descriptor_t& v2) const
{
	edge_list_map_const_iterator_t it;
	for(it = edge_list_map_.begin(); it != edge_list_map_.end(); it++){
		const vertex_descriptor_t v = it->first;
		if( (v != v1) && (v != v2) ){
			continue;
		}
		edge_const_iterator_t ite;
		
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			if( (*ite)->has(v1) && (*ite)->has(v2) )
				return true;
		}
	}
	return false;
}

/**
	Count the number of neighbors of a vertex.
*/
template<class Vertex, class EdgeWeight>
size_t Graph<Vertex,EdgeWeight>::number_of_neighbors(const vertex_descriptor_t& v) const
{
	edge_list_map_const_iterator_t it_map = edge_list_map_.find(v);
	const edge_list_t& elist = it_map->second;
	return elist.size();
}

template<class Vertex, class EdgeWeight>	
Graph<Vertex, EdgeWeight>& Graph<Vertex, EdgeWeight>::operator=(const Graph<Vertex,EdgeWeight>& rhs)
{
	edge_list_map_ = rhs.edge_list_map_;
	edge_weight_map_ = rhs.edge_weight_map_;
	return *this;
}

template<class Vertex, class EdgeWeight>	
template<class EdgeWeight2>
Graph<Vertex, EdgeWeight>& Graph<Vertex, EdgeWeight>::operator=(const Graph<Vertex,EdgeWeight2>& rhs)
{
	edge_list_map_ = rhs.edge_list_map();
	return *this;
}

template<class Vertex, class EdgeWeight>
ostream& operator<< (ostream& os, const Graph<Vertex,EdgeWeight>& rhs)
{
	typename Graph<Vertex,EdgeWeight>::edge_list_map_const_iterator_t it;
	for(it = rhs.edge_list_map_.begin(); it != rhs.edge_list_map_.end(); it++){
		const typename Graph<Vertex,EdgeWeight>::vertex_descriptor_t v = it->first;
		typename Graph<Vertex,EdgeWeight>::edge_const_iterator_t ite;
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			if((*ite)->first() != v)continue;				//avoid double counting
			os << *((*ite)->first()) << "----" << *((*ite)->second()) << "\n";
		}
	}
	
	return os;
}

template<class EdgeWeight>
ostream& operator<< (ostream& os, const Graph<int,EdgeWeight>& rhs)
{
	typename Graph<int,EdgeWeight>::edge_list_map_const_iterator_t it;
	for(it = rhs.edge_list_map_.begin(); it != rhs.edge_list_map_.end(); it++){
		const typename Graph<int,EdgeWeight>::vertex_descriptor_t v = it->first;
		typename Graph<int,EdgeWeight>::edge_const_iterator_t ite;
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			if((*ite)->first() != v)continue;				//avoid double counting
			os << (*ite)->first() << "----" << (*ite)->second() << "\n";
		}
	}
	
	return os;
}


template<class EdgeWeight>
ostream& operator<< (ostream& os, const Graph<string,EdgeWeight>& rhs)
{
	typename Graph<string,EdgeWeight>::edge_list_map_const_iterator_t it;
	for(it = rhs.edge_list_map_.begin(); it != rhs.edge_list_map_.end(); it++){
		const typename Graph<string,EdgeWeight>::vertex_descriptor_t v = it->first;
		typename Graph<string,EdgeWeight>::edge_const_iterator_t ite;
		for(ite = it->second.begin(); ite != it->second.end(); ite++){
			if((*ite)->first() != v)continue;				//avoid double counting
			os << (*ite)->first() << "----" << (*ite)->second() << "\n";
		}
	}
	
	return os;
}

}	//namespace bmpg_uncc_edu::algorithms
}	//namespace bmpg_uncc_edu

#endif
