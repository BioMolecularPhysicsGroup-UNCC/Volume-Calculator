/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_algorithms_Fast_Graph_hpp
#define bmpg_uncc_edu_algorithms_Fast_Graph_hpp

#include <list>
#include <vector>
#include <map>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <bmpg_uncc_edu/algorithms/GraphTraits.hpp>

namespace bmpg_uncc_edu {
namespace algorithms {
using namespace std;

/**
	The FastGraph class defines a generic, undirectional graph contaning 
	vertices and edges. Properties of vertices and edges are not stored in the graph. Instead,
	it's up to the user to set up vertex property map and edge property maps.
*/

class FastGraph
{
public:
	typedef size_t vertex_descriptor_t;
	
	class edge
	{
	public:
		edge();
		edge(vertex_descriptor_t v1,
		     vertex_descriptor_t v2)
		{
			if(v1 < v2){
				first_ = v1, second_ = v2;
			} else {
				first_ = v2, second_ = v1;
			};
		}
		
		const vertex_descriptor_t first() const {return first_;}
		const vertex_descriptor_t second() const {return second_;}
		vertex_descriptor_t target(const vertex_descriptor_t& v) const
		{
			if(v == first_){
				return second_;
			} else if(v == second_) {
				return first_;
			} else {
				vertex_descriptor_t empty;
				return empty;
			}
		}
		
		bool has(const vertex_descriptor_t& v) const {return (first_ == v) || (second_ == v);}
		
		bool operator==(const edge& rhs) const {return (first_ == rhs.first_ && second_ == rhs.second_) || (second_ == rhs.first_ && first_ == rhs.second_);}
		bool operator<(const edge& rhs) const {return (first_ < rhs.first_) || (first_ == rhs.first_ && second_ < rhs.second_);}
		
		vertex_descriptor_t first_;
		vertex_descriptor_t second_;
	};
	
	typedef edge edge_t;
	
	typedef std::list<vertex_descriptor_t> vertex_list_t;
	
	typedef std::list<vertex_descriptor_t> nbr_list_t;
	//iterators
	typedef nbr_list_t::iterator nbr_list_iterator_t;
	typedef nbr_list_t::const_iterator nbr_list_const_iterator_t;
	
	typedef boost::shared_ptr<nbr_list_t> nbr_list_descriptor_t;
	typedef std::vector<nbr_list_descriptor_t> nbr_list_vec_t;
	
	typedef nbr_list_vec_t::iterator nbr_list_vec_iterator_t;
	typedef nbr_list_vec_t::const_iterator nbr_list_vec_const_iterator_t;
	
	
	FastGraph();
	
	void add_edge(const edge_t& e);
	edge_t add_edge(const vertex_descriptor_t& v1,
		        const vertex_descriptor_t& v2);
	
	void remove_edge(const edge_t& e);
	void remove_edge(const vertex_descriptor_t& v1,
			 const vertex_descriptor_t& v2);

	vertex_descriptor_t add_vertex(const vertex_descriptor_t& v);

	bool is_neighbor(const vertex_descriptor_t& v1,
			 const vertex_descriptor_t& v2) const;
	
	inline vertex_list_t neighbors(const vertex_descriptor_t& v) const;
	inline vertex_list_t& neighbors(const vertex_descriptor_t& v);
	
	inline vertex_list_t neighbors(const vertex_descriptor_t& v,
				       const size_t& degree) const;
	
	inline size_t number_of_vertices() const;
	size_t number_of_neighbors(const vertex_descriptor_t& v) const;
	inline bool is_isolated(const vertex_descriptor_t& v) const{ return number_of_neighbors(v) < 1;}
	
	friend ostream& operator<<(ostream& os, FastGraph& rhs);

protected:
	nbr_list_vec_t nbr_list_vec_;
	vertex_list_t vertex_list_;
private:
	
	vertex_descriptor_t min_distance(const vertex_list_t& Q,
					 const std::map<vertex_descriptor_t,
					 size_t>& dist) const;

	bool look_up(const vertex_descriptor_t& v) const;
	bool look_up(const edge_t& e) const;
	bool look_up(const vertex_descriptor_t& v1,
		     const vertex_descriptor_t& v2) const;
};


}	//namespace bmpg_uncc_edu::algorithms
}	//namespace bmpg_uncc_edu

#endif
