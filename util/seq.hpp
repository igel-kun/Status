#ifndef SEQ_HPP
#define SEQ_HPP

#include "defs.hpp"
#include "graphs.hpp"
#include <unordered_map>
#include <unordered_set>
#include <list>

namespace status {

  // a status sequence is represented by a mapping of stati to nr. of occurances of this status
  typedef unordered_map<uint, uint> sequence_t;
  typedef list<uint> status_list_t;
  typedef unordered_set<uint> status_set_t;

  // compute the status of each vertex
  sequence_t compute_stati(const tree& t);

  // convert a status_sequence_t to a status_list
  list<uint> seq_to_list(const sequence_t& s);
  
  // return a list of stati occuring in the sequence_t
  list<uint> get_occuring_stati(const sequence_t& s);

  // return the number of vertices implied by the status sequence_t
  uint get_num_vertices(const sequence_t& s);

  // read a sequence_t from file
  sequence_t read_sequence_from_file(const string filename);
  void write_sequence_to_file(const sequence_t& s, const string filename);

  // compute one of the (at most two) vertices of minimum status
  vertex* compute_median_subtree(vertex* const v);
  inline vertex* compute_median(const tree& t){ return compute_median_subtree(t.get_root()); }

  // compare two sequence_ts
  bool equal(const sequence_t& s1, const sequence_t& s2);
};

ostream& operator<<(ostream& os, const status::sequence_t& s);


#endif
