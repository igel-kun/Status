
#include "seq.hpp"

namespace status {

  struct status_helper {
    uint paths_from_subtree;
    uint paths_to_subtree;
    uint subtree_size;
    uint status;
  };

  // compute the subtree sizes of each vertex in the subtree
  void compute_subtree_sizes(vertex* const v){
    assert(v != NULL);
    assert(v->get_data() == NULL);
    uint* const data((uint*)calloc(1, sizeof(uint)));
    const list<vertex*>& children(v->get_children());

    // recurse for all children u of v
    for(auto u = children.begin(); u != children.end(); ++u){
      // compute the values for u
      compute_subtree_sizes(*u);
      // and use them to update v's data
      *data += *((uint*)(*u)->get_data());
    }
    // add myself to the subtree_size
    ++(*data);
    v->set_data(data);
  }




  // compute the sum of lengths of paths in v's subtree going to v and write it to v's data (and return it)
  void compute_paths_from_subtree(vertex* const v){
    assert(v != NULL);
    assert(v->get_data() == NULL);
    status_helper* const data((status_helper*)calloc(1, sizeof(status_helper)));
    const list<vertex*>& children(v->get_children());

    // recurse for all children u of v
    for(auto u = children.begin(); u != children.end(); ++u){
      // compute the values for u
      compute_paths_from_subtree(*u);
      status_helper* const u_data((status_helper*)(*u)->get_data());
      // and use them to update v's data
      data->subtree_size += u_data->subtree_size;
      // each path in the subtree of u is one longer in the subtree of v
      data->paths_from_subtree += u_data->paths_from_subtree + u_data->subtree_size;
    }
    // add myself to the subtree_size
    data->subtree_size++;

    v->set_data(data);
  }

  // compute the sum of lengths of paths outside v's subtree going to v, add it to v's data (and return it)
  void compute_paths_to_subtree(vertex* const v, const uint tree_size){
    assert(v != NULL);
    status_helper* data((status_helper*)(v->get_data()));
    const list<vertex*>& children(v->get_children());
    // set the status of v to the total incoming paths of v
    data->status = data->paths_from_subtree + data->paths_to_subtree;

    // recurse for all children u of v
    for(auto u = children.begin(); u != children.end(); ++u){
      status_helper* const u_data((status_helper*)(*u)->get_data());
      // the size of paths to the subtree of u is computed from 3 parts: it is
      //    1. the status of v
      //    2. minus the paths from the subtree of u to v
      //    3. plus one for each path to v that is not in the subtree of u
      u_data->paths_to_subtree = data->status
                                - ( u_data->paths_from_subtree + u_data->subtree_size )
                                + ( tree_size - u_data->subtree_size );

      // recurse for u
      compute_paths_to_subtree(*u, tree_size);
    }
  }
  // read the previously computed stati
  void read_status_sequence(const vertex* const v, sequence_t& seq){
    assert(v != NULL);
    const list<vertex*>& children(v->get_children());

    seq[((status_helper*)(v->get_data()))->status]++;
    // recurse for all children u of v
    for(auto u = children.begin(); u != children.end(); ++u)
      read_status_sequence(*u, seq);
  }

  sequence_t compute_stati(const tree& t){
    sequence_t seq;
    if(! t.get_size()) return seq;

    // prepare the tree
    vertex* const root = t.get_root();
    DEBUG1(cout << "step1: paths from subtrees"<<endl);
    compute_paths_from_subtree(root);
    DEBUG1(cout << "step2: paths to subtrees"<<endl);
    compute_paths_to_subtree(root, t.get_size());
    DEBUG1(cout << "step3: stati"<<endl);
    read_status_sequence(root, seq);

    return seq;
  }


  // convert a status_sequence_t to a status_list
  list<uint> seq_to_list(const sequence_t& s){
    list<uint> result;
    for(auto i = s.begin(); i != s.end(); ++i)
      for(uint j = 0; j < i->second; ++j)
        result.push_back(i->first);
    return result;
  }

  list<uint> get_occuring_stati(const sequence_t& s){
    list<uint> result;
    for(auto i = s.begin(); i != s.end(); ++i)
      result.push_back(i->first);
    return result;
  }

  // return the number of vertices implied by the status sequence_t
  uint get_num_vertices(const sequence_t& s){
    uint result = 0;
    for(auto i = s.begin(); i != s.end(); ++i)
      result += i->second;
    return result;
  }

  // read a sequence_t from file
  sequence_t read_sequence_from_file(const string filename){
    ifstream f(filename);
    if(f.bad()) FAIL("unable to open "<<filename<<" for reading");
    
    sequence_t s;
    
    while(true){
      // read the next sequence
      string s_list;
      getline(f, s_list);
      if(f.eof()) break;

      list<string> items;

      // split by space (separating 1x38 and 5x32)
      split_no_empty(s_list, items);
      for(string i : items){
        list<string> components;
        // split by x (separating 1 and 38)
        split_no_empty(i, components, "x");
        // register in the sequence
        s[atoi(components.back().c_str())] = atoi(components.front().c_str());
      }
    }
    return s;
  }

  void write_sequence_to_file(const sequence_t& s, const string filename){
    ofstream f(filename);
    if(f.bad()) FAIL("unable to open "<<filename<<" for writing");
    f << s;
  }

  // compute one of the (at most two) vertices of minimum status
  vertex* compute_median_subtree(vertex* const v){
    const list<vertex*>& children(v->get_children());
    vertex* min_status = v;
    for(auto u = children.begin(); u != children.end(); ++u){
      vertex* tmp = compute_median_subtree(*u);
      // update min_status if the minimum status vertex in the subtree of u is smaller than the current min_status
      if(((status_helper*)tmp->get_data())->status < ((status_helper*)min_status->get_data())->status)
        min_status = tmp;
    }
    return min_status;
  }

  // compare two sequences
  bool equal(const sequence_t& s1, const sequence_t& s2){
    for(pair<uint,uint> entry: s1)
      if(s2.at(entry.first) != entry.second) return false;
    return true;
  }

};




ostream& operator<<(ostream& os, const status::sequence_t& s){
  for(auto i = s.begin(); i != s.end(); ++i)
    os << i->second << 'x' << i->first << ' ';
  return os << '\b';
}


