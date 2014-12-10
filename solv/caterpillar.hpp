#ifndef CATERPILLAR_HPP
#define CATERPILLAR_HPP

#include "../util/seq.hpp"
#include "../util/graphs.hpp"

namespace status{

  size_t get_set_list_max();

  struct cat_dynprog_input {
    uint status;
    uint subtree;
    uint influx;
  };

#define INVALID_INPUT ((cat_dynprog_input){0, UINT_MAX, 0})
#define NO_INPUT ((cat_dynprog_input){0, 0, 0})


  struct leaf_guess_t {
    unsigned char backbone;
    uint leaves;
  };

#define NEW_GUESS ((leaf_guess_t){(unsigned char)0, (uint)0})


  class input_hasher{
  public:
    // partition our uint into three parts and XOR them together
    uint operator()(const cat_dynprog_input& x, const uint bitwidth = (sizeof(uint) << 3)) const{
      return ( x.status << (bitwidth / 3) )
         ^   ( x.subtree << (bitwidth / 3) )
         ^   ( x.influx << (bitwidth / 3) );
    }
  };

  typedef pair<cat_dynprog_input, cat_dynprog_input> cat_dynprog_input_pair;

  class input_pair_hasher{
  public:
    // partition our uint into two parts and XOR them together
    input_hasher ih;
    uint operator()(const cat_dynprog_input_pair& x) const{
      return ih(x.first, sizeof(uint) << 2)
         ^   ih(x.second, sizeof(uint) << 2);
    }
  };

  // describes possible status-configurations tailing the backbone vertex adjacent to leaves of the current status
  // like so (o = adjacent to leaves of current status, O = tail configuration)
  //
  //  \|/
  // --o--O--O--O--O
  struct cat_dynprog_config {
    // the inner caterpillar with two docking vertices (far left & right on backbone)
    tree* inner_cat;
    pair<vertex*, vertex*> docks;
    // save which stati we used how often
    sequence_t stati_used;


    cat_dynprog_config():inner_cat(NULL),docks(NULL,NULL){}
    // TODO write desctructor!!
    // copy constructor, using copy constructor of status_set_t
    cat_dynprog_config(const cat_dynprog_config& config):stati_used(config.stati_used){
      // copy the tree and the docks
      unordered_map<const vertex*, vertex*> preserve;
      inner_cat = copy_tree_preserving(*config.inner_cat, preserve);
      if(config.docks.first) docks.first = preserve.at(config.docks.first);
      if(config.docks.second) docks.second = preserve.at(config.docks.second);
    }
    // equality means equality of used stati, we actually don't care about how the graph looks
    inline bool operator==(const cat_dynprog_config& conf) const {
      return equal(stati_used, conf.stati_used);
    }
  };

  class config_hasher{
  public:
    // config hasher: roll each status s to the left by x bits, where x is the #occurances of s mod sizeof(uint)
    uint operator()(const cat_dynprog_config& conf) const{
      uint result = 0;
      for(auto entry: conf.stati_used)
        result += rol(entry.first, (entry.second % sizeof(uint)));
      return result;
    }
  };


  
  typedef unordered_set<cat_dynprog_config, config_hasher> cat_dynprog_configs;

  bool operator==(const cat_dynprog_input& X, const cat_dynprog_input& Y);
  bool operator==(const cat_dynprog_input_pair& X, const cat_dynprog_input_pair& Y);


  // reconstruct a caterpillar from a given status sequence
  tree* stati_to_caterpillar(const sequence_t& s);

  // get the status of leaves attached to a backbone vertex with status s
  inline uint get_corresponding_leaf_status(const uint s, const uint num_vertices){
    return s + (num_vertices - 2);
  }
  // get the status of the backbone vertex attached to a leaf vertex with status s
  inline uint get_corresponding_backbone_status(const uint s, const uint num_vertices){
    return s - (num_vertices - 2);
  }


}

ostream& operator<<(ostream& os, const status::cat_dynprog_input& i);
ostream& operator<<(ostream& os, const status::cat_dynprog_config& c);
ostream& operator<<(ostream& os, const status::leaf_guess_t& g);

#endif

