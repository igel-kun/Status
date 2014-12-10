
#include "util/graphs.hpp"
#include "util/seq.hpp"
#include "solv/options.hpp"
#include "solv/caterpillar.hpp"
#include "math.h"
#include <cstdlib> // for std::rand()
#include <unistd.h> // for getpid
#include <algorithm> // for sort()

void usage(const char* progname, std::ostream& o){
  o << "usage: " << progname << " ftree <file to read> [more opts]\t- read tree from file" << std::endl;
  o << "       " << progname << " scat <file to read> [more opts]\t- read sequence from file" << std::endl;
  o << "       " << progname << " rtree <#vertices> [more opts]\t- create random tree"<< std::endl;
  o << "       " << progname << " rcat <#vertices> [more opts]\t- create random caterpillar"<< std::endl;
  o << "       " << progname << " rcats <#vertices> [more opts]\t- create random sparse caterpillar"<< std::endl;
  o << "       " << progname << " rscat <#vertices> <avg multiplicity> [more opts]\t- create random sequence & assume it's a caterpillar"<< std::endl;
  exit(1);
}

status::tree* get_random_tree(const size_t num_vertices){
  status::vertex* vertex_nr[num_vertices];
  status::vertex* root = new status::vertex();
  status::tree* t = new status::tree(root);

  std::srand(time(NULL) * getpid());

  vertex_nr[0] = root;
  for(uint i = 1; i < num_vertices; ++i)
    vertex_nr[i] = t->add_vertex(vertex_nr[std::rand() % i]);
  
  return t;
}


#define percentage_backbone 30
status::tree* get_random_caterpillar(const size_t num_vertices){
  status::vertex* vertex_nr[num_vertices];
  status::vertex* root = new status::vertex();
  status::tree* t = new status::tree(root);

  std::srand(time(NULL) * getpid());

  // create backbone
  size_t num_backbone = uint(num_vertices*percentage_backbone/100) + (std::rand() % uint(num_vertices*(100-percentage_backbone)/100));
  std::cout << "selected "<<num_backbone<<" vertices for backbone" << std::endl;
  vertex_nr[0] = root;
  for(uint i = 1; i < num_backbone; ++i)
    vertex_nr[i] = t->add_vertex(vertex_nr[i-1]);

  // create leaves
  for(uint i = num_backbone; i < num_vertices; ++i)
    vertex_nr[i] = t->add_vertex(vertex_nr[std::rand() % num_backbone]);
  
  return t;
}

// get a random caterpillar whose backbone vertices have between 0 and 1 leaves
status::tree* get_random_sparse_caterpillar(const size_t num_vertices){
  status::vertex* vertex_nr[num_vertices];
  status::vertex* root = new status::vertex();
  status::tree* t = new status::tree(root);

  std::srand(time(NULL) * getpid());

  // create backbone using at least half the vertices
  size_t num_backbone = uint((num_vertices+1)/2) + (std::rand() % uint(num_vertices/3));
  std::cout << "selected "<<num_backbone<<" vertices for backbone" << std::endl;
  vertex_nr[0] = root;
  for(uint i = 1; i < num_backbone; ++i)
    vertex_nr[i] = t->add_vertex(vertex_nr[i-1]);
  // shuffle backbone vertices
  std::random_shuffle(&vertex_nr[0], &vertex_nr[num_backbone]);

  // give one leaf to each of the first backbone vertices (since leaves < n/2, some bb vertices don't get leaves)
  for(uint i = num_backbone; i < num_vertices; ++i)
    vertex_nr[i] = t->add_vertex(vertex_nr[i - num_backbone]);
  
  return t;
}

// get a random sequence of given length and average multiplicity
status::sequence_t get_random_sequence(const uint num_vertices, const double avg_multi){
  std::srand(time(NULL) * getpid());
  // create a caterpillar of num_vertices vertices, compute its status sequence
  // and change stati into existing stati until the avg multiplicity is reached
  status::tree* t = get_random_caterpillar(num_vertices);
  status::sequence_t s(compute_stati(*t));
  delete t;

  // compute center status
  status::sequence_t::iterator center = s.begin();
  for(status::sequence_t::iterator i = s.begin(); i != s.end(); ++i)
    if(i->first < center->first) center = i;

  while(num_vertices < avg_multi * s.size()){
    // get two coordinates to merge
    const uint x(std::rand() % s.size());
    const uint y(std::rand() % s.size());
    if(x == y) continue;
    // transfer all occurances of the x'th status to the y'th status
    status::sequence_t::iterator from = s.begin();
    status::sequence_t::iterator to = s.begin();
    for(uint i = 0; i < x; ++i) ++from;
    for(uint i = 0; i < y; ++i) ++to;
    // don't touch the center!
    if((from == center) || (to == center)) continue;
    //std::cout << "merging "<<from->first<<" and "<<to->first<<std::endl;
    to->second += from->second;
    s.erase(from);
    //std::cout << "got "<<s.size()<<" stati, avg multi is "<<(double)num_vertices/s.size()<<std::endl;
  }

  return s;
}



const std::pair<string, int> _requires_params[] = {
  { "ftree", 1 },
  { "scat", 1 },
  { "rtree",  1 },
  { "rcat",  1 },
  { "rcats",  1 },
  { "rscat", 2},
};
// global arguments with their parameters
std::map<string, std::vector<string> > arguments;

void parse_args(int argc, char** argv, status::solv_options& opts){
  int arg_ptr;
  std::map<std::string, int>  requires_params(std::begin(_requires_params), std::end(_requires_params));
  arg_ptr = 1;
  while(arg_ptr < argc){
    const std::string arg(argv[arg_ptr++]);
    // if the argument is not registered in requires_args, then exit with usage
    if(requires_params.find(arg) == requires_params.end()) usage(argv[0], std::cerr);
    // if there are not enough parameters for this argument
    if(argc < arg_ptr + requires_params[arg]) usage(argv[0], std::cerr);
    // otherwise fill the argument map
    std::vector<string> params(requires_params[arg]);
    for(int i = 0; i < requires_params[arg]; ++i)
      params[i] = argv[arg_ptr++];
    arguments.insert(make_pair(arg, params));
  }
}

int main(int argc, char** argv)
{
  status::solv_options opts;
  status::sequence_t s;
  status::tree* t = NULL;
  bool is_caterpillar = false;
  status::tree* result;

  // parse the arguments, filling 'arguments'
  parse_args(argc, argv, opts);

  if(arguments.find("rtree") != arguments.end()){
    // create a random tree
    t = get_random_tree(atoi(arguments["rtree"][0].c_str()));
  } else if(arguments.find("rcat") != arguments.end()) {
    // create a random caterpillar
    t = get_random_caterpillar(atoi(arguments["rcat"][0].c_str()));
    is_caterpillar = true;
  } else if(arguments.find("rcats") != arguments.end()) {
    // create a random sparse caterpillar
    t = get_random_sparse_caterpillar(atoi(arguments["rcats"][0].c_str()));
    is_caterpillar = true;
  } else if(arguments.find("rscat") != arguments.end()) {
    // create a random sparse caterpillar
    s = get_random_sequence(atoi(arguments["rscat"][0].c_str()), atof(arguments["rscat"][1].c_str()));
    is_caterpillar = true;
  } else if(arguments.find("ftree") != arguments.end()) {
    // create a graph from file
    t = new status::tree();
    t->read_from_file(arguments["ftree"][0]);
    is_caterpillar = detect_caterpillar(*t);
  } else if(arguments.find("scat") != arguments.end()) {
    // create a graph from file
    status::sequence_t s(status::read_sequence_from_file(arguments["scat"][0]));
    std::cout << "status sequence: " << s << std::endl;
  } else usage(argv[0], std::cerr);

  if(t){ // if we're given a tree, convert it to a sequence, while writing it down to .tree
    cout << "writing tree to .tmp for reference" << endl;
    t->write_to_file(".tree");

    std::cout << *t << endl;
    cout << "computing stati"<< endl;
    const status::sequence_t s(status::compute_stati(*t));
    // reroot t at its median
    cout << "computing median"<<endl;
    status::vertex* median = status::compute_median(*t);
    t->clear_data();
    cout << "rerooting tree at median" <<endl;
    t->reroot(median);

    std::cout << *t << endl;
  } 

  std::cout << "status sequence: " << s << std::endl;
  // write the status sequence to .sequence
  status::write_sequence_to_file(s, ".sequence");

  if(is_caterpillar){
    result = status::stati_to_caterpillar(s);
  } else {
    std::cout << "this is not a caterpillar..."<<std::endl;
    result = NULL;
  }
  if(result){
    std::cout << "reconstructed:" << std::endl << *result << std::endl << "largest list: "<<status::get_set_list_max()<<std::endl;
    status::sequence_t check(status::compute_stati(*result));
    std::cout << "recheck stati "<<check<<": "<<(status::equal(s, check) ? "match! Good job :)" : "!!! NO MATCH !!!")<<std::endl;
  } else {
    std::cout << "could not reconstruct the graph" << std::endl << "largest list: "<<status::get_set_list_max()<<std::endl;
  }
}
