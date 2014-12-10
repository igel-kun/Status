#include "graphs.hpp"
#include <unordered_map>
#include <sstream>

namespace status{


  struct printv_helper {
    uint x;
    uint leaves;
  };

  void printv_compute_leaves(vertex* const v){
    assert(v->get_data() == NULL);

    const list<vertex*>& children(v->get_children());
    printv_helper* data((printv_helper*)calloc(1, sizeof(printv_helper*)));

    // if v is not a leaf, then recurse, otherwise, just set leaves to 1
    if(!v->is_leaf()){
      // first, compute the leaves
      for(auto u = children.begin(); u != children.end(); ++u){
        printv_compute_leaves(*u);
        data->leaves += ((printv_helper*)(*u)->get_data())->leaves;
      }
    } else data->leaves = 1;
    v->set_data(data);
  }

  void printv_compute_offset(vertex* const v, uint X_offset = 0){
    assert(v->get_data() != NULL);
    
    const list<vertex*>& children(v->get_children());
    printv_helper* data((printv_helper*)v->get_data());
    
    // according to the leaves and the XY_offset, decide an (x,y)-position of v:

    if(!v->is_leaf()){
      uint x_max = 0;
      uint x_min = INT_MAX;
      // and recurse again
      for(auto u = children.begin(); u != children.end(); ++u){
        printv_compute_offset(*u, X_offset);

        printv_helper* const u_data((printv_helper*)(*u)->get_data());
        X_offset += u_data->leaves;

        // save the boundaries to compute data->x
        if(u_data->x > x_max) x_max = u_data->x;
        if(u_data->x < x_min) x_min = u_data->x;
      }
      //  1. put v in the center of the children
      data->x = (x_min + x_max) >> 1;
    } else data->x = X_offset;
  }

  void vertex::print_to_stream_vertically(ostream& os) {
    // 1. get infos
    printv_compute_leaves(this);
    printv_compute_offset(this);

    // 2. go through the tree with BFS
    list<vertex*> l1, l2;
    l1.push_back(this);

    list<vertex*>* current_layer(&l1);
    list<vertex*>* next_layer(&l2);

    do{
      // save a list of start and endpoints of dashes for the next row
      list<pair<uint, uint> > dashes;
      uint current_x = 0;
      while(!current_layer->empty()){
        const vertex* const u(current_layer->front());
        current_layer->pop_front();

        // print indent
        while(current_x < ((printv_helper*)u->get_data())->x) { os << ' '; ++current_x; }
        os << 'o';
        ++current_x;

        // if u has children, we'll draw some dashes below u
        const list<vertex*>& children(u->get_children());
        if(!u->is_leaf()){
          uint x_max = 0;
          uint x_min = INT_MAX;
          for(auto w: children){
            const printv_helper* const w_data((printv_helper*)(w->get_data()));

            // save the boundaries for dashes in the next line
            if(w_data->x > x_max) x_max = w_data->x;
            if(w_data->x < x_min) x_min = w_data->x;

            // save the children for processing in next layer
            next_layer->push_back(w);
          }
          dashes.push_back(make_pair(x_min,x_max));
        }
      }
      // draw dashes
      os << endl;
      uint x_coord = 0;
      while(!dashes.empty()){
        pair<uint, uint> dash(dashes.front());
        dashes.pop_front();
        // print indent
        while(x_coord < dash.first) { os << ' '; ++x_coord; }
        // if there is just one child, draw a pipe instead of a dash
        if(dash.first == dash.second) {
          os << '|';
          ++x_coord;
        } else {
          os << '/'; ++x_coord;
          while(x_coord <= dash.second - 1) { os << '-'; ++x_coord; }
          os << '\\'; ++x_coord;
        }
      }
      os << endl;

      // exchange current and next layers
      list<vertex*>* tmp(next_layer);
      next_layer = current_layer;
      current_layer = tmp;
      // if the new current layer is empty, then break the loop
    } while(!current_layer->empty());

    // like my mom used to say: clean up after you!
    clear_data_subtree();
  }


  
  void tree::reroot(vertex* const new_root){
    assert(new_root != NULL);
    // make sure we're not already at the correct root
    vertex* const parent(new_root->get_parent());
    if(parent){
      // make the parent new root
      reroot(parent);
      // and change places with the parent
      parent->remove_child(new_root);
      parent->set_parent(new_root);
      new_root->add_child(parent);
      new_root->set_parent(NULL);
      root = new_root;
    }
  }

  // copy a tree and give a map resolving the old vertices to the new vertices
  tree* copy_tree_preserving(const tree& t, status::vertex_translator& translator){
    tree* result = new tree();
    translator.clear();

    // if we have the empty graph, just return the empty copy
    if(t.empty()) return result;

    // maintain a list of vertices to consider
    list<const vertex*> to_consider;
    to_consider.push_back(t.get_root());
    while(!to_consider.empty()){
      // get the first vertex
      const vertex* next(to_consider.front());
      to_consider.pop_front();

      // set the parent to NULL if we're considering the root and to the vertex corresponding to next's parent, otherwise
      vertex* parent(next->get_parent() ? translator.at(next->get_parent()) : NULL);

      // create the new vertex v in t
      vertex* v = result->add_vertex(parent, next->get_data());

      // relate next and v
      translator.insert(make_pair(next, v));

      // add next's children to to_consider
      for(vertex* i : next->get_children()) to_consider.push_back(i);
    }
    return result;
  }


  void tree::write_to_file(const string filename){
    list<vertex*> to_consider;
    ofstream f(filename);

    if(!root) return;
    if(f.bad()) FAIL("unable to open "<<filename<<" for writing");

    unordered_map<vertex*, uint> name;
    name[root] = 0;
    to_consider.push_back(root);
    uint last_used_name = 0;
    
    while(!to_consider.empty()){
      // get the first vertex
      vertex* v = to_consider.front();
      to_consider.pop_front();
      // add the edges to its children
      for(auto child: v->get_children()){
        // give the child a name
        name[child] = ++last_used_name;
        // print the edge to the outputfile
        f << name[v] << " " << name[child] << endl;
        // mark the child to_consider
        to_consider.push_back(child);
      }
    }
  }
  void tree::read_from_file(const string filename){
    ifstream f(filename);
    unordered_map<uint, vertex*> name;

    DEBUG2(cout << "reading tree from "<<filename<<endl);
    clear();
    if(f.bad()) FAIL("unable to open "<<filename<<" for reading");
    
    while(true){
      // read the next edge
      string next_edge;
      getline(f, next_edge);
      if(f.eof()) break;

      istringstream parse_edge(next_edge);
      uint u, v;
      parse_edge >> u;
      parse_edge >> v;
      // if we have not seen u yet, it's the root
      if(name.find(u) == name.end())
        name[u] = add_vertex(NULL);
      // we cannot have seen v before, so add his as a child of u
      name[v] = add_vertex(name[u]);
    }
  }


  // return whether the subtrees rooted at u and v are isomorph
  bool subtrees_isomorph(const vertex* const u, const vertex* const v){
    const list<vertex*>& u_childs(u->get_children());
    const list<vertex*>& v_childs(v->get_children());
    
    if(u_childs.size() == v_childs.size()){
//      auto ui = u_childs.begin();
//      TODO
    } else return false;
    return true;
  }

  void tree::merge_automorphism_classes(){
  }

  // return a list of non-leaf children
  inline const list<vertex*> vertex::get_non_leaf_children() const{
    list<vertex*> l;
    for(vertex* c : children) if(!c->is_leaf()) l.push_back(c);
    return l;
  }

  // return true iff replacing v's parent by a P2 still yields a caterpillar
  bool subtree_is_caterpillar(const vertex* v){
    assert(v != NULL);
    list<vertex*> l = v->get_non_leaf_children();
    if(l.empty()) return true;
    if(l.size() > 1) return false;
    return subtree_is_caterpillar(l.front());
  }
  // return true iff t is a caterpillar
  bool detect_caterpillar(const tree& t){
    // an empty tree is a caterpillar
    if(t.empty()) return true;
    // get the non-leaf neighbors of the root
    list<vertex*> l = t.get_root()->get_non_leaf_children();
    if(l.empty()) return true;
    if(l.size() > 2) return false;
    bool tmp = subtree_is_caterpillar(l.front());
    bool tmp2 = subtree_is_caterpillar(l.back());
    return tmp && tmp2;
  }

};

ostream& operator<<(ostream& os, const status::tree& t){
//    t.root->print_to_stream(os);
  if(!t.empty())
    t.root->print_to_stream_vertically(os);
  return os << "(size = "<<t.size<<")";
}



