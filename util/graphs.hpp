#ifndef GRAPHS_HPP
#define GRAPHS_HPP

#include <set>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include <fstream>
#include <iostream>
#include <sstream>

#include "defs.hpp"


using namespace std;


namespace status{

  class vertex;
  class tree;
};

ostream& operator<<(ostream& os, const status::tree& t);

namespace status {
  class vertex {
  private:
    // data
    void* data;
    // infrastructure
    vertex* parent;
    list<vertex*> children;

  public:
    // ==================== constructors =========================
    vertex(vertex* _parent = NULL):data(NULL),parent(_parent),children(){}
    vertex(vertex* _parent, void* _data):data(_data),parent(_parent),children(){}
    ~vertex(){ for(auto c: children) delete c; }

    // ================== data interaction =======================
    // get/set the data
    inline vertex* const get_parent() const { return parent; }
    inline void set_parent(vertex* const _parent) { parent = _parent; }
    inline void* get_data() const { return data; }
    inline void set_data(void* _data) { assert(data == NULL); data = _data; }
    // clear the data in this vertex
    void clear_data() { if(data) { free(data); data = NULL; } }
    // clear the data in this subtree
    void clear_data_subtree() {
      clear_data();
      for(auto i = children.begin(); i != children.end(); ++i) (*i)->clear_data_subtree();
    }

    // ================== infrastructure =========================
    inline bool is_root() const  { return parent == NULL; }
    inline uint degree() const { return children.size() + (parent == NULL ? 0 : 1); }
    inline bool is_leaf() const  { return children.empty(); }

    // get the children
    inline const list<vertex*>& get_children() const { return children; }
    // return a list of non-leaf children
    inline const list<vertex*> get_non_leaf_children() const;
    inline void add_child(vertex* const child) { children.push_back(child); }
    inline void remove_child(const list<vertex*>::iterator i) { children.erase(i); }
    inline void remove_child(const vertex* const child) {
      // first, find the child
      auto w = children.begin();
      for(; *w != child; ++w) assert(w != children.end());
      // then remove it
      remove_child(w);
    }

    // create new child and append it to children
    vertex* add_child(void* data = NULL){
      vertex* v = new vertex(this, data);
      add_child(v);
      return v;
    }

    void print_to_stream(ostream& os, const uint depth = 0, const string prepend = "") const{
      os << "o";
      // if we're a leaf, just return
      if(is_leaf()) return;
      // otherwise, print the first subtree
      auto i = children.begin();
      os << "-";
      (*i)->print_to_stream(os, depth + 1, children.size() == 1 ? prepend + "  " : prepend + "| ");
      // for each of the following subtrees, prepend a newline and depth x 2 spaces
      uint count = 1;
      while(++i != children.end()){
        os << std::endl;
        os << prepend << "|-";
        // if there are more children to come, add "| ", otherwise, add "  " to the prepend string
        if(++count < children.size())
          (*i)->print_to_stream(os, depth + 1, prepend + "| ");
        else
          (*i)->print_to_stream(os, depth + 1, prepend + "  ");
      }
    }

    // decide where to put this vertex (x,y), put it in the 2D buffer and return (the number of leaves, (x,y))
    void print_to_stream_vertically(ostream& os);
  };


  typedef unordered_map<const vertex*, vertex*> vertex_translator;

  class tree {
  private:
    uint size;
    vertex* root;

  public:
    // ==================== constructors =========================
    tree(vertex* _root = NULL):size(_root?1:0),root(_root){}
    
    // ================== data interaction =======================
    // get the size of the tree
    inline uint get_size() const { return size; }
    // get the root
    inline vertex* get_root() const { return root; }
    // clear the custom data in the tree
    inline void clear_data() const { if(root) root->clear_data_subtree(); }
    // clear the vertices
    inline void clear_vertices() { if(root) delete root; }
    // clear all the tree
    inline void clear() { clear_data(); clear_vertices(); }
    // return true iff the tree is empty
    inline bool empty() const { return size == 0; }

    // ================== infrastructure =========================
    vertex* add_vertex(vertex* const parent, void* data = NULL){
      if(!parent){
        // if no parent is given, but we have a root, make the root the parent...
        if(root) return add_vertex(root, data);
        size = 1;
        root = new vertex(NULL, data);
        return root;
      } else {
        ++size;
        return parent->add_child(data);
      }
    }
    inline vertex* add_vertex(){ return add_vertex(NULL); }

    // reroot the tree at new_root
    void reroot(vertex* const new_root);

    // merge equivalence classes under automorphism
    void merge_automorphism_classes();

    void write_to_file(const string filename);
    void read_from_file(const string filename);

    friend ostream& ::operator<<(ostream& os, const status::tree& t);
  
  };


  // copy a tree while preserving a set of vertices,
  // that is, translator will point to the vertices in the new tree that it pointed to in the old tree
  tree* copy_tree_preserving(const tree& t, status::vertex_translator& translator);

  // return true iff t is a caterpillar
  bool detect_caterpillar(const tree& t);

};




#endif
