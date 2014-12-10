#ifndef COMPACT_SETS_HPP
#define COMPACT_SETS_HPP

#include "defs.hpp"
#include <unordered_set>
#include <list>


using namespace std;

namespace status { class compact_multiset; }

ostream& operator<<(ostream& os, const status::compact_multiset& s);

namespace status{


  // a compact description of a list of sets in two parts:
  // 1. a set "fix"
  // 2. a set "variable"
  // the list then contains exactly the sets obtained by combining "fix" with each subset of "variable"
  class compact_multiset{
    unordered_multiset<uint> fix;
    unordered_multiset<uint> variable;
  public:
    compact_multiset(const list<uint>& _fix, const list<uint>& _variable){
      multiply_with_set(_fix);
      multiply_with_subsets(_variable);
    }
    // constructor
    compact_multiset():fix(), variable() {};
    // clear
    inline void clear(){
      fix.clear();
      variable.clear();
    }
    // multiply all sets in the list with one set X (add this set to all sets in the list)
    inline void multiply_with_set(const list<uint>& X){
      for(uint x: X) fix.insert(x);
    }
    // multiply all sets in the list with all subsets of X (add all subsets to all sets in the list)
    inline void multiply_with_subsets(const list<uint>& X){
      for(uint x: X) variable.insert(x);
    }
    // return whether there are some elements that are fix
    inline bool has_no_fixed_elements() const {
      return fix.empty();
    }

    // add or remove a (list of) items from all sets
    inline void add_to_all(const uint x){
      fix.insert(x);
    }
    inline void remove_all_from_all(const uint x){
      fix.erase(x);
      variable.erase(x);
    }
    inline void remove_once_from_all(const uint x){
      erase_once(fix, x);
      erase_once(variable, x);
    }
    inline void add_to_all(const list<uint>& X){
      for(auto x: X) add_to_all(x);
    }
    inline void remove_all_from_all(const list<uint>& X){
      for(auto x: X) remove_all_from_all(x);
    }

    // check containment of (lists of) items in all/some of the sets
    bool all_sets_contain(const uint x) const{
      return fix.find(x) != fix.end();
    }
    bool some_sets_contain(const uint x) const{
      return (fix.find(x) != fix.end()) || (variable.find(x) != variable.end());
    }
    bool all_sets_contain_all_of(const list<uint>& X) const{
      for(auto x: X) if(!(all_sets_contain(x))) return false;
      return true;
    }
    bool some_sets_contain_all_of(const list<uint>& X) const{
      // due to the nature of this list of sets, this is equal to: for each x in X, there is some set containing x
      for(auto x: X) if(!(some_sets_contain(x))) return false;
      return true;
    }

    // remove all sets that (do not) contain an item
    void remove_sets_containing(const uint x){
      if(fix.find(x) == fix.end()) // if x not in fix, then just remove all occurances of it from variable
        variable.erase(x);
      else clear(); // if x in fix, then all sets contain x, so clear the compact_set
    }
    void remove_sets_not_containing(const uint x){
      if(fix.find(x) == fix.end()){ // if x in fix, then x is already in all sets
        if(variable.find(x) != variable.end()){ // if x is optional, then make it madatory
          fix.insert(x);
          variable.erase(x);
        } else clear(); // if x is not in any set, then remove all sets
      }
    }
    // remove all sets that (do not) contain an item
    void remove_sets_containing_any_of(const list<uint>& X){
      for(auto x: X) remove_sets_containing(x);
    }
    void remove_sets_not_containing_all_of(const list<uint>& X){
      for(auto x: X) remove_sets_not_containing(x);
    }

    // return whether the list contains no element
    bool empty() const{
      return fix.empty() && variable.empty();
    }
    friend ostream& ::operator<<(ostream& os, const status::compact_multiset& s);
  };
}





#endif
