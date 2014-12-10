#include "compact_set.hpp"

namespace status {


}

ostream& operator<<(ostream& os, const status::compact_multiset& s){
  // the set is somewhat well described by the fix set times the powerset of variable
  return os << s.fix << "xP" << s.variable;
}

