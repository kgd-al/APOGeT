#include "treetypes.h"

namespace phylogeny {

std::ostream& operator<< (std::ostream &os, SID sid) {
  return os << std::underlying_type<SID>::type(sid);
}

} // end of namespace phylogeny
