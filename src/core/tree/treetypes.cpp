#include "treetypes.h"

namespace phylogeny {

std::ostream& operator<< (std::ostream &os, GID gid) {
  return os << std::underlying_type<GID>::type(gid);
}

std::ostream& operator<< (std::ostream &os, SID sid) {
  return os << std::underlying_type<SID>::type(sid);
}

} // end of namespace phylogeny
