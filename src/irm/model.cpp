#include <microscopes/irm/model.hpp>

namespace microscopes {
namespace irm {

template class state<-1>;
template class state<2>;
template class state<3>;
template class state<4>;

template class model<-1>;
template class model<2>;
template class model<3>;
template class model<4>;

} // namespace irm
} // namespace microscopes

using namespace std;
using namespace microscopes::common;
using namespace microscopes::common::relation;
using namespace microscopes::irm;
using namespace microscopes::irm::detail;
using namespace microscopes::io;
using namespace microscopes::models;

model_definition::model_definition(
    const vector<size_t> &domains,
    const vector<relation_definition> &relations)
  : domains_(domains), relations_(relations)
{
  MICROSCOPES_DCHECK(domains.size(), "no domains given");
  MICROSCOPES_DCHECK(relations.size(), "no relations given");
#ifdef DEBUG_BUILD
  for (auto s : domains)
    MICROSCOPES_DCHECK(s, "empty domain given");
  for (const auto &r : relations) {
    for (auto d : r.domains())
      MICROSCOPES_DCHECK(d < domains.size(), "invalid domain given");
  }
#endif
}
