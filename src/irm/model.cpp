#include <microscopes/irm/model.hpp>
#include <microscopes/common/util.hpp>
#include <distributions/special.hpp>

#include <algorithm>

using namespace std;
using namespace distributions;
using namespace microscopes::common;
using namespace microscopes::common::sparse_ndarray;
using namespace microscopes::irm;
using namespace microscopes::irm::detail;
using namespace microscopes::models;

state::state(const vector<size_t> &domains,
             const vector<relation_t> &relations)
  : domains_(), relations_()
{
  MICROSCOPES_DCHECK(domains.size(), "no domains given");
  domains_.reserve(domains.size());
  for (auto n : domains)
    domains_.emplace_back(n);
  for (auto s : domains)
    MICROSCOPES_DCHECK(s, "empty domain given");
  for (const auto &r : relations) {
    for (auto d : r.domains_)
      MICROSCOPES_DCHECK(d < domains.size(), "invalid domain given");
    relations_.emplace_back(r);
  }
  domain_relations_.reserve(domains_.size());
  for (size_t i = 0; i < domains_.size(); i++)
    domain_relations_.emplace_back(domain_relations(i));
}

static inline void
AssertAllEntitiesAccounted(const vector<set<size_t>> &c, size_t n)
{
  vector<bool> ents(n, false);
  for (const auto &s : c)
    for (auto eid : s) {
      MICROSCOPES_DCHECK(eid < ents.size(), "bad eid given");
      MICROSCOPES_DCHECK(!ents[eid], "eid given twice");
      ents[eid] = true;
    }
  for (auto b : ents)
    MICROSCOPES_DCHECK(b, "ent unaccounted for");
}

void
state::eids_to_gids_under_relation(
    vector<size_t> &gids,
    const vector<size_t> &eids,
    const relation_t &desc) const
{
  gids.clear();
  gids.reserve(desc.domains_.size());
  for (size_t i = 0; i < desc.domains_.size(); i++) {
    MICROSCOPES_DCHECK(
        domains_[desc.domains_[i]].assignments()[eids[i]] != -1,
        "eid is not assigned to a valid group");
    gids.push_back(domains_[desc.domains_[i]].assignments()[eids[i]]);
  }
}

void
state::add_value_to_feature_group(
    const vector<size_t> &gids,
    const value_accessor &value,
    relation_container_t &relation,
    rng_t &rng,
    float *acc_score)
{
  MICROSCOPES_ASSERT(!value.anymasked());
  shared_ptr<feature_group> group;
  auto it = relation.suffstats_table_.find(gids);
  if (it == relation.suffstats_table_.end()) {
    group = relation.desc_.model_->create_feature_group(rng);
    auto &ss = relation.suffstats_table_[gids];
    ss.ident_ = relation.ident_gen_++;
    ss.count_ = 1;
    MICROSCOPES_ASSERT(!ss.ss_);
    ss.ss_ = group;
    MICROSCOPES_ASSERT(relation.ident_table_.find(ss.ident_) == relation.ident_table_.end());
    relation.ident_table_[ss.ident_] = gids;
  } else {
    it->second.count_++;
    group = it->second.ss_;
  }
  MICROSCOPES_ASSERT(group);
  if (acc_score)
    *acc_score += group->score_value(*relation.desc_.model_, value, rng);
  group->add_value(*relation.desc_.model_, value, rng);
}

void
state::remove_value_from_feature_group(
    const vector<size_t> &gids,
    const value_accessor &value,
    relation_container_t &relation,
    rng_t &rng)
{
  auto it = relation.suffstats_table_.find(gids);
  MICROSCOPES_ASSERT(!value.anymasked());
  MICROSCOPES_ASSERT(it != relation.suffstats_table_.end());
  MICROSCOPES_ASSERT(it->second.count_);
  MICROSCOPES_ASSERT(it->second.ss_);
  MICROSCOPES_ASSERT(
      relation.ident_table_.find(it->second.ident_) != relation.ident_table_.end() &&
      relation.ident_table_[it->second.ident_] == gids);
  it->second.ss_->remove_value(*relation.desc_.model_, value, rng);
  it->second.count_--;
  // XXX: unfortunately, we cannot clean this up now!! this is because for
  // non-conjugate models, score_value() depends on the randomness we sampled
  // in add_value_to_feature_group() for correctness. yes, this is quite hacky
  // (for conjugate models, it is safe to delete here)
  //
  // note that, the point at which the suffstat can be GC-ed for non-conj
  // models is when >= of the gids associated with it is no longer a valid gid
  // (which should imply the count is zero also)
  //if (!--it->second.count_) {
  //  relation.ident_table_.erase(it->second.ident_);
  //  relation.suffstats_table_.erase(it);
  //}
}

void
state::iterate_over_entity_data(
    size_t domain,
    size_t eid,
    const dataset_t &d,
    function<void(size_t, const vector<size_t> &, const value_accessor &)> callback) const
{
  for (const auto &dr : domain_relations_[domain]) {
    auto &relation = relations_[dr.rel_];
    auto &data = d[dr.rel_];
    vector<size_t> ignore_idxs;
    for (size_t i = 0; i < dr.pos_; i++)
      if (relation.desc_.domains_[i] == domain)
        ignore_idxs.push_back(i);
    for (const auto &p : data->slice(dr.pos_, eid)) {
      // don't double count
      bool skip = false;
      for (auto idx : ignore_idxs) {
        if (p.first[idx] == eid) {
          skip = true;
          break;
        }
      }
      if (skip)
        continue;
      callback(dr.rel_, p.first, p.second);
    }
  }
}

vector< vector<size_t> >
state::entity_data_positions(size_t domain, size_t eid, const dataset_t &d) const
{
  vector< vector<size_t> > ret;
  iterate_over_entity_data(
      domain, eid, d,
      [&ret](size_t, const vector<size_t> &eids, const value_accessor &) {
        ret.emplace_back(eids);
      });
  return ret;
}

void
state::initialize(const vector<vector<set<size_t>>> &clusters, const dataset_t &d, rng_t &rng)
{
  MICROSCOPES_DCHECK(clusters.size() == domains_.size(), "invalid number of clusterings");
  assert_correct_shape(d);

  for (size_t i = 0; i < domains_.size(); i++) {
    MICROSCOPES_DCHECK(!domains_[i].ngroups(), "domain not empty");
    const auto &c = clusters[i];
    MICROSCOPES_DCHECK(c.size() > 0, "no clusters given!");
    for (size_t j = 0; j < c.size(); j++)
      domains_[i].create_group();
    AssertAllEntitiesAccounted(c, domains_[i].nentities());
    for (size_t j = 0; j < c.size(); j++)
      for (auto eid : c[j])
        domains_[i].add_value(j, eid);
  }
#ifdef DEBUG_MODE
  for (auto &d : domains_)
    for (auto s : d.assignments())
      MICROSCOPES_DCHECK(s != -1, "assignments should all be filled");
#endif

  vector<size_t> gids;
  for (size_t i = 0; i < relations_.size(); i++) {
    auto &relation = relations_[i];
    for (const auto &p : *d[i]) {
      eids_to_gids_under_relation(gids, p.first, relation.desc_);
      add_value_to_feature_group(gids, p.second, relation, rng, nullptr);
    }
  }
}

void
state::random_initialize(const dataset_t &d, rng_t &rng)
{
  assert_correct_shape(d);
  vector<vector<set<size_t>>> clusters;
  clusters.reserve(domains_.size());
  for (auto &d : domains_) {
    vector<set<size_t>> cluster;
    // create min(100, n/2) + 1 groups
    const size_t ngroups = min(size_t(100), d.nentities()) + 1;
    cluster.resize(ngroups);
    const auto groups = util::range(ngroups);
    for (size_t i = 0; i < d.nentities(); i++) {
      const auto choice = util::sample_choice(groups, rng);
      cluster[choice].insert(i);
    }
    clusters.emplace_back(cluster);
  }
  initialize(clusters, d, rng);
}

void
state::add_value0(size_t domain, size_t gid, size_t eid, const dataset_t &d, rng_t &rng, float *acc_score)
{
  domains_[domain].add_value(gid, eid);
  vector<size_t> gids;
  iterate_over_entity_data(
      domain, eid, d,
      [this, &gids, &rng, acc_score](
         size_t rid,
         const vector<size_t> &eids,
         const value_accessor &value) {
          auto &relation = this->relations_[rid];
          this->eids_to_gids_under_relation(gids, eids, relation.desc_);
          this->add_value_to_feature_group(gids, value, relation, rng, acc_score);
      });
}

size_t
state::remove_value0(size_t domain, size_t eid, const dataset_t &d, rng_t &rng)
{
  vector<size_t> gids;
  iterate_over_entity_data(
      domain, eid, d,
      [this, &gids, &rng](
         size_t rid,
         const vector<size_t> &eids,
         const value_accessor &value) {
          auto &relation = this->relations_[rid];
          this->eids_to_gids_under_relation(gids, eids, relation.desc_);
          this->remove_value_from_feature_group(gids, value, relation, rng);
      });
  return domains_[domain].remove_value(eid).first;
}

pair<vector<size_t>, vector<float>>
state::score_value0(size_t did, size_t eid, const dataset_t &d, rng_t &rng) const
{
  const auto &domain = domains_[did];
  MICROSCOPES_DCHECK(!domain.empty_groups().empty(), "no empty groups");

  pair<vector<size_t>, vector<float>> ret;
  ret.first.reserve(domain.ngroups());
  ret.second.reserve(domain.ngroups());

  float pseudocounts = 0;
  for (const auto &g : domain) {
    const float pseudocount = domain.pseudocount(g.first, g.second);
    float sum = fast_log(pseudocount);
    const_cast<state *>(this)->add_value0(did, g.first, eid, d, rng, &sum);
    const size_t gid = const_cast<state *>(this)->remove_value0(did, eid, d, rng);
    if (unlikely(gid != g.first))
      MICROSCOPES_ASSERT(false);
    ret.first.push_back(g.first);
    ret.second.push_back(sum);
    pseudocounts += pseudocount;
  }
  const float lgnorm = fast_log(pseudocounts);
  for (auto &s : ret.second)
    s -= lgnorm;
  return ret;
}

float
state::score_likelihood(size_t relation, ident_t id, rng_t &rng) const
{
  auto &ss = get_suffstats_t(relation, id);
  return ss.ss_->score_data(*relations_[relation].desc_.model_, rng);
}

float
state::score_likelihood(size_t relation, rng_t &rng) const
{
  MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
  float score = 0.;
  auto &m = relations_[relation].desc_.model_;
  for (auto &p : relations_[relation].suffstats_table_)
    score += p.second.ss_->score_data(*m, rng);
  return score;
}
