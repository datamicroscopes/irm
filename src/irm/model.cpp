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

// XXX: don't copy
float
domain::score_assignment() const
{
  map<size_t, size_t> counts;
  MICROSCOPES_ASSERT(assignments_[0] != -1);
  counts[assignments_[0]] = 1;
  float sum = 0.;
  for (size_t i = 1; i < assignments_.size(); i++) {
    const ssize_t gid = assignments_[i];
    MICROSCOPES_ASSERT(gid != -1);
    const auto it = counts.find(gid);
    const bool found = (it != counts.end());
    const float numer = (!found) ? alpha_ : it->second;
    const float denom = float(i) + alpha_;
    sum += fast_log(numer / denom);
    if (found)
      it->second++;
    else
      counts[gid] = 1;
  }
  return sum;
}

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
    MICROSCOPES_ASSERT(domains_[desc.domains_[i]].assignments()[eids[i]] != -1);
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
  shared_ptr<feature_group> group;
  auto it = relation.suffstats_table_.find(gids);
  if (it == relation.suffstats_table_.end()) {
    group = relation.desc_.model_->create_feature_group(rng);
    auto &ss = relation.suffstats_table_[gids];
    ss.ident_ = relation.ident_gen_++;
    ss.count_ = 1;
    ss.ss_ = group;
    MICROSCOPES_ASSERT(relation.ident_table_.find(ss.ident_) == relation.ident_table_.end());
    relation.ident_table_[ss.ident_] = gids;
  } else {
    MICROSCOPES_ASSERT(it->second.count_);
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
  MICROSCOPES_ASSERT(it != relation.suffstats_table_.end());
  MICROSCOPES_ASSERT(it->second.count_);
  MICROSCOPES_ASSERT(it->second.ss_);
  MICROSCOPES_ASSERT(
      relation.ident_table_.find(it->second.ident_) != relation.ident_table_.end() &&
      relation.ident_table_[it->second.ident_] == gids);
  it->second.ss_->remove_value(*relation.desc_.model_, value, rng);
  if (!--it->second.count_) {
    relation.ident_table_.erase(it->second.ident_);
    relation.suffstats_table_.erase(it);
  }
}

void
state::initialize(const vector<vector<set<size_t>>> &clusters, const dataset_t &d, rng_t &rng)
{
  MICROSCOPES_DCHECK(clusters.size() == domains_.size(), "invalid number of clusterings");
  assert_correct_shape(d);

  reset();

  for (size_t i = 0; i < domains_.size(); i++) {
    const auto &c = clusters[i];
    MICROSCOPES_DCHECK(c.size() > 0, "no clusters given!");
    for (size_t j = 0; j < c.size(); j++)
      domains_[i].create_group();
    AssertAllEntitiesAccounted(c, domains_[i].nentities());
    for (size_t j = 0; j < c.size(); j++)
      for (auto eid : c[j])
        domains_[i].add_value(j, eid);
  }

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
  for (const auto &dr : domain_relations_[domain]) {
    auto &relation = relations_[dr.rel_];
    auto &data = d[dr.rel_];
    vector<size_t> ignore_idxs;
    for (size_t i = 0; i < dr.pos_; i++)
      if (relation.desc_.domains_[i] == domain)
        ignore_idxs.push_back(i);
    vector<size_t> gids;
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
      eids_to_gids_under_relation(gids, p.first, relation.desc_);
      add_value_to_feature_group(gids, p.second, relation, rng, acc_score);
    }
  }
}

size_t
state::remove_value0(size_t domain, size_t eid, const dataset_t &d, rng_t &rng)
{
  // XXX: we duplicate this iteration code twice (also in add_value0)
  for (const auto &dr : domain_relations_[domain]) {
    auto &relation = relations_[dr.rel_];
    auto &data = d[dr.rel_];
    vector<size_t> ignore_idxs;
    for (size_t i = 0; i < dr.pos_; i++)
      if (relation.desc_.domains_[i] == domain)
        ignore_idxs.push_back(i);
    vector<size_t> gids;
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
      eids_to_gids_under_relation(gids, p.first, relation.desc_);
      remove_value_from_feature_group(gids, p.second, relation, rng);
    }
  }
  return domains_[domain].remove_value(eid);
}

pair<vector<size_t>, vector<float>>
state::score_value0(size_t did, size_t eid, const dataset_t &d, rng_t &rng) const
{
  const auto &domain = domains_[did];
  MICROSCOPES_DCHECK(!domain.empty_groups().empty(), "no empty groups");
  const auto &groups = domain.groups_;

  pair<vector<size_t>, vector<float>> ret;
  ret.first.reserve(groups.size());
  ret.second.reserve(groups.size());

  const float empty_group_alpha = domain.alpha_ / float(domain.empty_groups().size());
  size_t count = 0;
  for (const auto &g : groups) {
    float sum = g.second ? float(g.second) : empty_group_alpha;
    const_cast<state *>(this)->add_value0(did, g.first, eid, d, rng, &sum);
    const size_t gid = const_cast<state *>(this)->remove_value0(did, eid, d, rng);
    if (unlikely(gid != g.first))
      MICROSCOPES_ASSERT(false);
    ret.first.push_back(g.first);
    ret.second.push_back(sum);
    count += g.second;
  }
  const float lgnorm = fast_log(float(count) + domain.alpha_);
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

void
state::reset()
{
  for (auto &d : domains_)
    d.reset();
  for (auto &r : relations_) {
    r.suffstats_table_.clear();
    r.ident_table_.clear();
    r.ident_gen_ = 0;
  }
}
