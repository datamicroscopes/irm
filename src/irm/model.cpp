#include <microscopes/irm/model.hpp>
#include <microscopes/common/util.hpp>
#include <distributions/special.hpp>

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

void
state::random_initialize(const dataset_t &d, rng_t &rng)
{
  // XXX: generate initial assignments

  MICROSCOPES_DCHECK(is_correct_shape(d), "not presented with the relations");
  for (size_t i = 0; i < relations_.size(); i++) {
    auto &relation = relations_[i];
    MICROSCOPES_DCHECK(
        relation.suffstats_table_.empty() &&
        relation.ident_table_.empty() &&
        !relation.ident_gen_, "data already present");
    for (const auto &p : *d[i]) {
      vector<size_t> groups;
      groups.reserve(relation.desc_.domains_.size());
      for (size_t i = 0; i < relation.desc_.domains_.size(); i++) {
        MICROSCOPES_DCHECK(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]] != -1, "unassigned");
        groups.push_back(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]]);
      }
      shared_ptr<feature_group> group;
      auto it = relation.suffstats_table_.find(groups);
      if (it == relation.suffstats_table_.end()) {
        group = relation.desc_.model_->create_feature_group(rng);
        auto &ss = relation.suffstats_table_[groups];
        ss.ident_ = relation.ident_gen_++;
        ss.count_ = 1;
        ss.ss_ = group;
        MICROSCOPES_ASSERT(relation.ident_table_.find(ss.ident_) == relation.ident_table_.end());
        relation.ident_table_[ss.ident_] = groups;
      } else {
        MICROSCOPES_ASSERT(it->second.count_);
        it->second.count_++;
        group = it->second.ss_;
      }
      MICROSCOPES_ASSERT(group);
      group->add_value(*relation.desc_.model_, p.second, rng);
    }
  }
}

// XXX: all the *_value methods have too much code duplication

void
state::add_value(size_t domain, size_t gid, size_t eid, const dataset_t &d, rng_t &rng)
{
  MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
  MICROSCOPES_DCHECK(domains_[domain].isactivegroup(gid), "invalid gid");
  MICROSCOPES_DCHECK(eid < domains_[domain].nentities(), "invalid eid");
  MICROSCOPES_DCHECK(is_correct_shape(d), "not presented with the relations");

  for (const auto &dr : domain_relations_[domain]) {
    auto &relation = relations_[dr.rel_];
    auto &data = d[dr.rel_];

    vector<size_t> ignore_idxs;
    for (size_t i = 0; i < dr.pos_; i++)
      if (relation.desc_.domains_[i] == domain)
        ignore_idxs.push_back(i);

    vector<size_t> groups;
    for (const auto &p : data->slice(dr.pos_, eid)) {
      // don't double count
      bool skip = false;
      for (auto idx : ignore_idxs)
        if (p.first[idx] == eid) {
          skip = true;
          break;
        }
      if (skip)
        continue;

      groups.clear();
      for (size_t i = 0; i < relation.desc_.domains_.size(); i++) {
        if (relation.desc_.domains_[i] == domain && p.first[i] == eid)
          groups.push_back(gid);
        else {
          MICROSCOPES_DCHECK(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]] != -1, "unassigned");
          groups.push_back(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]]);
        }
      }
      shared_ptr<feature_group> group;
      auto it = relation.suffstats_table_.find(groups);
      if (it == relation.suffstats_table_.end()) {
        group = relation.desc_.model_->create_feature_group(rng);
        auto &ss = relation.suffstats_table_[groups];
        ss.ident_ = relation.ident_gen_++;
        ss.count_ = 1;
        ss.ss_ = group;
        MICROSCOPES_ASSERT(relation.ident_table_.find(ss.ident_) == relation.ident_table_.end());
        relation.ident_table_[ss.ident_] = groups;
      } else {
        MICROSCOPES_ASSERT(it->second.count_);
        it->second.count_++;
        group = it->second.ss_;
      }
      MICROSCOPES_ASSERT(group);
      group->add_value(*relation.desc_.model_, p.second, rng);
    }
  }

  domains_[domain].add_value(gid, eid);
}

size_t
state::remove_value(size_t domain, size_t eid, const dataset_t &d, rng_t &rng)
{
  MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
  MICROSCOPES_DCHECK(eid < domains_[domain].nentities(), "invalid eid");
  MICROSCOPES_DCHECK(is_correct_shape(d), "not presented with the relations");

  for (const auto &dr : domain_relations_[domain]) {
    auto &relation = relations_[dr.rel_];
    auto &data = d[dr.rel_];

#ifdef DEBUG_MODE
    vector<size_t> expected_shape;
    for (auto d : relation.desc_.domains_)
      expected_shape.push_back(domains_[d].nentities());
    MICROSCOPES_DCHECK(expected_shape == data->shape(), "improper relation shape given");
#endif

    vector<size_t> ignore_idxs;
    for (size_t i = 0; i < dr.pos_; i++)
      if (relation.desc_.domains_[i] == domain)
        ignore_idxs.push_back(i);

    vector<size_t> groups;
    for (const auto &p : data->slice(dr.pos_, eid)) {
      // don't double count
      bool skip = false;
      for (auto idx : ignore_idxs)
        if (p.first[idx] == eid) {
          skip = true;
          break;
        }
      if (skip)
        continue;

      groups.clear();
      for (size_t i = 0; i < relation.desc_.domains_.size(); i++) {
        MICROSCOPES_DCHECK(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]] != -1, "unassigned");
        groups.push_back(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]]);
      }
      auto it = relation.suffstats_table_.find(groups);
      MICROSCOPES_ASSERT(it != relation.suffstats_table_.end());
      MICROSCOPES_ASSERT(it->second.count_);
      MICROSCOPES_ASSERT(it->second.ss_);
      MICROSCOPES_ASSERT(
          relation.ident_table_.find(it->second.ident_) != relation.ident_table_.end() &&
          relation.ident_table_[it->second.ident_] == groups);

      it->second.ss_->remove_value(*relation.desc_.model_, p.second, rng);
      if (!--it->second.count_) {
        relation.ident_table_.erase(it->second.ident_);
        relation.suffstats_table_.erase(it);
      }
    }
  }

  return domains_[domain].remove_value(eid);
}

pair<vector<size_t>, vector<float>>
state::score_value(size_t did, size_t eid, const dataset_t &d, rng_t &rng) const
{
  MICROSCOPES_DCHECK(is_correct_shape(d), "not presented with the relations");
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
    for (const auto &dr : domain_relations_[did]) {
      auto &relation = relations_[dr.rel_];
      auto &data = d[dr.rel_];

      vector<size_t> ignore_idxs;
      for (size_t i = 0; i < dr.pos_; i++)
        if (relation.desc_.domains_[i] == did)
          ignore_idxs.push_back(i);

      vector<size_t> groups;
      for (const auto &p : data->slice(dr.pos_, eid)) {
        // don't double count
        bool skip = false;
        for (auto idx : ignore_idxs)
          if (p.first[idx] == eid) {
            skip = true;
            break;
          }
        if (skip)
          continue;

        groups.clear();
        for (size_t i = 0; i < relation.desc_.domains_.size(); i++) {
          if (relation.desc_.domains_[i] == did && p.first[i] == eid)
            groups.push_back(g.first);
          else {
            MICROSCOPES_DCHECK(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]] != -1, "unassigned");
            groups.push_back(domains_[relation.desc_.domains_[i]].assignments()[p.first[i]]);
          }
        }
        shared_ptr<feature_group> group;
        auto it = relation.suffstats_table_.find(groups);
        if (it == relation.suffstats_table_.end()) {
          group = relation.desc_.model_->create_feature_group(rng);
        } else {
          MICROSCOPES_ASSERT(it->second.count_);
          group = it->second.ss_;
        }
        MICROSCOPES_ASSERT(group);
        sum += group->score_value(*relation.desc_.model_, p.second, rng);
      }
    }
    ret.first.push_back(g.first);
    ret.second.push_back(sum);
    count += g.second;
  }
  const float lgnorm = fast_log(float(count) + domain.alpha_);
  for (auto &s : ret.second)
    s -= lgnorm;
  return ret;
}
