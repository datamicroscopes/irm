#include <microscopes/irm/model.hpp>
#include <distributions/special.hpp>

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
    relations_.emplace_back(r, map<tuple_t, group_with_count_t>());
  }
  domain_relations_.reserve(domains_.size());
  for (size_t i = 0; i < domains_.size(); i++)
    domain_relations_.emplace_back(domain_relations(i));
}

void
state::add_value(size_t domain, size_t gid, size_t eid, const dataset_t &d, rng_t &rng)
{
  MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
  MICROSCOPES_DCHECK(domains_[domain].isactivegroup(gid), "invalid gid");
  MICROSCOPES_DCHECK(eid < domains_[domain].nentities(), "invalid eid");
  MICROSCOPES_DCHECK(d.size() == relations_.size(), "not presented with the relations");

  for (const auto &dr : domain_relations_[domain]) {
    auto &relation = relations_[dr.first];
    auto &data = d[dr.first];

#ifdef DEBUG_MODE
    vector<size_t> expected_shape;
    for (auto d : relation.first.domains_)
      expected_shape.push_back(domains_[d].nentities());
    MICROSCOPES_DCHECK(expected_shape == data->shape(), "improper relation shape given");
#endif

    vector<size_t> groups;
    for (const auto &p : data->slice(dr.second, eid)) {
      groups.clear();
      for (size_t i = 0; i < relation.first.domains_.size(); i++) {
        if (relation.first.domains_[i] == domain && p.first[i] == eid)
          groups.push_back(gid);
        else {
          MICROSCOPES_DCHECK(domains_[relation.first.domains_[i]].assignments()[p.first[i]] != -1, "unassigned");
          groups.push_back(domains_[relation.first.domains_[i]].assignments()[p.first[i]]);
        }
      }
      shared_ptr<feature_group> group;
      auto it = relation.second.find(groups);
      if (it == relation.second.end()) {
        group = relation.first.model_->create_feature_group(rng);
        auto &g = relation.second[groups];
        g.count_ = 1;
        g.group_ = group;
      } else {
        MICROSCOPES_ASSERT(it->second.count_);
        it->second.count_++;
        group = it->second.group_;
      }
      MICROSCOPES_ASSERT(group);
      group->add_value(*relation.first.model_, p.second, rng);
    }
  }

  domains_[domain].add_value(gid, eid);
}

size_t
state::remove_value(size_t domain, size_t eid, const dataset_t &d, rng_t &rng)
{
  MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
  MICROSCOPES_DCHECK(eid < domains_[domain].nentities(), "invalid eid");
  MICROSCOPES_DCHECK(d.size() == relations_.size(), "not presented with the relations");

  for (const auto &dr : domain_relations_[domain]) {
    auto &relation = relations_[dr.first];
    auto &data = d[dr.first];

#ifdef DEBUG_MODE
    vector<size_t> expected_shape;
    for (auto d : relation.first.domains_)
      expected_shape.push_back(domains_[d].nentities());
    MICROSCOPES_DCHECK(expected_shape == data->shape(), "improper relation shape given");
#endif

    vector<size_t> groups;
    for (const auto &p : data->slice(dr.second, eid)) {
      groups.clear();
      for (size_t i = 0; i < relation.first.domains_.size(); i++) {
        MICROSCOPES_DCHECK(domains_[relation.first.domains_[i]].assignments()[p.first[i]] != -1, "unassigned");
        groups.push_back(domains_[relation.first.domains_[i]].assignments()[p.first[i]]);
      }
      auto it = relation.second.find(groups);
      MICROSCOPES_ASSERT(it != relation.second.end());
      MICROSCOPES_ASSERT(it->second.count_);
      MICROSCOPES_ASSERT(it->second.group_);

      it->second.group_->remove_value(*relation.first.model_, p.second, rng);
      if (!--it->second.count_)
        relation.second.erase(it);
    }
  }

  return domains_[domain].remove_value(eid);
}

pair<vector<size_t>, vector<float>>
state::score_value(size_t did, size_t eid, const dataset_t &d, rng_t &rng) const
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
    for (const auto &dr : domain_relations_[did]) {
      auto &relation = relations_[dr.first];
      auto &data = d[dr.first];

#ifdef DEBUG_MODE
      vector<size_t> expected_shape;
      for (auto d : relation.first.domains_)
        expected_shape.push_back(domains_[d].nentities());
      MICROSCOPES_DCHECK(expected_shape == data->shape(), "improper relation shape given");
#endif

      vector<size_t> groups;
      for (const auto &p : data->slice(dr.second, eid)) {
        groups.clear();
        for (size_t i = 0; i < relation.first.domains_.size(); i++) {
          if (relation.first.domains_[i] == did && p.first[i] == eid)
            groups.push_back(g.first);
          else {
            MICROSCOPES_DCHECK(domains_[relation.first.domains_[i]].assignments()[p.first[i]] != -1, "unassigned");
            groups.push_back(domains_[relation.first.domains_[i]].assignments()[p.first[i]]);
          }
        }
        shared_ptr<feature_group> group;
        auto it = relation.second.find(groups);
        if (it == relation.second.end()) {
          group = relation.first.model_->create_feature_group(rng);
        } else {
          MICROSCOPES_ASSERT(it->second.count_);
          group = it->second.group_;
        }
        MICROSCOPES_ASSERT(group);
        sum += group->score_value(*relation.first.model_, p.second, rng);
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
