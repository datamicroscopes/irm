#pragma once

#include <microscopes/common/sparse_ndarray/dataview.hpp>
#include <microscopes/common/typedefs.hpp>
#include <microscopes/common/macros.hpp>
#include <microscopes/common/assert.hpp>
#include <microscopes/models/base.hpp>
#include <microscopes/io/schema.pb.h>

#include <cmath>
#include <vector>
#include <set>
#include <functional>
#include <map>
#include <memory>
#include <sstream>
#include <utility>
#include <stdexcept>

namespace microscopes {
namespace irm {

/**
 * XXX: we are duplicating the group management code found in
 * microscopes::mixture::model::state. abstract that away
 */
class domain {
public:
  typedef io::CRP message_type;

  domain(size_t n)
    : alpha_(),
      gcount_(),
      gempty_(),
      assignments_(n, -1),
      groups_()
  {}

  inline common::hyperparam_bag_t
  get_hp() const
  {
    message_type m;
    m.set_alpha(alpha_);
    std::ostringstream out;
    m.SerializeToOstream(&out);
    return out.str();
  }

  inline void
  set_hp(const common::hyperparam_bag_t &hp)
  {
    std::istringstream inp(hp);
    message_type m;
    m.ParseFromIstream(&inp);
    alpha_ = m.alpha();
  }

  inline void *
  get_hp_raw_ptr(const std::string &key)
  {
    if (key == "alpha")
      return &alpha_;
    throw std::runtime_error("unknown key: " + key);
  }

  inline const std::vector<ssize_t> &
  assignments() const
  {
    return assignments_;
  }

  inline const std::set<size_t> &
  empty_groups() const
  {
    return gempty_;
  }

  inline size_t nentities() const { return assignments_.size(); }
  inline size_t ngroups() const { return groups_.size(); }

  inline bool
  isactivegroup(size_t gid) const
  {
    return groups_.find(gid) != groups_.end();
  }

  inline size_t
  groupsize(size_t gid) const
  {
    const auto it = groups_.find(gid);
    MICROSCOPES_ASSERT(it != groups_.end());
    return it->second;
  }

  inline std::vector<size_t>
  groups() const
  {
    std::vector<size_t> ret;
    ret.reserve(ngroups());
    for (auto &g : groups_)
      ret.push_back(g.first);
    return ret;
  }

  inline size_t
  create_group()
  {
    const size_t gid = gcount_++;
    groups_[gid] = 0;
    MICROSCOPES_ASSERT(!gempty_.count(gid));
    gempty_.insert(gid);
    return gid;
  }

  inline void
  delete_group(size_t gid)
  {
    auto it = groups_.find(gid);
    MICROSCOPES_ASSERT(it != groups_.end());
    MICROSCOPES_ASSERT(!it->second);
    MICROSCOPES_ASSERT(gempty_.count(gid));
    groups_.erase(it);
    gempty_.erase(gid);
  }

  inline void
  add_value(size_t gid, size_t eid)
  {
    MICROSCOPES_ASSERT(assignments_.at(eid) == -1);
    auto it = groups_.find(gid);
    MICROSCOPES_ASSERT(it != groups_.end());
    if (!it->second++) {
      MICROSCOPES_ASSERT(gempty_.count(gid));
      gempty_.erase(gid);
      MICROSCOPES_ASSERT(!gempty_.count(gid));
    } else {
      MICROSCOPES_ASSERT(!gempty_.count(gid));
    }
    assignments_[eid] = gid;
  }

  inline size_t
  remove_value(size_t eid)
  {
    MICROSCOPES_ASSERT(assignments_.at(eid) != -1);
    const size_t gid = assignments_[eid];
    auto it = groups_.find(gid);
    MICROSCOPES_ASSERT(it != groups_.end());
    MICROSCOPES_ASSERT(!gempty_.count(gid));
    if (!--it->second)
      gempty_.insert(gid);
    assignments_[eid] = -1;
    return gid;
  }

private:
  float alpha_;
  size_t gcount_;
  std::set<size_t> gempty_;
  std::vector<ssize_t> assignments_;
  std::map<size_t, size_t> groups_;
};

typedef std::vector<std::shared_ptr<common::sparse_ndarray::dataview>> dataset_t;

class state {
public:

  typedef domain::message_type message_type;

  typedef std::vector<size_t> tuple_t;

  struct group_with_count_t {
    group_with_count_t() : count_(), group_() {}
    unsigned count_; // a ref count, so we know when to remove
    std::shared_ptr<models::feature_group> group_;
  };

  struct relation_t {
    relation_t() = default;
    relation_t(const std::vector<size_t> &domains,
               const std::shared_ptr<models::model> &model)
      : domains_(domains), model_(model) {}
    std::vector<size_t> domains_;
    std::shared_ptr<models::model> model_;
  };

  state(const std::vector<size_t> &domains,
        const std::vector<relation_t> &relations);

  inline common::hyperparam_bag_t
  get_domain_hp(size_t domain) const
  {
    return domains_[domain].get_hp();
  }

  inline void
  set_domain_hp(size_t domain, const common::hyperparam_bag_t &hp)
  {
    domains_[domain].set_hp(hp);
  }

  inline common::hyperparam_bag_t
  get_relation_hp(size_t relation) const
  {
    return relations_[relation].first.model_->get_hp();
  }

  inline void
  set_relation_hp(size_t relation, const common::hyperparam_bag_t &hp)
  {
    relations_[relation].first.model_->set_hp(hp);
  }

  inline size_t
  create_group(size_t domain)
  {
    return domains_[domain].create_group();
  }

  inline void
  delete_group(size_t domain, size_t gid)
  {
    domains_[domain].delete_group(gid);
  }

  // assigns a value to a group w/o associating it with any particular piece of
  // data; should only be invoked during bootstrapping phases
  inline void
  assign_value(size_t domain, size_t gid, size_t eid)
  {
    domains_[domain].add_value(gid, eid);
  }

  void add_value(size_t domain, size_t gid, size_t eid, const dataset_t &d, common::rng_t &rng);
  size_t remove_value(size_t domain, size_t eid, const dataset_t &d, common::rng_t &rng);

  std::pair<std::vector<size_t>, std::vector<float>>
  score_value(size_t domain, size_t eid, const dataset_t &d, common::rng_t &rng);

  inline const std::map<tuple_t, group_with_count_t> &
  get_suff_stats(size_t relation) const
  {
    return relations_[relation].second;
  }

private:

  // for the given domain, return a list of (relation, position) pairs
  // for each relation it is a part of
  inline std::vector<std::pair<size_t, size_t>>
  domain_relations(size_t d) const
  {
    std::vector<std::pair<size_t, size_t>> ret;
    for (size_t i = 0; i < relations_.size(); i++) {
      const auto &ds = relations_[i].first.domains_;
      const auto it = std::find(ds.begin(), ds.end(), d);
      if (it != ds.end())
        ret.emplace_back(i, (it - ds.begin()));
    }
    return ret;
  }

  std::vector<domain> domains_;
  std::vector<std::vector<std::pair<size_t, size_t>>> domain_relations_;

  std::vector<
    std::pair<relation_t, std::map<tuple_t, group_with_count_t>>
  > relations_;

  float running_score_;
};

} // namespace irm
} // namespace microscopes
