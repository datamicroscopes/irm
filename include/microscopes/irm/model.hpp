#pragma once

#include <microscopes/common/relation/dataview.hpp>
#include <microscopes/common/entity_state.hpp>
#include <microscopes/common/group_manager.hpp>
#include <microscopes/common/typedefs.hpp>
#include <microscopes/common/macros.hpp>
#include <microscopes/common/assert.hpp>
#include <microscopes/common/util.hpp>
#include <microscopes/common/static_vector.hpp>
#include <microscopes/models/base.hpp>
#include <microscopes/io/schema.pb.h>

#include <distributions/special.hpp>

#include <cmath>
#include <vector>
#include <set>
#include <functional>
#include <map>
#include <memory>
#include <sstream>
#include <utility>
#include <stdexcept>
#include <algorithm>

namespace microscopes {
namespace irm {

class relation_definition {
public:
  relation_definition() = default;
  relation_definition(const std::vector<size_t> &domains,
                      const std::shared_ptr<models::model> &model)
    : domains_(domains), model_(model)
  {
    MICROSCOPES_DCHECK(domains.size(), "arity impossible");
    MICROSCOPES_DCHECK(model.get(), "nullptr model");
  }
  inline const std::vector<size_t> & domains() const { return domains_; }
  inline const std::shared_ptr<models::model> & model() const { return model_; }
  inline size_t arity() const { return domains_.size(); }
private:
  std::vector<size_t> domains_;
  std::shared_ptr<models::model> model_;
};

class model_definition {
public:
  model_definition(const std::vector<size_t> &domains,
                   const std::vector<relation_definition> &relations);
  inline const std::vector<size_t> & domains() const { return domains_; }
  inline const std::vector<relation_definition> & relations() const { return relations_; }
private:
  std::vector<size_t> domains_;
  std::vector<relation_definition> relations_;
};

namespace detail {
struct _empty {};
typedef common::group_manager<_empty> domain;

template <typename T, ssize_t M, bool UseStatic>
struct vector_type_selector_impl {};

template <typename T, ssize_t M>
struct vector_type_selector_impl<T, M, true> {
  static_assert(M >= 0, "");
  typedef common::static_vector<T, size_t(M)> type;

  static inline type
  from_variadic(const std::vector<T> &v)
  {
    return type(v.begin(), v.end());
  }
};

template <typename T, ssize_t M>
struct vector_type_selector_impl<T, M, false> {
  typedef std::vector<T> type;

  static inline type from_variadic(const std::vector<T> &t) { return t; }
};

template <typename T, ssize_t M>
struct vector_type_selector {
  static const bool using_static = (M > 0);
  typedef typename vector_type_selector_impl<T, M, using_static>::type type;

  static inline type
  from_variadic(const std::vector<T> &t)
  {
    return vector_type_selector_impl<T, M, using_static>::from_variadic(t);
  }
};


} // namespace detail

typedef std::vector<const common::relation::dataview *> dataset_t;
typedef detail::domain domain;

template <ssize_t MaxRelationArity = -1>
class state {

  static_assert(MaxRelationArity == -1 || MaxRelationArity >= 2,
                "Invalid MaxRelationArity, either -1 or >= 2");

  template <ssize_t> friend class model;

public:

  typedef typename detail::vector_type_selector<size_t, MaxRelationArity>::type tuple_t;
  typedef std::vector<size_t> variadic_tuple_t;

  struct suffstats_t {
    suffstats_t() : ident_(), count_(), ss_() {}
    common::ident_t ident_; // an identifier for outside naming
    unsigned count_; // a ref count, so we know when to remove
    // XXX: unique_ptr instead?
    std::shared_ptr<models::group> ss_;
  };

  struct relation_container_t {
    relation_container_t()
      : desc_(), hypers_(),
        suffstats_table_(), ident_table_(), ident_gen_() {}
    relation_container_t(const relation_definition &desc)
      : desc_(desc), hypers_(desc.model()->create_hypers()),
        suffstats_table_(), ident_table_(), ident_gen_()
    {
      if (MaxRelationArity != -1)
        MICROSCOPES_DCHECK(
            desc.arity() <= MaxRelationArity,
            "cannot handle arity");
    }

    void
    dump(io::IrmRelation &r) const
    {
      r.set_hypers(hypers_->get_hp());
      for (const auto &p : suffstats_table_) {
        io::IrmSuffstat &ss = *r.add_suffstats();
        for (auto gid : p.first)
          ss.add_gids(gid);
        ss.set_id(p.second.ident_);
        ss.set_count(p.second.count_);
        ss.set_suffstat(p.second.ss_->get_ss());
      }
    }

    relation_definition desc_;
    // XXX: unique_ptr instead?
    std::shared_ptr<models::hypers> hypers_;
    std::map<tuple_t, suffstats_t> suffstats_table_;
    std::map<common::ident_t, tuple_t> ident_table_;
    common::ident_t ident_gen_;
  };

  state(const std::vector<domain> &domains,
        const std::vector<relation_container_t> &relations)
    : domains_(domains), relations_(relations)
  {
    domain_relations_.reserve(domains_.size());
    for (size_t i = 0; i < domains_.size(); i++)
      domain_relations_.emplace_back(domain_relations(i));
  }

  inline size_t
  ndomains() const
  {
    return domains_.size();
  }

  inline size_t nrelations() const { return relations_.size(); }

  inline size_t
  nentities(size_t domain) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].nentities();
  }

  inline size_t
  ngroups(size_t domain) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].ngroups();
  }

  inline const std::vector<ssize_t> &
  assignments(size_t domain) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].assignments();
  }

  inline std::vector<size_t>
  groups(size_t domain) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].groups();
  }

  inline const std::set<size_t> &
  empty_groups(size_t domain) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].empty_groups();
  }

  inline bool
  isactivegroup(size_t domain, size_t gid) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].isactivegroup(gid);
  }

  inline size_t
  groupsize(size_t domain, size_t gid) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].groupsize(gid);
  }

  inline common::hyperparam_bag_t
  get_domain_hp(size_t domain) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain id");
    return domains_[domain].get_hp();
  }

  inline void
  set_domain_hp(size_t domain, const common::hyperparam_bag_t &hp)
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain id");
    domains_[domain].set_hp(hp);
  }

  inline common::value_mutator
  get_domain_hp_mutator(size_t domain, const std::string &key)
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain id");
    return domains_[domain].get_hp_mutator(key);
  }

  inline common::hyperparam_bag_t
  get_relation_hp(size_t relation) const
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    return relations_[relation].hypers_->get_hp();
  }

  inline void
  set_relation_hp(size_t relation, const common::hyperparam_bag_t &hp)
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    relations_[relation].hypers_->set_hp(hp);
  }

  inline void
  set_relation_hp(size_t relation, const models::hypers &proto)
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    relations_[relation].hypers_->set_hp(proto);
  }

  inline common::value_mutator
  get_relation_hp_mutator(size_t relation, const std::string &key)
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    return relations_[relation].hypers_->get_hp_mutator(key);
  }

  inline std::vector<common::ident_t>
  suffstats_identifiers(size_t relation)
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    const auto &tab = relations_[relation].ident_table_;
    std::vector<common::ident_t> ret;
    ret.reserve(tab.size());
    for (const auto &p : tab)
      ret.push_back(p.first);
    return ret;
  }

  inline common::suffstats_bag_t
  get_suffstats(size_t relation, common::ident_t id) const
  {
    return get_suffstats_t(relation, id).ss_->get_ss();
  }

  inline bool
  get_suffstats(size_t relation, const variadic_tuple_t &gids, common::suffstats_bag_t &ss) const
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    const auto &gids1 =
      detail::vector_type_selector<size_t, MaxRelationArity>::from_variadic(gids);
    const auto it = relations_[relation].suffstats_table_.find(gids1);
    if (it == relations_[relation].suffstats_table_.end())
      return false;
    ss = it->second.ss_->get_ss();
    return true;
  }

  inline void
  set_suffstats(size_t relation, common::ident_t id, const common::suffstats_bag_t &ss)
  {
    get_suffstats_t(relation, id).ss_->set_ss(ss);
  }

  inline common::value_mutator
  get_suffstats_mutator(size_t relation, common::ident_t id, const std::string &key)
  {
    return get_suffstats_t(relation, id).ss_->get_ss_mutator(key);
  }

  inline size_t
  get_suffstats_count(size_t relation, common::ident_t id) const
  {
    return get_suffstats_t(relation, id).count_;
  }

  inline size_t
  create_group(size_t domain)
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain id");
    return domains_[domain].create_group().first;
  }

  inline void
  delete_group(size_t domain, size_t gid)
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain id");
    // XXX: see note in remove_value_from_feature_group()
    // XXX: our GC strategy is currently inefficient
    for (const auto &dr : domain_relations_[domain]) {
      auto &relation = relations_[dr.rel_];
      auto it = relation.suffstats_table_.begin();
      while (it != relation.suffstats_table_.end()) {
        if (it->first[dr.pos_] != gid) {
          ++it;
          continue;
        }
        MICROSCOPES_ASSERT(!it->second.count_);
        relation.ident_table_.erase(it->second.ident_);
        relation.suffstats_table_.erase(it++); // must use postfix add
      }
    }
    domains_[domain].delete_group(gid);
  }

  inline void
  add_value(size_t domain, size_t gid, size_t eid, const dataset_t &d, common::rng_t &rng)
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    assert_correct_shape(d);
    add_value0(domain, gid, eid, d, rng, nullptr);
  }

  inline size_t
  remove_value(size_t domain, size_t eid, const dataset_t &d, common::rng_t &rng)
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    assert_correct_shape(d);
    return remove_value0(domain, eid, d, rng);
  }

  inline std::pair<std::vector<size_t>, std::vector<float>>
  score_value(size_t domain, size_t eid, const dataset_t &d, common::rng_t &rng) const
  {
    std::pair<std::vector<size_t>, std::vector<float>> ret;
    inplace_score_value(ret, domain, eid, d, rng);
    return ret;
  }

  inline void
  inplace_score_value(
      std::pair<std::vector<size_t>, std::vector<float>> &scores,
      size_t domain, size_t eid, const dataset_t &d, common::rng_t &rng) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    assert_correct_shape(d);
    inplace_score_value0(scores, domain, eid, d, rng);
  }

  inline float
  score_assignment(size_t domain) const
  {
    MICROSCOPES_DCHECK(domain < domains_.size(), "invalid domain");
    return domains_[domain].score_assignment();
  }

  inline float
  score_assignment() const
  {
    float score = 0.;
    for (auto &d : domains_)
      score += d.score_assignment();
    return score;
  }

  inline float
  score_likelihood(size_t relation, common::ident_t id, common::rng_t &rng) const
  {
    auto &ss = get_suffstats_t(relation, id);
    return ss.ss_->score_data(*relations_[relation].hypers_, rng);
  }

  inline float
  score_likelihood(size_t relation, common::rng_t &rng) const
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    float score = 0.;
    auto &m = relations_[relation].hypers_;
    for (auto &p : relations_[relation].suffstats_table_)
      score += p.second.ss_->score_data(*m, rng);
    return score;
  }

  inline float
  score_likelihood(common::rng_t &rng) const
  {
    float score = 0.;
    for (size_t i = 0; i < relations_.size(); i++)
      score += score_likelihood(i, rng);
    return score;
  }

  inline void
  assert_correct_shape(const dataset_t &d) const
  {
    MICROSCOPES_DCHECK(d.size() == relations_.size(), "#s dont match");
    for (size_t i = 0; i < d.size(); i++) {
      MICROSCOPES_DCHECK(relations_[i].desc_.domains().size() == d[i]->dims(), "arity does not match");
      for (size_t j = 0; j < d[i]->dims(); j++)
        MICROSCOPES_DCHECK(domains_[relations_[i].desc_.domains()[j]].nentities() == d[i]->shape()[j], "shape does not match");
    }
  }

  // for testing
  std::vector<variadic_tuple_t>
  entity_data_positions(size_t domain, size_t eid, const dataset_t &d) const
  {
    std::vector<variadic_tuple_t> ret;
    iterate_over_entity_data(
        domain, eid, d,
        [&ret](size_t, const variadic_tuple_t &eids, const common::value_accessor &) {
          ret.emplace_back(eids);
        });
    return ret;
  }

  std::string
  serialize() const
  {
    io::IrmState m;
    for (const auto &d : domains_)
      m.add_domains(d.serialize([](const detail::_empty &) {return "";}));
    for (const auto &r : relations_) {
      io::IrmRelation &msg = *m.add_relations();
      r.dump(msg);
    }
    return common::util::protobuf_to_string(m);
  }

  /**
   * initialized to an **invalid** point in the state space!
   *
   * within each domain:
   *   (A) no entities assigned
   *   (B) no groups
   *   (C) no hypers initialized
   *
   * useful primarily for testing purposes
   */
  static std::shared_ptr<state>
  unsafe_initialize(const model_definition &defn)
  {
    std::vector<domain> domains;
    domains.reserve(defn.domains().size());
    for (auto n : defn.domains())
      domains.emplace_back(n);
    std::vector<relation_container_t> relations;
    relations.reserve(defn.relations().size());
    for (const auto &r : defn.relations()) {
      relations.emplace_back();
      auto &reln = relations.back();
      reln.desc_ = r;
      reln.hypers_ = r.model()->create_hypers();
    }
    return std::make_shared<state>(domains, relations);
  }

  static std::shared_ptr<state>
  initialize(const model_definition &defn,
             const std::vector<common::hyperparam_bag_t> &cluster_inits,
             const std::vector<common::hyperparam_bag_t> &relation_inits,
             const std::vector<variadic_tuple_t> &domain_assignments,
             const dataset_t &data,
             common::rng_t &rng);

  static std::shared_ptr<state>
  deserialize(const model_definition &defn,
              const common::serialized_t &s);

protected:
  // the *_value0 methods do no error checking

  inline void
  add_value0(
      size_t domain, size_t gid, size_t eid,
      const dataset_t &d, common::rng_t &rng, float *acc_score)
  {
    domains_[domain].add_value(gid, eid);
    tuple_t gids;
    iterate_over_entity_data(
        domain, eid, d,
        [this, &gids, &rng, acc_score](
          size_t rid,
          const variadic_tuple_t &eids,
          const common::value_accessor &value)
        {
          auto &relation = this->relations_[rid];
          this->eids_to_gids_under_relation(gids, eids, relation.desc_);
          this->add_value_to_feature_group(gids, value, relation, rng, acc_score);
        });
  }

  inline size_t
  remove_value0(
      size_t domain, size_t eid,
      const dataset_t &d, common::rng_t &rng)
  {
    tuple_t gids;
    iterate_over_entity_data(
        domain, eid, d,
        [this, &gids, &rng](
           size_t rid,
           const variadic_tuple_t &eids,
           const common::value_accessor &value)
        {
            auto &relation = this->relations_[rid];
            this->eids_to_gids_under_relation(gids, eids, relation.desc_);
            this->remove_value_from_feature_group(gids, value, relation, rng);
        });
    return domains_[domain].remove_value(eid).first;
  }

  void
  inplace_score_value0(
      std::pair<std::vector<size_t>, std::vector<float>> &scores,
      size_t did,
      size_t eid,
      const dataset_t &d,
      common::rng_t &rng) const
  {
    using distributions::fast_log;

    const auto &domain = domains_[did];
    MICROSCOPES_DCHECK(!domain.empty_groups().empty(), "no empty groups");

    scores.first.clear();
    scores.second.clear();
    scores.first.reserve(domain.ngroups());
    scores.second.reserve(domain.ngroups());

    float pseudocounts = 0;
    for (const auto &g : domain) {
      const float pseudocount = domain.pseudocount(g.first, g.second);
      float sum = fast_log(pseudocount);
      const_cast<state *>(this)->add_value0(did, g.first, eid, d, rng, &sum);
      const size_t gid = const_cast<state *>(this)->remove_value0(did, eid, d, rng);
      if (unlikely(gid != g.first))
        MICROSCOPES_ASSERT(false);
      scores.first.push_back(g.first);
      scores.second.push_back(sum);
      pseudocounts += pseudocount;
    }

    const float lgnorm = fast_log(pseudocounts);
    for (auto &s : scores.second)
      s -= lgnorm;
  }

private:

  struct rel_pos_t {
    rel_pos_t() : rel_(), pos_() {}
    rel_pos_t(size_t rel, size_t pos) : rel_(rel), pos_(pos) {}
    size_t rel_;
    size_t pos_;
  };

  inline void
  eids_to_gids_under_relation(
      tuple_t &gids,
      const variadic_tuple_t &eids,
      const relation_definition &desc) const
  {
    gids.clear();
    gids.reserve(desc.domains().size());
    for (size_t i = 0; i < desc.domains().size(); i++) {
      MICROSCOPES_DCHECK(
          domains_[desc.domains()[i]].assignments()[eids[i]] != -1,
          "eid is not assigned to a valid group");
      gids.push_back(domains_[desc.domains()[i]].assignments()[eids[i]]);
    }
  }

  void
  add_value_to_feature_group(
      const tuple_t &gids,
      const common::value_accessor &value,
      relation_container_t &relation,
      common::rng_t &rng,
      float *acc_score)
  {
    MICROSCOPES_ASSERT(!value.anymasked());
    std::shared_ptr<models::group> group;
    auto it = relation.suffstats_table_.find(gids);
    if (it == relation.suffstats_table_.end()) {
      group = relation.hypers_->create_group(rng);
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
      *acc_score += group->score_value(*relation.hypers_, value, rng);
    group->add_value(*relation.hypers_, value, rng);
  }

  void
  remove_value_from_feature_group(
      const tuple_t &gids,
      const common::value_accessor &value,
      relation_container_t &relation,
      common::rng_t &rng)
  {
    auto it = relation.suffstats_table_.find(gids);
    MICROSCOPES_ASSERT(!value.anymasked());
    MICROSCOPES_ASSERT(it != relation.suffstats_table_.end());
    MICROSCOPES_ASSERT(it->second.count_);
    MICROSCOPES_ASSERT(it->second.ss_);
    MICROSCOPES_ASSERT(
        relation.ident_table_.find(it->second.ident_) != relation.ident_table_.end() &&
        relation.ident_table_[it->second.ident_] == gids);
    it->second.ss_->remove_value(*relation.hypers_, value, rng);
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

  template <typename T>
  void
  iterate_over_entity_data(
      size_t domain,
      size_t eid,
      const dataset_t &d,
      T callback) const
  {
    tuple_t ignore_idxs;
    for (const auto &dr : domain_relations_[domain]) {
      auto &relation = relations_[dr.rel_];
      auto &data = d[dr.rel_];
      ignore_idxs.clear();
      for (size_t i = 0; i < dr.pos_; i++)
        if (relation.desc_.domains()[i] == domain)
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

  inline suffstats_t &
  get_suffstats_t(size_t relation, common::ident_t id)
  {
    MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
    auto &rel = relations_[relation];
    auto &tab = rel.ident_table_;
    auto it = tab.find(id);
    MICROSCOPES_DCHECK(it != tab.end(), "invalid ident");
    auto it1 = rel.suffstats_table_.find(it->second);
    MICROSCOPES_ASSERT(it1 != rel.suffstats_table_.end());
    return it1->second;
  }

  inline const suffstats_t &
  get_suffstats_t(size_t relation, common::ident_t id) const
  {
    return const_cast<state *>(this)->get_suffstats_t(relation, id);
  }

  // for the given domain, return a list of (relation, position) pairs
  // for each relation it is a part of
  inline std::vector<rel_pos_t>
  domain_relations(size_t d) const
  {
    std::vector<rel_pos_t> ret;
    for (size_t r = 0; r < relations_.size(); r++) {
      const auto &ds = relations_[r].desc_.domains();
      for (size_t pos = 0; pos < ds.size(); pos++)
        if (ds[pos] == d)
          ret.emplace_back(r, pos);
    }
    return ret;
  }

  std::vector<domain> domains_;
  std::vector<std::vector<rel_pos_t>> domain_relations_;
  std::vector<relation_container_t> relations_;
};

template <ssize_t MaxRelationArity>
std::shared_ptr<state<MaxRelationArity>>
state<MaxRelationArity>::initialize(
    const model_definition &defn,
    const std::vector<common::hyperparam_bag_t> &cluster_inits,
    const std::vector<common::hyperparam_bag_t> &relation_inits,
    const std::vector<variadic_tuple_t> &domain_assignments,
    const dataset_t &data,
    common::rng_t &rng)
{
  MICROSCOPES_DCHECK(cluster_inits.size() == defn.domains().size(),
      "# domains mismatch");
  MICROSCOPES_DCHECK(relation_inits.size() == defn.relations().size(),
      "# relations mismatch");
  MICROSCOPES_DCHECK(domain_assignments.size() == defn.domains().size(),
      "# domains mismatch");
  MICROSCOPES_DCHECK(data.size() == defn.relations().size(),
      "# relations mismatch");
  auto p = unsafe_initialize(defn);
  for (size_t i = 0; i < cluster_inits.size(); i++)
    p->domains_[i].set_hp(cluster_inits[i]);
  for (size_t i = 0; i < relation_inits.size(); i++)
    p->relations_[i].hypers_->set_hp(relation_inits[i]);

  for (size_t i = 0; i < domain_assignments.size(); i++) {
    auto assignment = domain_assignments[i];
    if (assignment.empty())
      assignment = common::util::random_assignment_vector(defn.domains()[i], rng);
    const size_t ngroups =
      *std::max_element(assignment.begin(), assignment.end()) + 1;
    for (size_t j = 0; j < ngroups; j++)
      p->domains_[i].create_group();
    for (size_t j = 0; j < assignment.size(); j++)
      p->domains_[i].add_value(assignment[j], j);
  }

#ifdef DEBUG_MODE
  for (auto &d : p->domains_)
    for (auto s : d.assignments())
      MICROSCOPES_DCHECK(s != -1, "assignments should all be filled");
#endif

  tuple_t gids;
  for (size_t i = 0; i < p->relations_.size(); i++) {
    auto &relation = p->relations_[i];
    for (size_t outer = 0; outer < data[i]->shape().front(); outer++)
      for (const auto &pp : data[i]->slice(0, outer)) {
        p->eids_to_gids_under_relation(gids, pp.first, relation.desc_);
        p->add_value_to_feature_group(gids, pp.second, relation, rng, nullptr);
      }
  }
  return p;
}

template <ssize_t MaxRelationArity>
std::shared_ptr<state<MaxRelationArity>>
state<MaxRelationArity>::deserialize(
    const model_definition &defn,
    const common::serialized_t &s)
{
  common::rng_t rng; // XXX: hack

  // some attempt made to validate inputs, but not foolproof
  io::IrmState m;
  common::util::protobuf_from_string(m, s);

  MICROSCOPES_DCHECK((size_t)m.domains_size() == defn.domains().size(),
      "# domains mismatch");
  std::vector<domain> domains;
  domains.reserve(defn.domains().size());
  for (size_t i = 0; i < defn.domains().size(); i++)
    domains.emplace_back(
        m.domains(i), [](const std::string &) { return detail::_empty(); });

  MICROSCOPES_DCHECK((size_t)m.relations_size() == defn.relations().size(),
      "# relations mismatch");
  std::vector<relation_container_t> relations;
  relations.reserve(defn.relations().size());
  for (size_t i = 0; i < defn.relations().size(); i++) {
    const auto &rdef = defn.relations()[i];
    relation_container_t reln;
    reln.desc_ = rdef;
    const auto &r = m.relations(i);
    reln.hypers_ = rdef.model()->create_hypers();
    reln.hypers_->set_hp(r.hypers());

    for (size_t j = 0; j < size_t(r.suffstats_size()); j++) {
      suffstats_t suffstat;
      const auto &ss = r.suffstats(j);
      MICROSCOPES_DCHECK((size_t)ss.gids_size() == rdef.domains().size(),
          "arity mismatch");
      tuple_t gids;
      gids.reserve(ss.gids_size());
      for (size_t k = 0; k < rdef.domains().size(); k++)
        gids.push_back(ss.gids(k));
      // see note in schema.proto: not validated
      suffstat.ident_ = ss.id();
      suffstat.count_ = ss.count();
      suffstat.ss_ = reln.hypers_->create_group(rng);
      suffstat.ss_->set_ss(ss.suffstat());

      reln.suffstats_table_[gids] = suffstat;
      reln.ident_table_[ss.id()] = gids;
      reln.ident_gen_ = std::max<size_t>(reln.ident_gen_, ss.id() + 1);
    }

    relations.emplace_back(reln);
  }

  return std::make_shared<state<MaxRelationArity>>(domains, relations);
}

/**
 * The binds happen on a per-domain basis
 */
template <ssize_t MaxRelationArity = -1>
class model : public common::entity_based_state_object {
public:
  model(const std::shared_ptr<state<MaxRelationArity>> &impl,
        size_t domain,
        const std::vector<std::shared_ptr<common::relation::dataview>> &data)
    : impl_(impl), domain_(domain), data_(data), data_raw_()
  {
    MICROSCOPES_DCHECK(impl.get(), "nullptr impl");
    data_raw_.reserve(data_.size());
    for (auto &p : data_) {
      MICROSCOPES_DCHECK(p.get(), "nullptr dataview");
      data_raw_.push_back(p.get());
    }
    impl_->assert_correct_shape(data_raw_);
  }

  size_t nentities() const override { return impl_->nentities(domain_); }
  size_t ngroups() const override { return impl_->ngroups(domain_); }
  size_t ncomponents() const override { return impl_->nrelations(); }

  std::vector<ssize_t> assignments() const override { return impl_->assignments(domain_); }
  std::vector<size_t> groups() const override { return impl_->groups(domain_); }

  std::vector<size_t>
  empty_groups() const override
  {
    const auto &egs = impl_->empty_groups(domain_);
    return std::vector<size_t>(egs.begin(), egs.end());
  }

  size_t groupsize(size_t gid) const override { return impl_->groupsize(domain_, gid); }

  common::hyperparam_bag_t get_cluster_hp() const override { return impl_->get_domain_hp(domain_); }
  void set_cluster_hp(const common::hyperparam_bag_t &hp) override { impl_->set_domain_hp(domain_, hp); }
  common::value_mutator get_cluster_hp_mutator(const std::string &key) override { return impl_->get_domain_hp_mutator(domain_, key); }

  common::hyperparam_bag_t get_component_hp(size_t component) const override { return impl_->get_relation_hp(component); }
  void set_component_hp(size_t component, const common::hyperparam_bag_t &hp) override { impl_->set_relation_hp(component, hp); }
  void set_component_hp(size_t component, const models::hypers &proto) override { impl_->set_relation_hp(component, proto); }
  common::value_mutator get_component_hp_mutator(size_t component, const std::string &key) override { return impl_->get_relation_hp_mutator(component, key); }

  std::vector<common::ident_t>
  suffstats_identifiers(size_t component) const override
  {
    return impl_->suffstats_identifiers(component);
  }

  common::suffstats_bag_t
  get_suffstats(size_t component, common::ident_t id) const override
  {
    return impl_->get_suffstats(component, id);
  }

  void
  set_suffstats(size_t component, common::ident_t id, const common::suffstats_bag_t &ss) override
  {
    impl_->set_suffstats(component, id, ss);
  }

  common::value_mutator get_suffstats_mutator(size_t component, common::ident_t id, const std::string &key) override { return impl_->get_suffstats_mutator(component, id, key); }

  void
  add_value(size_t gid, size_t eid, common::rng_t &rng) override
  {
    impl_->add_value0(domain_, gid, eid, data_raw_, rng, nullptr);
  }

  size_t
  remove_value(size_t eid, common::rng_t &rng) override
  {
    return impl_->remove_value0(domain_, eid, data_raw_, rng);
  }

  std::pair<std::vector<size_t>, std::vector<float>>
  score_value(size_t eid, common::rng_t &rng) const override
  {
    return impl_->score_value(domain_, eid, data_raw_, rng);
  }

  void
  inplace_score_value(
      std::pair<std::vector<size_t>, std::vector<float>> &scores,
      size_t eid, common::rng_t &rng) const override
  {
    impl_->inplace_score_value0(scores, domain_, eid, data_raw_, rng);
  }

  float score_assignment() const override { return impl_->score_assignment(domain_); }

  float
  score_likelihood(size_t component, common::ident_t id, common::rng_t &rng) const override
  {
    return impl_->score_likelihood(component, id, rng);
  }

  float
  score_likelihood(size_t component, common::rng_t &rng) const override
  {
    return impl_->score_likelihood(component, rng);
  }

  size_t create_group(common::rng_t &rng) override { return impl_->create_group(domain_); }

  void delete_group(size_t gid) override { impl_->delete_group(domain_, gid); }

private:
  std::shared_ptr<state<MaxRelationArity>> impl_;
  size_t domain_;
  std::vector<std::shared_ptr<common::relation::dataview>> data_;
  std::vector<const common::relation::dataview *> data_raw_;
};

// template instantiations
extern template class state<-1>;
extern template class state<2>;
extern template class state<3>;
extern template class state<4>;

extern template class model<-1>;
extern template class model<2>;
extern template class model<3>;
extern template class model<4>;

} // namespace irm
} // namespace microscopes
