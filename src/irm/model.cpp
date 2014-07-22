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
using namespace microscopes::io;
using namespace microscopes::models;

//static inline void
//AssertAllEntitiesAccounted(const vector<set<size_t>> &c, size_t n)
//{
//  vector<bool> ents(n, false);
//  for (const auto &s : c)
//    for (auto eid : s) {
//      MICROSCOPES_DCHECK(eid < ents.size(), "bad eid given");
//      MICROSCOPES_DCHECK(!ents[eid], "eid given twice");
//      ents[eid] = true;
//    }
//  for (auto b : ents)
//    MICROSCOPES_DCHECK(b, "ent unaccounted for");
//}

static void
AssertRelationDefinitionPossible(
  const model_definition &model,
  const relation_definition &relation)
{
  for (auto d : relation.domains())
    MICROSCOPES_DCHECK(d < model.domains().size(),
        "invalid domain given");
}

model_definition::model_definition(
    const vector<size_t> &domains,
    const vector<relation_definition> &relations)
  : domains_(domains), relations_(relations)
{
  MICROSCOPES_DCHECK(domains.size(), "no domains given");
  MICROSCOPES_DCHECK(relations.size(), "no relations given");
  for (auto s : domains)
    MICROSCOPES_DCHECK(s, "empty domain given");
  for (const auto &r : relations) {
    for (auto d : r.domains())
      MICROSCOPES_DCHECK(d < domains.size(), "invalid domain given");
  }
}

state::state(const model_definition &defn,
             const vector<domain> &domains,
             const vector<relation_container_t> &relations)
  : domains_(domains), relations_(relations)
{
  // very little effort made to validate inputs
  MICROSCOPES_DCHECK(defn.domains().size() == domains.size(),
      "# domains mismatch");
  MICROSCOPES_DCHECK(defn.relations().size() == relations.size(),
      "# relations mismatch");
  domain_relations_.reserve(domains_.size());
  for (size_t i = 0; i < domains_.size(); i++)
    domain_relations_.emplace_back(domain_relations(i));
  for (size_t i = 0; i < defn.domains().size(); i++)
    MICROSCOPES_DCHECK(defn.domains()[i] == domains[i].nentities(),
      "# entities mismatch");
  for (const auto &r : relations)
    AssertRelationDefinitionPossible(defn, r.desc_);
}

void
state::eids_to_gids_under_relation(
    vector<size_t> &gids,
    const vector<size_t> &eids,
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
state::add_value_to_feature_group(
    const vector<size_t> &gids,
    const value_accessor &value,
    relation_container_t &relation,
    rng_t &rng,
    float *acc_score)
{
  MICROSCOPES_ASSERT(!value.anymasked());
  shared_ptr<group> group;
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

//void
//state::initialize(const vector<vector<set<size_t>>> &clusters, const dataset_t &d, rng_t &rng)
//{
//  MICROSCOPES_DCHECK(clusters.size() == domains_.size(), "invalid number of clusterings");
//  assert_correct_shape(d);
//
//  for (size_t i = 0; i < domains_.size(); i++) {
//    MICROSCOPES_DCHECK(!domains_[i].ngroups(), "domain not empty");
//    const auto &c = clusters[i];
//    MICROSCOPES_DCHECK(c.size() > 0, "no clusters given!");
//    for (size_t j = 0; j < c.size(); j++)
//      domains_[i].create_group();
//    AssertAllEntitiesAccounted(c, domains_[i].nentities());
//    for (size_t j = 0; j < c.size(); j++)
//      for (auto eid : c[j])
//        domains_[i].add_value(j, eid);
//  }
//#ifdef DEBUG_MODE
//  for (auto &d : domains_)
//    for (auto s : d.assignments())
//      MICROSCOPES_DCHECK(s != -1, "assignments should all be filled");
//#endif
//
//  vector<size_t> gids;
//  for (size_t i = 0; i < relations_.size(); i++) {
//    auto &relation = relations_[i];
//    for (const auto &p : *d[i]) {
//      eids_to_gids_under_relation(gids, p.first, relation.desc_);
//      add_value_to_feature_group(gids, p.second, relation, rng, nullptr);
//    }
//  }
//}
//
//void
//state::random_initialize(const dataset_t &d, rng_t &rng)
//{
//  assert_correct_shape(d);
//  vector<vector<set<size_t>>> clusters;
//  clusters.reserve(domains_.size());
//  for (auto &d : domains_) {
//    vector<set<size_t>> cluster;
//    // create min(100, n/2) + 1 groups
//    const size_t ngroups = min(size_t(100), d.nentities()) + 1;
//    cluster.resize(ngroups);
//    const auto groups = util::range(ngroups);
//    for (size_t i = 0; i < d.nentities(); i++) {
//      const auto choice = util::sample_choice(groups, rng);
//      cluster[choice].insert(i);
//    }
//    clusters.emplace_back(cluster);
//  }
//  initialize(clusters, d, rng);
//}

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
  return ss.ss_->score_data(*relations_[relation].hypers_, rng);
}

float
state::score_likelihood(size_t relation, rng_t &rng) const
{
  MICROSCOPES_DCHECK(relation < relations_.size(), "invalid relation id");
  float score = 0.;
  auto &m = relations_[relation].hypers_;
  for (auto &p : relations_[relation].suffstats_table_)
    score += p.second.ss_->score_data(*m, rng);
  return score;
}

void
state::relation_container_t::dump(IrmRelation &r) const
{
  r.set_hypers(hypers_->get_hp());
  for (const auto &p : suffstats_table_) {
    IrmSuffstat &ss = *r.add_suffstats();
    for (auto gid : p.first)
      ss.add_gids(gid);
    ss.set_id(p.second.ident_);
    ss.set_count(p.second.count_);
    ss.set_suffstat(p.second.ss_->get_ss());
  }
}

string
state::serialize() const
{
  IrmState m;
  for (const auto &d : domains_)
    m.add_domains(d.serialize([](const _empty &) {return "";}));
  for (const auto &r : relations_) {
    IrmRelation &msg = *m.add_relations();
    r.dump(msg);
  }
  return util::protobuf_to_string(m);
}

shared_ptr<state>
state::unsafe_initialize(const model_definition &defn)
{
  vector<domain> domains;
  domains.reserve(defn.domains().size());
  for (auto n : defn.domains())
    domains.emplace_back(n);
  vector<relation_container_t> relations;
  relations.reserve(defn.relations().size());
  for (const auto &r : defn.relations()) {
    relations.emplace_back();
    auto &reln = relations.back();
    reln.desc_ = r;
    reln.hypers_ = r.model()->create_hypers();
  }
  return make_shared<state>(defn, domains, relations);
}

shared_ptr<state>
state::initialize(const model_definition &defn,
                  const vector<hyperparam_bag_t> &cluster_inits,
                  const vector<hyperparam_bag_t> &relation_inits,
                  const vector<vector<size_t>> &domain_assignments,
                  const dataset_t &data,
                  rng_t &rng)
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
      assignment = util::random_assignment_vector(defn.domains()[i], rng);
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

  vector<size_t> gids;
  for (size_t i = 0; i < p->relations_.size(); i++) {
    auto &relation = p->relations_[i];
    for (const auto &pp : *data[i]) {
      p->eids_to_gids_under_relation(gids, pp.first, relation.desc_);
      p->add_value_to_feature_group(gids, pp.second, relation, rng, nullptr);
    }
  }
  return p;
}

shared_ptr<state>
state::deserialize(const model_definition &defn, const serialized_t &s)
{
  rng_t rng; // XXX: hack

  // some attempt made to validate inputs, but not foolproof
  IrmState m;
  util::protobuf_from_string(m, s);

  MICROSCOPES_DCHECK(m.domains_size() == defn.domains().size(),
      "# domains mismatch");
  vector<domain> domains;
  domains.reserve(defn.domains().size());
  for (size_t i = 0; i < defn.domains().size(); i++)
    domains.emplace_back(
        m.domains(i), [](const string &) { return _empty(); });

  MICROSCOPES_DCHECK(m.relations_size() == defn.relations().size(),
      "# relations mismatch");
  vector<relation_container_t> relations;
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
      MICROSCOPES_DCHECK(ss.gids_size() == rdef.domains().size(),
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

      reln.suffstats_table_[gids] = move(suffstat);
      reln.ident_table_[ss.id()] = move(gids);
      reln.ident_gen_ = max<size_t>(reln.ident_gen_, ss.id() + 1);
    }

    relations.emplace_back(reln);
  }

  return make_shared<state>(defn, domains, relations);
}
