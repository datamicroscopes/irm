#include <microscopes/irm/model.hpp>
#include <microscopes/models/distributions.hpp>
#include <microscopes/common/relation/dataview.hpp>

#include <distributions/models/bb.hpp>

#include <memory>
#include <iostream>

using namespace std;
using namespace microscopes;
using namespace distributions;

static common::hyperparam_bag_t
crp_hp_messsage(float alpha)
{
  io::CRP m;
  m.set_alpha(alpha);
  return common::util::protobuf_to_string(m);
}

static common::hyperparam_bag_t
bb_hp_messsage(float alpha, float beta)
{
  distributions::protobuf::BetaBernoulli::Shared m;
  m.set_alpha(alpha);
  m.set_beta(beta);
  return common::util::protobuf_to_string(m);
}

static pair<
  shared_ptr<irm::state<4>>,
  vector<shared_ptr<common::relation::dataview>>
>
make_irm(size_t groups,
         size_t entities_per_group,
         size_t relations,
         common::rng_t &r)
{
  cout << "entities_per_group: " << entities_per_group << ", "
       << "relations: " << relations << ", "
       << "groups: " << groups << endl;

  const size_t n = groups * entities_per_group;

  irm::relation_definition reldef({0, 0},
      make_shared<models::distributions_model<BetaBernoulli>>());
  vector<irm::relation_definition> reldefs(relations, reldef);

  irm::model_definition defn({n}, reldefs);

  vector<common::hyperparam_bag_t> cluster_inits({crp_hp_messsage(1.)});
  vector<common::hyperparam_bag_t> relation_inits({bb_hp_messsage(1., 1.)});

  vector<size_t> assignment0;
  for (size_t i = 0; i < groups; i++)
    for (size_t j = 0; j < entities_per_group; j++)
      assignment0.push_back(i);

  // memory leaks
  bool * data = new bool[n * n];
  for (size_t i = 0; i < (n * n); i++)
    data[i] = bernoulli_distribution()(r);

  auto view = shared_ptr<common::relation::row_major_dense_dataview>(
      new common::relation::row_major_dense_dataview(
        reinterpret_cast<const uint8_t *>(data),
        nullptr,
        {n, n},
        common::runtime_type(TYPE_B)));

  auto latent = irm::state<4>::initialize(
      defn,
      cluster_inits,
      relation_inits,
      {assignment0},
      {view.get()},
      r);

  vector<shared_ptr<common::relation::dataview>> dataset({view});

  return make_pair(latent, dataset);
}

// for performance debugging purposes
// doesn't change the group assignments
void
perftest(common::entity_based_state_object &s, rng_t &rng)
{
  pair<vector<size_t>, vector<float>> scores;
  for (auto i : common::util::permute(s.nentities(), rng)) {
    const size_t gid = s.remove_value(i, rng);
    s.inplace_score_value(scores, i, rng);
    const auto choice = scores.first[common::util::sample_discrete_log(scores.second, rng)];
    (void)choice; // XXX: make sure compiler does not optimize this out
    s.add_value(gid, i, rng);
  }
}

int
main(int argc, char **argv)
{
  random_device rd;
  common::rng_t r(rd());
  auto p = make_irm(100, 100, 1, r);
  irm::model<4> m(p.first, 0, p.second);

  for (;;)
    perftest(m, r);

  return 0;
}
