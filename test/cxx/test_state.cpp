#include <microscopes/irm/model.hpp>
#include <microscopes/common/sparse_ndarray/dataview.hpp>
#include <microscopes/common/random_fwd.hpp>
#include <microscopes/models/distributions.hpp>

#include <random>
#include <iostream>

using namespace std;
using namespace distributions;
using namespace microscopes::common;
using namespace microscopes::common::sparse_ndarray;
using namespace microscopes::models;
using namespace microscopes::irm;

template <typename T>
static string
protobuf_to_string(const T &t)
{
  ostringstream out;
  t.SerializeToOstream(&out);
  return out.str();
}

static inline hyperparam_bag_t
beta_bernoulli_hp(float alpha, float beta)
{
  distributions_model<BetaBernoulli>::message_type m_feature_hp;
  m_feature_hp.set_alpha(1.0);
  m_feature_hp.set_beta(1.0);
  return protobuf_to_string(m_feature_hp);
}

static void
test1()
{
  random_device rd;
  rng_t r(rd());

  // 1 domain
  const vector<size_t> domains({10});

  // CRP priors
  const vector<float> domain_alphas({2.0});

  // 1 binary relation
  const vector<state::relation_t> relations({
      state::relation_t({0, 0}, distributions_factory<BetaBernoulli>().new_instance())
  });

  const vector<hyperparam_bag_t> relation_hps({
      beta_bernoulli_hp(2.0, 2.0)
  });

  state s(domains, relations);
  for (size_t i = 0; i < domain_alphas.size(); i++) {
    state::message_type hp;
    hp.set_alpha(domain_alphas[i]);
    s.set_domain_hp(i, protobuf_to_string(hp));
  }

  for (size_t i = 0; i < relation_hps.size(); i++)
    s.set_relation_hp(i, relation_hps[i]);

  // create fake, absolutely meaningless data for our relation

  bool *friends = new bool[domains[0]*domains[0]];
  bool *masks = new bool[domains[0]*domains[0]];

  size_t present = 0;
  for (size_t u = 0; u < domains[0]; u++) {
    for (size_t m = 0; m < domains[0]; m++) {
      const size_t idx = u*domains[0] + m;
      // coin flip to see if this data is present
      if (bernoulli_distribution(0.2)(r)) {
        masks[idx] = false;
        // coin flip to see if friends
        friends[idx] = bernoulli_distribution(0.8)(r);
        present++;
      } else {
        masks[idx] = true;
      }
    }
  }

  cout << "present: " << present << endl;

  unique_ptr<dataview> view(
    new row_major_dense_dataview(
        reinterpret_cast<uint8_t*>(friends), masks,
        {domains[0], domains[0]}, runtime_type(TYPE_B)));

  s.random_initialize({view.get()}, r);

  // peek @ suffstats
  size_t sum = 0;
  for (auto ident : s.suffstats_identifiers(0))
    sum += s.get_suffstats_count(0, ident);
  MICROSCOPES_DCHECK(sum == present, "suff stats don't match up");

  const size_t gid = s.remove_value(0, 0, {view.get()}, r);
  s.add_value(0, gid, 0, {view.get()}, r);

  sum = 0;
  for (auto ident : s.suffstats_identifiers(0))
    sum += s.get_suffstats_count(0, ident);
  MICROSCOPES_DCHECK(sum == present, "suff stats don't match up");
}

static void
test2()
{
  random_device rd;
  rng_t r(rd());

  // 10 users, 100 movies
  const vector<size_t> domains({10, 100});

  // CRP priors
  const vector<float> domain_alphas({2.0, 20.0});

  // a single binary relation between users x movies
  const vector<state::relation_t> relations({
      state::relation_t({0, 1}, distributions_factory<BetaBernoulli>().new_instance())
  });

  const vector<hyperparam_bag_t> relation_hps({
      beta_bernoulli_hp(2.0, 2.0)
  });

  state s(domains, relations);
  for (size_t i = 0; i < domain_alphas.size(); i++) {
    state::message_type hp;
    hp.set_alpha(domain_alphas[i]);
    s.set_domain_hp(i, protobuf_to_string(hp));
  }

  for (size_t i = 0; i < relation_hps.size(); i++)
    s.set_relation_hp(i, relation_hps[i]);

  // create fake, absolutely meaningless data for our relation

  bool *likes = new bool[domains[0]*domains[1]];
  bool *masks = new bool[domains[0]*domains[1]];

  size_t present = 0;
  for (size_t u = 0; u < domains[0]; u++) {
    for (size_t m = 0; m < domains[1]; m++) {
      const size_t idx = u*domains[1] + m;
      // coin flip to see if this data is present
      if (bernoulli_distribution(0.2)(r)) {
        masks[idx] = false;
        // coin flip to see if user likes
        likes[idx] = bernoulli_distribution(0.8)(r);
        present++;
      } else {
        masks[idx] = true;
      }
    }
  }

  cout << "present: " << present << endl;

  unique_ptr<dataview> view(
    new row_major_dense_dataview(
        reinterpret_cast<uint8_t*>(likes), masks,
        domains, runtime_type(TYPE_B)));

  s.random_initialize({view.get()}, r);

  size_t sum = 0;
  for (auto ident : s.suffstats_identifiers(0))
    sum += s.get_suffstats_count(0, ident);
  MICROSCOPES_DCHECK(sum == present, "suff stats don't match up");

  // score the 1st data point
  const size_t gid = s.remove_value(0, 0, {view.get()}, r);
  s.create_group(0);
  const auto scores = s.score_value(0, 0, {view.get()}, r);
  //cout << "scores: " << scores << endl;
  s.add_value(0, gid, 0, {view.get()}, r);

  // remove the data
  for (size_t u = 0; u < domains[0]; u++)
    s.remove_value(0, u, {view.get()}, r);

  // peek @ suffstats
  sum = 0;
  for (auto ident : s.suffstats_identifiers(0))
    sum += s.get_suffstats_count(0, ident);
  MICROSCOPES_DCHECK(!sum, "suff stats don't match up");

  delete [] likes;
  delete [] masks;
}

static pair<unique_ptr<bool[]>, unique_ptr<bool[]>>
binary_relation_generate(
    size_t domain0, size_t domain1,
    float p0, float p1,
    rng_t &r)
{
  unique_ptr<bool[]> data(new bool[domain0*domain1]);
  unique_ptr<bool[]> mask(new bool[domain0*domain1]);

  for (size_t x = 0; x < domain0; x++) {
    for (size_t y = 0; y < domain1; y++) {
      const size_t idx = x*domain1 + y;
      if (bernoulli_distribution(p0)(r)) {
        mask[idx] = false;
        data[idx] = bernoulli_distribution(p1)(r);
      } else {
        mask[idx] = true;
      }
    }
  }

  return make_pair(move(data), move(mask));
}

static void
test3()
{
  random_device rd;
  rng_t r(rd());

  const vector<size_t> domains({10, 5});
  const vector<float> domain_alphas({2.0, 20.0});

  const vector<state::relation_t> relations({
      state::relation_t({0, 1}, distributions_factory<BetaBernoulli>().new_instance()),
      state::relation_t({1, 1}, distributions_factory<BetaBernoulli>().new_instance())
  });

  const vector<hyperparam_bag_t> relation_hps({
      beta_bernoulli_hp(2.0, 2.0),
      beta_bernoulli_hp(2.0, 3.0)
  });

  shared_ptr<state> s = make_shared<state>(domains, relations);

  for (size_t i = 0; i < domain_alphas.size(); i++) {
    state::message_type hp;
    hp.set_alpha(domain_alphas[i]);
    s->set_domain_hp(i, protobuf_to_string(hp));
  }

  for (size_t i = 0; i < relation_hps.size(); i++)
    s->set_relation_hp(i, relation_hps[i]);

  auto rel0 = binary_relation_generate(
      domains[0], domains[1],
      0.2, 0.8, r);

  auto rel1 = binary_relation_generate(
      domains[1], domains[1],
      0.3, 0.7, r);

  shared_ptr<dataview> rel0view(
    new row_major_dense_dataview(
        reinterpret_cast<uint8_t*>(rel0.first.get()),
        rel0.second.get(),
        {domains[0], domains[1]},
        runtime_type(TYPE_B)));

  shared_ptr<dataview> rel1view(
    new row_major_dense_dataview(
        reinterpret_cast<uint8_t*>(rel1.first.get()),
        rel1.second.get(),
        {domains[1], domains[1]},
        runtime_type(TYPE_B)));

  const vector<const dataview *> d({rel0view.get(), rel1view.get()});
  s->random_initialize(d, r);

  bound_state s0(s, 0, {rel0view, rel1view});
  bound_state s1(s, 1, {rel0view, rel1view});

  for (size_t i = 0; i < s0.nentities(); i++) {
    s0.remove_value(i, r);
    auto scores = s0.score_value(i, r);
    const auto choice = scores.first[util::sample_discrete_log(scores.second, r)];
    s0.add_value(choice, i, r);
  }

  for (size_t i = 0; i < s1.nentities(); i++) {
    s1.remove_value(i, r);
    auto scores = s1.score_value(i, r);
    const auto choice = scores.first[util::sample_discrete_log(scores.second, r)];
    s1.add_value(choice, i, r);
  }
}

int
main(void)
{
  test1();
  test2();
  test3();
  return 0;
}
