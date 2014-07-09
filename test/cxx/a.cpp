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
CheckDataview2DArray(
    const bool *data,
    const bool *mask,
    size_t n,
    size_t m,
    const dataview &d)
{
  // check each slice
  for (size_t i = 0; i < n; i++) {
    vector<bool> seen(m, false);
    for (auto &p : d.slice(0, i)) {
      MICROSCOPES_DCHECK(p.first[0] == i, "not a valid slice");
      MICROSCOPES_DCHECK(p.first[1] < m, "out of bounds");
      const size_t idx = i*m + p.first[1];
      const bool value = p.second.get<bool>(0);
      MICROSCOPES_DCHECK(data[idx] == value, "values don't match");
      MICROSCOPES_DCHECK(!mask[idx], "data is masked");
      seen[p.first[1]] = true;

      //cout << p.first << " : " << p.second.debug_str() << endl;
    }
    for (size_t j = 0; j < seen.size(); j++) {
      if (seen[j])
        continue;
      const size_t idx = i*m + j;
      //if (!mask[idx])
      //  cout << "should have seen (" << i << ", " << j << ")" << endl;
      MICROSCOPES_DCHECK(mask[idx], "data is not masked");
    }
  }
}

int
main(void)
{
  //random_device rd;
  //rng_t r(rd());

  rng_t r(34);

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

  // create the dataview
  //shared_ptr<dataview> view =
  //  make_shared<row_major_dense_dataview>(
  //      reinterpret_cast<uint8_t*>(likes), masks,
  //      domains, runtime_type(TYPE_B));

  const dataview *view =
    new row_major_dense_dataview(
        reinterpret_cast<uint8_t*>(likes), masks,
        domains, runtime_type(TYPE_B));

  for (auto &p : view->slice(0, 1)) {
    cout << p.first << " => " << p.second.debug_str() << endl;
  }

  CheckDataview2DArray(likes, masks, domains[0], domains[1], *view);


  // create groups for initial assignment; 2 user groups, 3 movie groups
  s.create_group(0);
  s.create_group(0);

  s.create_group(1);
  s.create_group(1);
  s.create_group(1);

  // bootstrap the movies first
  for (size_t m = 0; m < domains[1]; m++) {
    size_t gid;
    if (m < domains[1]/3)
      gid = 0;
    else if (m < (2*domains[1])/3)
      gid = 1;
    else
      gid = 2;
    s.assign_value(1, gid, m);
  }

  // now add the data
  for (size_t u = 0; u < domains[0]; u++) {
    size_t gid;
    if (u < domains[0]/2)
      gid = 0;
    else
      gid = 1;
    s.add_value(0, gid, u, {view}, r);
  }

  // peek @ suffstats
  for (const auto &p : s.get_suff_stats(0)) {
    cout << p.first << " : " << p.second.count_ << endl;
  }

  // score the 1st data point
  s.remove_value(0, 0, {view}, r);
  s.create_group(0);
  const auto scores = s.score_value(0, 0, {view}, r);
  cout << "scores: " << scores << endl;
  s.add_value(0, 0, 0, {view}, r);

  // remove the data
  for (size_t u = 0; u < domains[0]; u++) {
    s.remove_value(0, u, {view}, r);
  }

  // peek @ suffstats
  cout << "--" << endl;
  for (const auto &p : s.get_suff_stats(0)) {
    cout << p.first << " : " << p.second.count_ << endl;
  }

  delete [] likes;
  delete [] masks;
  delete view;
  return 0;
}
