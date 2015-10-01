/*
 * agent.h
 *
 *  Created on: Dec 3, 2011
 *      Author: baj
 */

#ifndef AGENT_H_
#define AGENT_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/math/distributions.hpp>

#include "utils.h"
#include "bandit.h"
#include "statistic.h"

class MetaAgent {
public:
	MetaAgent() {

	}

	virtual ~MetaAgent() {

	}

	virtual std::pair<double, int> run(const Bandit& bandit) = 0;
};


typedef boost::shared_ptr<MetaAgent> agent_ptr;

agent_ptr CreateAgent(const Params &params, const int arms);

template <typename Data>
class Agent {
public:
	Agent(int arms): mArms(arms), mT(0) {
		mN.resize(arms);
		mQ.resize(arms);

		fill(mN.begin(), mN.end(), 0);
	}

	virtual ~Agent() {
		for (uint i = 0; i < mQ.size(); ++i) {
			std::cout << "#Arm " << i << ": n=" << mN[i] << " ";
			mQ[i].Print("", std::cout);
		}
		std::cout << std::endl;
	}

	int arms() const {
		return mArms;
	}

	int& T() { return mT; }
	std::vector<int>& N() { return mN; }

	Data& Q(int arm) { return mQ[arm]; }
	const Data& Q(int arm) const { return mQ[arm]; }

	void UpdateCount(int arm) {
		N()[arm] += 1;
		T() += 1;
	}

private:
	int mArms;
	int mT;
	std::vector<int> mN;
	std::vector<Data> mQ;
};

class QAgent: public Agent<STATISTIC>, public MetaAgent {
public:
	QAgent(int arms, double epsilon = 0.1): Agent(arms), mEpsilon(epsilon) {
	}

	virtual ~QAgent() {
	}

	virtual std::pair<double, int> run(const Bandit & bandit) {
		int arm = 0;

		if (UniformGenerator(0.0, 1.0)() < mEpsilon * pow(1.0 - 1.0e-3, T())) {
			arm = random_selection();
		}
		else {
			arm = greedy_selection();
		}

		double reward = bandit.reward(arm).first;

		Q(arm).Add(reward);

		UpdateCount(arm);
		return std::make_pair(reward, arm);
	}

	int random_selection() {
		return rand() % arms();
	}

	int greedy_selection() {
		double max = -1.0e6;
		int best = 0;

		for (int i = 0; i < arms(); ++i) {
			if (Q(i).GetValue() > max) {
				max = Q(i).GetValue();
				best = i;
			}
		}

		return best;
	}

private:
	double mEpsilon;
};

class UCBAgent: public Agent<STATISTIC>, public MetaAgent {
public:
	UCBAgent(int arms, double c): Agent(arms), mC(c) {
	}

	virtual ~UCBAgent() {
	}

	virtual std::pair<double, int> run(const Bandit & bandit) {
		int arm = UCB1();
		double reward = bandit.reward(arm).first;

		Q(arm).Add(reward);

		UpdateCount(arm);
		return std::make_pair(reward, arm);
	}

	int UCB1() {
		double max = -1.0e6;
		int best = 0;

		for (int i = 0; i < arms(); ++i) {
			double c = mC;
			double upper_boud = 1.0e6;

			if (c < 0.0) {
				c = Q(i).GetValue();
			}

			if (N()[i] > 0) {
				double bonous = c * Confidence(T(), N()[i]);
				upper_boud = Q(i).GetValue() + bonous;
			}

			if (upper_boud > max) {
				max = upper_boud;
				best = i;
			}
		}

		return best;
	}

	double Confidence(int t, int n) {
		if (n > 0) {
			return sqrt(2.0 * log(t) / n);
		}
		else {
			return 1.0e6;
		}
	}

private:
	const double mC;
};

template <typename Data>
class ThompsonSamplingAgent: public Agent<Data>, public MetaAgent {
public:
	ThompsonSamplingAgent(int arms): Agent<Data>(arms) {

	}

	virtual ~ThompsonSamplingAgent() {

	}

	int thomas_sampling() {
		double max = -1.0e6;
		int best = 0;

		for (int i = 0; i < this->arms(); ++i) {
			if (this->N()[i] < 1) {
				return i;
			}

			double sampled_expected_reward = this->Q(i).ThompsonSampling();

			if (sampled_expected_reward > max) {
				max = sampled_expected_reward;
				best = i;
			}
		}

		return best;
	}
};

class NormalGammaAgent: public ThompsonSamplingAgent<NormalGammaInfo> {
public:
	NormalGammaAgent(int arms): ThompsonSamplingAgent(arms) {
	}

	virtual ~NormalGammaAgent() {
	}

	virtual std::pair<double, int> run(const Bandit & bandit) {
		int arm = thomas_sampling();

		double reward = bandit.reward(arm).first; //observed reward

		Q(arm).Add(std::vector<double>(1, reward));

		UpdateCount(arm);
		return std::make_pair(reward, arm);
	}
};

class DirichletNormalGammaAgent: public ThompsonSamplingAgent<DirichletNormalGammaInfo> {
public:
	DirichletNormalGammaAgent(int arms, int num_basis): ThompsonSamplingAgent(arms) {
		for (int i = 0; i < arms; ++i) {
            Q(i).Initialise(num_basis);
		}
	}

	virtual ~DirichletNormalGammaAgent() {

	}

	virtual std::pair<double, int> run(const Bandit & bandit) {
		int arm = thomas_sampling();

		std::pair<double, int> ret = bandit.reward(arm); //observed reward

		Q(arm).Add(std::vector<double>(1, ret.first), ret.second);

		UpdateCount(arm);
		return std::make_pair(ret.first, arm);
	}
};

#endif /* AGENT_H_ */
