/*
 * bandit.h
 *
 *  Created on: Dec 3, 2011
 *      Author: baj
 */

#ifndef BANDIT_H_
#define BANDIT_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>

#include "utils.h"
#include "random.h"

class Bandit {
public:
	Bandit(int arms):
		mArms(arms),
		mOptimalArm(0),
		mMaxExpextedReward(0),
		mLow(0),
		mHigh(0)
	{

	}

	Bandit(int arms, double low, double high):
		mArms(arms),
		mLow(low),
		mHigh(high)
	{
		mExpectedReward.resize(mArms);
		std::generate(mExpectedReward.begin(), mExpectedReward.end(), UniformGenerator(mLow, mHigh));

		update();

		std::cout << "#Arms: " << mArms << std::endl;
		std::cout << "#MaxQ: " << mMaxExpextedReward << std::endl;
		std::cout << "#Optimal: " << mOptimalArm << std::endl;

		std::cout << std::endl;
	}

	virtual ~Bandit() {

	}

	void update() {
		mMaxExpextedReward = -1.0;

		for (int i = 0; i < mArms; ++i) {
			std::cout << "#Arm " << i << ": " << expected_reward(i) << std::endl;
			double eva = mExpectedReward[i];

			if (eva > mMaxExpextedReward) {
				mMaxExpextedReward = eva;
				mOptimalArm = i;
			}
		}
	}

	virtual std::pair<double, int> reward(int arm) const = 0;

	int arms() const {
		return mArms;
	}

	int optimal_arm() {
		return mOptimalArm;
	}

	double expected_reward(int arm) const {
		return mExpectedReward[arm];
	}

	double max_expected_reward() const {
		return mMaxExpextedReward;
	}

private:
	friend class MixtureBandit;

	int mArms;
	int mOptimalArm;
	double mMaxExpextedReward;

	std::vector<double> mExpectedReward;
	const double mLow;
	const double mHigh;
};


typedef boost::shared_ptr<Bandit> bandit_ptr;

bandit_ptr CreateBandit(const Params &params, const int arms);

class NormalBandit: public Bandit {
public:
	NormalBandit(int arms): Bandit(arms, -10.0, 10.0) {
		mSigma.resize(arms);
		std::generate(mSigma.begin(), mSigma.end(), UniformGenerator(1.0, 10.0));
	}

	virtual ~NormalBandit() {

	}

	std::pair<double, int> reward(int arm) const {
		return std::make_pair(boost::normal_distribution<>(expected_reward(arm), mSigma[arm])(RNG), 0);
	}

private:
	std::vector<double> mSigma;
};

class UniformBandit: public Bandit {
public:
	UniformBandit(int arms): Bandit(arms, -10.0, 10.0) {
	}

	virtual ~UniformBandit() {

	}

	std::pair<double, int> reward(int arm) const {
		return std::make_pair(UniformGenerator(expected_reward(arm) - 10.0, expected_reward(arm) + 10.0)(), 0);
	}
};

class BernoulliBandit: public Bandit {
public:
	BernoulliBandit(int arms): Bandit(arms, 0.0, 1.0) {

	}

	virtual ~BernoulliBandit() {

	}

	std::pair<double, int> reward(int arm) const {
		return std::make_pair(UniformGenerator(0.0, 1.0)() < expected_reward(arm)? 1.0: 0.0, 0);
	}
};

class MixtureBandit: public Bandit {
public:
	MixtureBandit(int arms, int num_basis): Bandit(arms) {
		mWeights.resize(arms);
		for (int i = 0; i < arms; ++i) {
	        std::vector<double> tmp;
	        for (int j = 0; j < num_basis - 1; ++j) {
	            tmp.push_back(drand48());
	        }
	        tmp.push_back(0.0);
	        tmp.push_back(1.0);

	        std::sort(tmp.begin(), tmp.end());

	        std::cout << "#Weights " << i << ": ";
	        mWeights[i].resize(num_basis);
	        for (int j = 0; j < num_basis; ++j) {
	        	mWeights[i][j] = tmp[j+1] - tmp[j];

	        	std::cout << mWeights[i][j];
	        	if (j != num_basis - 1) {
	        		std::cout << ", ";
	        	}
	        	else {
	        		std::cout << std::endl;
	        	}
	        }
		}

        std::cout << "#Basis: " << num_basis << std::endl;
		for (int i = 0; i < num_basis; ++i) {
			mBandits.push_back(bandit_ptr(new NormalBandit(arms)));
		}

		mExpectedReward.resize(mArms, 0.0);
		for (int i = 0; i < mArms; ++i) {
			for (int j = 0; j < num_basis; ++j) {
				mExpectedReward[i] += mBandits[j]->expected_reward(i) * mWeights[i][j];
			}
		}

		update();

		std::cout << "#Arms: " << mArms << std::endl;
		std::cout << "#MaxQ: " << mMaxExpextedReward << std::endl;
		std::cout << "#Optimal: " << mOptimalArm << std::endl;

		std::cout << std::endl;
	}

	virtual ~MixtureBandit() {
		for (int i = 0; i < arms(); ++i) {
			std::stringstream ss;
			ss << "MixtureBandit-" << i;
			DumpDistribution(MixtureGenerator(mBandits, mWeights[i], i), 1 << 22, ss.str().c_str(), false);
		}
	}

	std::pair<double, int> reward(int arm) const {
		return MixtureGeneratorImp(mBandits, mWeights[arm], arm)();
	}

	class MixtureGeneratorImp {
	public:
		MixtureGeneratorImp(const std::vector<bandit_ptr >& bandits, const std::vector<double>& weights, int arm):
			mBandits(bandits),
			mWeights(weights),
			mArm(arm)
		{

		}

		std::pair<double, int> operator()() {
			double prob = drand48();

			for (uint i = 0; i < mBandits.size(); ++i) {
				if (prob < mWeights[i]) {
					return std::make_pair(mBandits[i]->reward(mArm).first, i);
				}

				prob -= mWeights[i];
			}

			assert(0);
			return std::make_pair(0.0, 0);
		}

	private:
		const std::vector<bandit_ptr>& mBandits;
		const std::vector<double>& mWeights;
		const int mArm;
	};

	class MixtureGenerator {
		public:
			MixtureGenerator(const std::vector<bandit_ptr >& bandits, const std::vector<double>& weights, int arm):
				mBandits(bandits),
				mWeights(weights),
				mArm(arm)
			{

			}

			double operator()() {
				return MixtureGeneratorImp(mBandits, mWeights, mArm)().first;
			}

		private:
			const std::vector<bandit_ptr>& mBandits;
			const std::vector<double>& mWeights;
			const int mArm;
		};

private:
	std::vector<std::vector<double> > mWeights;
	std::vector<bandit_ptr> mBandits;
};

#endif /* BANDIT_H_ */
