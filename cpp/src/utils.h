/*
 * utils.h
 *
 *  Created on: May 17, 2013
 *      Author: baj
 */

#ifndef UTILS_H_
#define UTILS_H_


enum AgentType {
	QLearning = 1,
	UCB1 = 2,
	NormalGamma = 3,
	DirichletNormalGamma = 4
};

enum BanditType {
	Bernoulli = 1,
	Uniform = 2,
	Normal = 3,
	Mixture = 4
};

struct Params {
	Params(int seed):
		mBanditType(0),
		mAgentType(0),
		mConstant(-1.0),
		mEplison(0.1),
		mBeta(1000.0),
		mNumBasis(2),
		mSeed(seed),
		mArmsPower(0),
		mTrialsPower(0),
		mRunsPower(0)
	{

	}

	int mBanditType;
	int mAgentType;
	double mConstant;
	double mEplison;
	double mBeta;
	int mNumBasis;
	int mSeed;
	int mArmsPower;
	int mTrialsPower;
	int mRunsPower;
};

#endif /* UTILS_H_ */
