/*
 * bandit.cpp
 *
 *  Created on: May 17, 2013
 *      Author: baj
 */

#include "bandit.h"

bandit_ptr CreateBandit(const Params &params, const int arms)
{
	switch (params.mBanditType) {
	case Bernoulli: return bandit_ptr(new BernoulliBandit(arms));
	case Uniform: return bandit_ptr(new UniformBandit(arms));
	case Normal: return bandit_ptr(new NormalBandit(arms));
	case Mixture: return bandit_ptr(new MixtureBandit(arms, params.mNumBasis));
	default: return bandit_ptr();
	}
}


