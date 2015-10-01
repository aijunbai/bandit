/*
 * agent.cpp
 *
 *  Created on: Dec 3, 2011
 *      Author: baj
 */

#include "agent.h"

double NormalGammaInfo::ALPHA = 0.01;
double NormalGammaInfo::BETA = 100.0;

agent_ptr CreateAgent(const Params &params, const int arms)
{
	switch (params.mAgentType) {
	case QLearning: return agent_ptr(new QAgent(arms, params.mEplison));
	case UCB1: return agent_ptr(new UCBAgent(arms, params.mConstant));
	case NormalGamma: return agent_ptr(new NormalGammaAgent(arms));
	case DirichletNormalGamma: return agent_ptr(new DirichletNormalGammaAgent(arms, params.mNumBasis));
	default: return agent_ptr();
	}
}
