#include <iostream>
#include <algorithm>

#include "random.h"
#include "agent.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {
	Params params(getpid());

	try {
		po::options_description desc("Allowed options");
		desc.add_options()
				("help", "produce help message")
				("seed", po::value<int>(&params.mSeed), "random seed")
				("bandit", po::value<int>(&params.mBanditType), "bandit type")
				("agent", po::value<int>(&params.mAgentType), "agent type")
				("arms", po::value<int>(&params.mArmsPower), "arms (power of 2)")
				("trials", po::value<int>(&params.mTrialsPower), "trials (power of 2)")
				("runs", po::value<int>(&params.mRunsPower), "runs (power of 2)")
				("eplison", po::value<double>(&params.mEplison), "eplison for Q agent")
				("exploration", po::value<double>(&params.mConstant), "exploration constant for UCB agent")
				("beta", po::value<double>(&params.mBeta), "beta intializer for NormalGamma agent")
				("basis", po::value<int>(&params.mNumBasis), "number of basis of MixtureBandit")
				;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (params.mBanditType == 0 || params.mAgentType == 0 || vm.count("help")) {
			cout << desc << "\n";
			return 1;
		}
	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}
	catch(...) {
		cerr << "Exception of unknown type!\n";
		return 1;
	}

	RandomSeeding(params.mSeed);
	NormalGammaInfo::SetBETA(params.mBeta);

	const int trials = pow(2, params.mTrialsPower),
			arms = pow(2, params.mArmsPower),
			runs = pow(2, params.mRunsPower);

	bandit_ptr bandit = CreateBandit(params, arms);
	vector<STATISTIC> regret(runs);

	for (int i = 0; i < trials; ++i) {
		agent_ptr agent = CreateAgent(params, arms);

		double sum = 0.0;
		for (int j = 1; j < runs; ++j) {
			sum += agent->run(*bandit).first;
			regret[j].Add(j * bandit->max_expected_reward() - sum);
		}
	}

	for (int i = 1; i < runs; ++i) {
		cout << i << " " << regret[i].GetValue() << endl;
	}

	return 0;
}
