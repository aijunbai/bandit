/*
 * random.cpp
 *
 *  Created on: Dec 3, 2011
 *      Author: baj
 */

#include "random.h"
#include "statistic.h"

boost::mt19937 RNG;

void RandomSeeding(int seed)
{
	srand48(seed);
	srand(seed);

	RNG.seed(seed);

	std::cout << "#Seed: " << seed << std::endl;
}
