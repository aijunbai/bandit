/*
 * random.h
 *
 *  Created on: Dec 3, 2011
 *      Author: baj
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions.hpp>

#include "statistic.h"

void RandomSeeding(int seed);

#endif /* RANDOM_H_ */
