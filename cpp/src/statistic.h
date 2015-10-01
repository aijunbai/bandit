#ifndef STATISTIC_H
#define STATISTIC_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/math/distributions.hpp>

#define MULTI_ROLLOUTS 0
#define DIRICHLET_INIT 1
#define IMMEDIATE_BETA 1

class GammaGenerator;

extern boost::mt19937 RNG;
extern GammaGenerator GAMMA;

template <typename Generator>
inline void DumpDistribution(Generator gen, int samples, const char *file_name, bool append) {
    const int scale = 100;

    std::map<int, int> samples_count;

    for (int i = 0; i < samples; ++i) {
        samples_count[rint(gen() * scale)] += 1;
    }

    std::ofstream fout;
    if (append) {
    	fout.open(file_name, std::ios_base::out | std::ios_base::app);
    }
    else {
    	fout.open(file_name, std::ios_base::out | std::ios_base::trunc);
    }

    if (fout.good()) {
        for (std::map<int, int>::iterator it = samples_count.begin(); it != samples_count.end(); ++it) {
        	double p = it->first / double(scale);

        	if (p > -50 && p < 50) {
        		fout << p << " " << it->second << std::endl;
        	}
        }

        fout.close();
    }
}

template<typename _Tp>
inline const _Tp&
Max(const _Tp& x, const _Tp& y)
{
    return std::max(x, y);
}

template<typename _Tp>
inline const _Tp&
Min(const _Tp& x, const _Tp& y)
{
    return std::min(x, y);
}

template<typename _Tp>
inline const _Tp&
MinMax(const _Tp& min, const _Tp& x, const _Tp& max)
{
    return Min(Max(min, x), max);
}


template<class COUNT>
class VALUE
{
public:
	VALUE() {
		Count = 0;
		Total = 0;
	}

    void Set(double count, double value)
    {
        Count = count;
        Total = value * count;
    }

    void Add(double x)
    {
        Count += 1.0;
        Total += x;
    }

    void Add(double x, COUNT weight)
    {
        Count += weight;
        Total += x * weight;
    }

    double GetValue() const
    {
        return Count == 0 ? Total : Total / double(Count);
    }

    COUNT GetCount() const
    {
        return Count;
    }

private:

    COUNT Count;
    double Total;
};


class STATISTIC
{
public:

    STATISTIC();

    ~STATISTIC() {

    }

    int GetCount() const;
    double GetValue() const; //merge from VALUE
    void Set(int count, double value); //merge from VALUE

    void Initialise();
    double GetTotal() const;
    double GetMean() const;
    double GetVariance() const;
    double GetStdDev() const; //标准差
    double GetStdErr() const; //标准误差（平均值的标准差）
    double GetMax() const;
    double GetMin() const;

    void SetMin(double min) {
    	Min = min;
    }

    void SetMax(double max) {
    	Max = max;
    }

    void AdjustRange(double min, double max) {
    	Min = std::min(min, Min);
    	Max = std::max(max, Max);
    }

    void Print(const std::string& name, std::ostream& ostr) const;
    void Add(double val);

private:
    double Count;
    double Mean; //平均值
    double Variance; //方差
    double Min, Max;
};

inline STATISTIC::STATISTIC()
{
	Initialise();
}

inline void STATISTIC::Set(int count, double value)
{
	Initialise();

	Count = count;
	Mean = value;
}

inline void STATISTIC::Add(double val)
{
    double meanOld = Mean;
    int countOld = Count;

    ++Count;
    assert(Count > 0); // overflow
    Mean += (val - Mean) / Count;
    Variance = (countOld * (Variance + meanOld * meanOld) + val * val) / Count - Mean * Mean;

    if (Variance < 0.0)
    	Variance = 0.0;
    if (val > Max)
        Max = val;
    if (val < Min)
        Min = val;
}

inline void STATISTIC::Initialise()
{
    Count = 0;
    Mean = 0;
    Variance = 0;
    Min = +10e6;
    Max = -10e6;
}

inline int STATISTIC::GetCount() const
{
    return Count;
}

inline double STATISTIC::GetTotal() const
{
    return Mean * Count;
}

inline double STATISTIC::GetValue() const
{
	return GetMean();
}

inline double STATISTIC::GetMean() const
{
    return Mean;
}

inline double STATISTIC::GetVariance() const
{
	return Variance;
}

inline double STATISTIC::GetStdDev() const
{
    return sqrt(Variance);
}

inline double STATISTIC::GetStdErr() const
{
    return sqrt(Variance / Count);
}

inline double STATISTIC::GetMax() const
{
    return Max;
}

inline double STATISTIC::GetMin() const
{
    return Min;
}
    
inline void STATISTIC::Print(const std::string& name, std::ostream& ostr) const
{
    ostr << name
    		<< ": " << Mean
    		<< " (" << GetCount()
    		<< ") [" << Min
    		<< ", " << Max
    		<< "] +- " << GetStdErr()
    		<< ", sigma=" << GetStdDev()
    		<< std::endl;
}

class UniformGenerator {
public:
	UniformGenerator(double low, double high): mLow(low), mHigh(high) {

	}

	double operator()() {
		return mLow + drand48() * (mHigh - mLow);
	}

	double mLow;
	double mHigh;
};

class NormalGammaGenerator {
public:
	NormalGammaGenerator(double mu, double lambda, double alpha, double beta): Mu(mu), Lambda(lambda), Alpha(alpha), Beta(beta) {

	}

	double operator()() {
		const double t = boost::gamma_distribution<>(Alpha, 1.0 / Beta)(RNG);
		const double p = std::max(Lambda * t, 1.0e-6);
		const double m = boost::normal_distribution<>(Mu, sqrt(1.0 / p))(RNG);

		return m;
	}

private:
	const double Mu;
	const double Lambda;
	const double Alpha;
	const double Beta;
};

class BetaInfo {
public:
	BetaInfo() {
		Initialise();
	}

    ~BetaInfo() {

    }

    void Initialise() {
    	Alpha = ALPHA;
    	Beta = BETA;

    	Mean.Set(0, 0.0);
    }

    void Set(int count, double value) {
    	Alpha = ALPHA;
    	Beta = BETA;

    	Mean.Set(count, value);
    }

	void Add(double value) { //add a new sample
#if IMMEDIATE_BETA
		double prob = (value - MIN) / (MAX - MIN);
		double trial = drand48(); //boost::uniform_real<>(0.0, 1.0)(RNG);

		if (trial < prob) { //sucess
			Alpha += 1.0;
		}
		else {
			Beta += 1.0;
		}
#endif

		Mean.Add(value);
	}

	double GetExpectation() const {
		return ThompsonSampling(false);
	}

	double ThompsonSampling(bool sampling = true) const { //Two Step: 采样一个模型参数，并计算出该模型参数对应的期望收益
#if IMMEDIATE_BETA
		if (sampling) {
			double x = boost::gamma_distribution<>(Alpha)(RNG);
			double y = boost::gamma_distribution<>(Beta)(RNG);

			return x / (x + y) * (MAX - MIN) + MIN;
		}
		else {
			return Mean.GetValue();
		}
#else
		(void) sampling;
		return Mean.GetValue();
#endif
	}

	void Print(const std::string& name, std::ostream& ostr) const
	{
	    ostr << name
	    		<< "{"
#if IMMEDIATE_BETA
	    		<< "a=" << Alpha
	    		<< " b=" << Beta
	    		<< " p=" << Alpha / (Alpha + Beta) << " "
#endif
	    		<< "m=" << GetExpectation()
	    		<< "}";
	}

	static void setMinMax(double min, double max) {
		MIN = min;
		MAX = max;
	}

private:
	static double MIN;
	static double MAX;
	static double ALPHA;
	static double BETA;

	double Alpha;
	double Beta;

	VALUE<double> Mean;
};

class NormalGammaInfo {
public:
	NormalGammaInfo():
		Mu(0.0),
		Lambda(0.0),
		Alpha(ALPHA),
		Beta(BETA) {
		Initialise();
	}

	~NormalGammaInfo() {

	}

	void Initialise() {
		Mu = 0.0;
		Lambda = 0.0;
		Alpha = ALPHA;
		Beta = BETA;
	}

	double GetValue() const {
		return Mu;
	}

	double GetCount() const {
		return Lambda;
	}

	double GetAlpha() const {
		return Alpha;
	}

	double GetBeta() const {
		return Beta;
	}

	double GetExpectation() const {
		return ThompsonSampling(false);
	}

	void Set(int count, double value) {
		Mu = value;
		Lambda = count;
		Alpha = ALPHA;
		Beta = BETA;
	}

	void Add(const std::vector<double>& values) { //add a new sample
		STATISTIC samples;
		for (uint i = 0; i < values.size(); ++i) {
			samples.Add(values[i]);
		}

		double n = samples.GetCount();
		double m = samples.GetMean();
		double s = samples.GetVariance();

		double mu = (Lambda * Mu + n * m) / (Lambda + n);
		double lambda = Lambda + n;
		double alpha = Alpha + 0.5 * n;
		double beta = Beta + 0.5 * n * s + 0.5 * (Lambda * n / (Lambda + n)) * (m - Mu) * (m - Mu) ;

		Mu = mu;
		Lambda = lambda;
		Alpha = alpha;
		Beta = beta;
	}

	double ThompsonSampling(bool sampling = true) const { //Two Step: 采样一个模型参数，并计算出该模型参数对应的期望收益
		return sampling? NormalGammaGenerator(Mu, Lambda, Alpha, Beta)(): Mu;
	}

	void Print(const std::string& name, std::ostream& ostr) const
	{
		ostr << name << ":"
			<< " mu=" << Mu
			<< " lambda=" << Lambda
			<< " alpha=" << Alpha
			<< " beta=" << Beta
			<< " error=" << sqrt(Beta / (Lambda * (Alpha - 1)))
			<< " sigma=" << sqrt(Beta / (Alpha - 1))
			<< std::endl;
	}

	static void SetALPHA(double alpha) {
		ALPHA = alpha;
	}

	static void SetBETA(double beta) {
		BETA = beta;
	}

private:
	double Mu;
	double Lambda;
	double Alpha;
	double Beta;

	static double ALPHA;
	static double BETA;
};

class DirichletNormalGammaInfo {
	typedef boost::unordered_map<int, double> map_a;
	typedef boost::unordered_map<int, NormalGammaInfo> map_ng;

public:
	DirichletNormalGammaInfo() {
		Initialise();
	}

	~DirichletNormalGammaInfo() {

	}

	void Initialise() {
//		Alpha.resize(Order);
//		NormalGamma.resize(Order);
//
//		for (int i = 0; i < Order; ++i) {
//			Alpha[i] = 1;
//			NormalGamma[i].Initialise();
//		}
	}

	double GetExpectation() const {
		double ret = 0.0;
		double sum = 0.0;

		for (map_a::const_iterator it = Alpha.begin(); it != Alpha.end(); ++it) {
			ret += it->second * NormalGamma[it->first].GetExpectation();
			sum += it->second;
		}

		return ret / sum;
	}

//	void Set(int count, double value) {
//		Alpha.resize(Order);
//
//		for (int i = 0; i < Order; ++i) {
//		    Alpha[i] = count;
//			NormalGamma[i].Set(count, value);
//		}
//	}

	void Initialise(int k) {
		for (int i = 0; i < k; ++i) {
			Alpha[k] = 1;
			NormalGamma[k].Initialise();
		}
	}

	void Add(const std::vector<double>& values, int k) { //add a new sample
		Alpha[k] += 1;
		NormalGamma[k].Add(values);
	}

	double ThompsonSampling() const { //Two Step: 采样一个模型参数，并计算出该模型参数对应的期望收益
		map_a p;

		double sum = 0.0;
		for (map_a::const_iterator it = Alpha.begin(); it != Alpha.end(); ++it) {
			p[it->first] = boost::gamma_distribution<>(it->second)(RNG);
			sum += p[it->first];
		}

		double ret = 0.0;
		for (map_ng::const_iterator it = NormalGamma.begin(); it != NormalGamma.end(); ++it) {
			p[it->first] /= sum;
			ret += p[it->first] * it->second.ThompsonSampling();
		}

		return ret;
	}

	void Print(const std::string& name, std::ostream& ostr) const
	{
		double sum = 0.0;
		for (map_a::const_iterator it = Alpha.begin(); it != Alpha.end(); ++it) {
			sum += it->second;
		}

		ostr << name << ": Alpha=(";
		for (map_a::const_iterator it = Alpha.begin(); it != Alpha.end(); ++it) {
			ostr << it->second / sum << ", ";
		}
		ostr << ")" << std::endl;

		for (map_ng::const_iterator it = NormalGamma.begin(); it != NormalGamma.end(); ++it) {
			std::stringstream ss;
			ss << "#NormalGamma " << it->first;
			it->second.Print(ss.str().c_str(), ostr);
		}
	}

private:
    mutable map_a Alpha;
    mutable map_ng NormalGamma;
};

template <class T, class H>
class DirichletInfo {
public:
    const DirichletInfo<T, H> &operator=(const DirichletInfo<T, H> &o) {
        Alpha = o.Alpha;

        return *this;
    }

	const std::vector<std::pair<T, double> > &GetExpectation() const {
		return ThompsonSampling(false);
	}

	void Initial(const T &s) {
		Alpha[s] = 0.5; //XXX
	}

	void Clear() {
		Alpha.clear();
	}

	void Add(const T &s) {
		Alpha[s] += 1.0;
	}

	const std::vector<std::pair<T, double> > &ThompsonSampling(bool sampling = true) const { //Two Step: 采样一个模型参数，并计算出该模型参数对应的期望收益
		outcomes_.clear();

		double sum = 0.0;
        for(typename H::iterator it = Alpha.begin(); it != Alpha.end(); ++it ) {
        	outcomes_.push_back(std::make_pair(it->first, 0));
        	outcomes_.back().second = sampling? boost::gamma_distribution<>(it->second)(RNG): it->second;
        	sum += outcomes_.back().second;
        }

        for (typename std::vector<std::pair<T, double> >::iterator it = outcomes_.begin(); it != outcomes_.end(); ++it) {
        	it->second /= sum;
        }

		return outcomes_;
	}

	void Print(const std::string& name, std::ostream& ostr) const
	{
		const std::vector<std::pair<T, double> > &outcomes = GetExpectation();

		ostr << name << ": Alpha=(";
		for (typename std::vector<std::pair<T, double> >::const_iterator it = outcomes.begin(); it != outcomes.end(); ++it) {
			ostr << "(" << it->second << "), ";
		}
		ostr << ")" << std::endl;
	}

private:
	mutable H Alpha;
	mutable std::vector<std::pair<T, double> > outcomes_;
};

class DirichletInfo_POMCP: public DirichletInfo<int, boost::unordered_map<int, double> > {
public:
	void Initialise(int k) {
		for (int i = 0; i < k; ++i) {
			Initial(i);
		}
	}

	void Print(const std::string& name, std::ostream& ostr) const
	{
		const std::vector<std::pair<int, double> > &outcomes = GetExpectation();

		ostr << name << "{";
		for (std::vector<std::pair<int, double> >::const_iterator it = outcomes.begin(); it != outcomes.end(); ++it) {
			ostr << "#" << it->first << "(" << it->second << "), ";
		}
		ostr << "}";
	}
};

#endif // STATISTIC
