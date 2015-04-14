
#ifndef RANDOMCLASS
#define RANDOMCLASS
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random.hpp>
#include <iostream>
using namespace std;


typedef boost::mt19937 base_generator_type;
//typedef boost::ecuyer1988 base_generator_type;
//typedef boost::minstd_rand base_generator_type;
class RandomClass
{
private:
    base_generator_type m_generator;
    boost::random_number_generator<base_generator_type > m_stl_rand;
    boost::uniform_real<> m_uni_dist;
    boost::normal_distribution<> m_normal_dist;
    boost::variate_generator<base_generator_type&, boost::uniform_real<> > m_randu;
    boost::variate_generator<base_generator_type&, boost::normal_distribution<> > m_randn;
    unsigned long timesCalled;
public:
    RandomClass();
    void initialize(unsigned int seed);
    double randn();
    double randu();
};




inline RandomClass::RandomClass(): m_generator(0), m_stl_rand(m_generator), m_uni_dist(0,1),
				   m_normal_dist(0,1), m_randu(m_generator, m_uni_dist), m_randn(m_generator, m_normal_dist)
{
    timesCalled=0;
    return;
}

inline void RandomClass::initialize(unsigned int seed)
{
    m_generator.seed(seed);
}
inline double RandomClass::randn()
{
    return m_randn();
}
inline double RandomClass::randu()
{
    timesCalled++;
    double ret = m_randu();
    m_randu();
    return ret;;
}
#endif
