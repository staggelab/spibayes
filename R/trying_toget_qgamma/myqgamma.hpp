#include <boost/math/distributions.hpp> 

double myqgamma(double p, double shape, double scale, std::ostream* pstream__) {
  boost::math::gamma_distribution<> dist(shape, scale); 
  return quantile(dist, p);
}
