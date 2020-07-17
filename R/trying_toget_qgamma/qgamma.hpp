inline double qgamma(double p, double shape, double scale, std::ostream* pstream) {
   #include <boost/math/distributions/gamma.hpp> 

  boost::math::gamma_distribution<> dist(shape, scale); 
  double q;
  q = quantile(dist, 0.4);


  return q;
}
