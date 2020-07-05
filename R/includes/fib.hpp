//#include <boost/math/distributions/gamma.hpp> 

int fib(const int&n, std::ostream* pstream__) {
  if (n <= 0) {
    stringstream errmsg;
    errmsg << "n must be positive";
    throw std::domain_error(errmsg.str());
  }

//   boost::math::gamma_distribution<> dist(0.5, 5); 
//   return quantile(dist, 0.4);
  return n <= 2 ? 1 : fib(n - 1, 0) + fib(n - 2, 0);
}
