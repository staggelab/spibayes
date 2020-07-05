
int ka(int&n, std::ostream* pstream__) {
  if (n <= 0) {
    stringstream errmsg;
    errmsg << "n must be positive";
    throw std::domain_error(errmsg.str());
  }

//   boost::math::gamma_distribution<> dist(0.5, 5); 
//   return quantile(dist, 0.4);
//   int q;
//   q = quantile(dist, 0.4);
//   q = 4 + 8;

//   return q;
  return n <= 2 ? 1 : fib(n - 1, 0) + fib(n - 2, 0);
}


double
sinc(const double& x, std::ostream* pstream__) {
  return x != 0.0 ? sin(x) / x : 1.0;
}

var
sinc(const var& x, std::ostream* pstream__) {
  double x_ = x.val();
  double f = x_ != 0.0 ? sin(x_) / x_ : 1.0;
  double dfdx_ = x_ != 0.0 ? (cos(x_) - sin(x_)) / x_ : 0.0;
  return var(new precomp_v_vari(f, x.vi_, dfdx_));
}


#include <vector>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>

double
qgamma(const double& p, double& shape, std::ostream* pstream__) {
    boost::math::gamma_distribution<> gam(0.5, 1);
	Q0 = boost::math::quantile(gam, p);
    std::vector<> Q;
    Q = boost::math::tools::real_cast<WP>(Q0));
  return Q0;
}



     boost::math::gamma_distribution<HP, my_policy> gam(alpha);
            HP Q0 = boost::math::quantile(gam, normcdf);
            std::vector<WP> Q;
            Q.push_back(boost::math::tools::real_cast<WP>(Q0));
