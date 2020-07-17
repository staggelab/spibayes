functions {
vector qgamma(vector y,        // unknowns
              vector theta,    // parameters
              real[] x_r,      // data (real)
              int[] x_i) {     // data (integer)
  vector[1] z;
  z[1] = gamma_cdf(y, theta[1], 1/theta[2]) - theta[3];
  return z;
}
}
data {
int<lower=0> N;
vector[N] y;
}
parameters {
real<lower=0> shape;
real<lower=0> scale;
}
transformed parameters {
real<lower = 0> p_goal;

p_goal = shape /2;

}
model {
shape ~ gamma(0.5,0.2);
scale ~ gamma(0.5,0.2);

for (n in 1:N)
  y[n] ~ gamma(shape, 1/scale);
}
generated quantities {
vector[1] qgamma_result;
//vector[4] theta = [shape, scale, y[1], p_val]';
vector[3] theta ;
//vector[2] theta;
vector[1] y_guess = [0.5]';

int x_i[0];
real x_r[0]; // Return the 40th percentile

theta = [shape, scale, p_goal]';
//theta = [shape, scale]' ;

qgamma_result = algebra_solver(qgamma, y_guess, theta, x_r, x_i);

//qgamma_result = algebra_solver(qgamma, [1.3]', [shape, scale, p_goal]', x_r, x_i);

}
