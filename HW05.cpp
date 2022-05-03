#include <iostream>
#include <cmath>

//using fptr = double (double); // puntero de funcion con forma : double function(double x);

const double M = 120.0;
const double G = 9.81;
const double VT = 40;
const double T = 10;

double f(double x);

template <typename func_t>
double newton(double x0, double eps, func_t func, int nmax, int & nsteps);

int main (int argc, char *argv[])
{
  double X0 = std::atof(argv[1]);//10
  double eps = std::atof(argv[2]);//0.0001
int NMAX = 1000;
std::cout.precision(15); std::cout.setf(std::ios::scientific);
int steps = 0;
double root = newton(X0, eps, f, NMAX, steps);
 std::cout << "\t" << root << "\t" << f(root) << "\t" << steps << "\n";
 
}

double f(double x)
{
return M*G*(1 - std::exp(-x*T/M))/x - VT;
}

template <typename func_t>
double newton(double x0, double eps, func_t func, int nmax, int & nsteps)
{
nsteps = 0;
double xr = x0;
while(nsteps <= nmax) {
if (std::fabs(func(xr)) < eps) {
break;
} else {
double h = 0.001;
double deriv = (func(xr+h/2) - func(xr-h/2))/h;
xr = xr - func(xr)/deriv;
}
nsteps++;
}

return xr;
}
