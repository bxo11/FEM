#include <cstdlib>
using namespace std;

double f(double x, double y);
double default_3points_integral_2D();
double default_2points_integral_2D();

double f(double x, double y) {
	return -5 * x * x * y + 2 * x * y * y + 10;
	//return 2 * x * x * y + 2 * x * y + 4;
}

double default_3points_integral_2D()
{
	double wynik = 0.;
	double w1 = 5. / 9.;
	double w2 = 8. / 9.;
	double x = sqrt(3. / 5.);
	double y = 0.;

	wynik += w1 * w1 * f(x, x);
	wynik += w1 * w1 * f(-x, -x);
	wynik += w1 * w1 * f(-x, x);
	wynik += w1 * w1 * f(x, -x);

	wynik += w1 * w2 * f(x, y);
	wynik += w1 * w2 * f(-x, -y);
	wynik += w1 * w2 * f(-x, y);
	wynik += w1 * w2 * f(x, -y);

	wynik += w2 * w2 * f(y, y);

	return wynik;
}

double default_2points_integral_2D()
{
	double wynik = 0.;
	double w1 = 1.;
	double w2 = 1.;
	double x = 1. / sqrt(3.);

	wynik += w1 * w1 * f(x, x);
	wynik += w1 * w2 * f(-x, -x);
	wynik += w2 * w1 * f(-x, x);
	wynik += w2 * w2 * f(x, -x);

	return wynik;
}
