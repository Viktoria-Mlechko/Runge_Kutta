#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

double function1(double y, double x, double z)
{
	return(-2.4 * x * y * y + z * z - x - 1);
}
double function2(double y, double x, double z)
{
	return(1 / (1.27 * 1.27) - y - (x / z));
}

void runge_kutta_2(double x0, double z0, double yn, double zn, double y0, double h, double xn, int n)
{		
	ofstream file_to_write;


	double k[2], m[2];
	int i;
	for (i = 0; i < n; i++)
	{

		k[0] = h * (function1(x0, y0, z0));
		m[0] = h * (function2(x0, y0, z0));

		k[1] = h * (function1((x0 + h), (y0 + k[0]), (z0 + m[0])));
		m[1] = h * (function2((x0 + h), (y0 + k[0]), (z0 + m[0])));

		yn = y0 + (k[0] + k[1]) / 2;
		zn = z0 + (m[0] + m[1]) / 2;

		file_to_write.open("y runge-kutta 2.txt", ios::app);
		if (file_to_write.is_open())
		{
			file_to_write << x0 << " " << y0 << endl;
		}
		file_to_write.close();


		file_to_write.open("z runge-kutta 2.txt", ios::app);
		if (file_to_write.is_open())
		{
			file_to_write << x0 << " " << z0 << endl;
		}
		file_to_write.close();
		cout << x0 << "\t" << y0 << "\t" << yn << endl;
		cout << x0 << "\t" << z0 << "\t" << zn << endl << endl;
		x0 = x0 + h;
		y0 = yn;
		z0 = zn;
	}
	cout << "\nValue of y at x = " << xn << " is " << yn << endl << "\nValue of z at x = " << xn << " is " << zn << endl;

}
// y(0)=1/1.2; z(0)=1; [0;1]

void runge_kutta_4(double x0, double z0, double yn, double zn, double y0, double h, double xn, int n)
{
	ofstream file_to_write;

	double k[4], m[4];
	int i;
	for (i = 0; i < n; i++)
	{
		k[0] = h * (function1(x0, y0, z0));
		m[0] = h * (function2(x0, y0, z0));

		k[1] = h * (function1((x0 + h / 2), (y0 + k[0] / 2), (z0 + m[0] / 2)));
		m[1] = h * (function2((x0 + h / 2), (y0 + k[0] / 2), (z0 + m[0] / 2)));

		k[2] = h * (function1((x0 + h / 2), (y0 + k[1] / 2), (z0 + m[1] / 2)));
		m[2] = h * (function2((x0 + h / 2), (y0 + k[1] / 2), (z0 + m[1] / 2)));

		k[3] = h * (function1((x0 + h), (y0 + k[2]), (z0 + m[2])));
		m[3] = h * (function2((x0 + h), (y0 + k[2]), (z0 + m[2])));


	file_to_write.open("y runge-kutta 4.txt", ios::app);
	if (file_to_write.is_open())
	{
		file_to_write << x0 << " " << y0 << endl;
	}
	file_to_write.close();


	file_to_write.open("z runge-kutta 4.txt", ios::app);
	if (file_to_write.is_open())
	{
		file_to_write << x0 << " " << z0 << endl;
	}
	file_to_write.close();

		yn = y0 + (k[0] + 2 * k[1] + 2 * k[2] + k[3]) / 6;
		zn = z0 + (m[0] + 2 * m[1] + 2 * m[2] + m[3]) / 6;

		cout << x0 << "\t" << y0 << "\t" << yn << endl;
		cout << x0 << "\t" << z0 << "\t" << zn << endl << endl;
		x0 = x0 + h;
		y0 = yn;
		z0 = zn;
	}
	cout << "\nValue of y at x = " << xn << " is " << yn << endl << "\nValue of z at x = " << xn << " is " << zn << endl;

}
int main()
{
	double x0, z0, y0, xn, zn{}, h, yn{}, k[4], m[4];
	int i, n;

	cout << "Enter Initial Condition" << endl;
	cout << "y0 = ";
	cin >> y0;
	cout << "z0 = ";
	cin >> z0;
	cout << "x0 = ";
	cin >> x0;
	cout << "Enter calculation point xn = ";
	cin >> xn;
	cout << "Enter number of steps: ";
	cin >> n;
	

	h = (xn - x0) / n;
	cout << "\nx0\ty0\tyn\n\n\n";
	cout << "\n\tUsing Runge-Kutta second order method\n\n";
	runge_kutta_2(x0, z0, yn, zn, y0, h, xn, n);
	
	cout << "\n\tUsing Runge-Kutta forth order method\n\n";
	runge_kutta_4(x0, z0, yn, zn, y0, h, xn, n);
	return 0;
}