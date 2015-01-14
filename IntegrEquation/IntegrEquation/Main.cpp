#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include "Matrix.h"
#include "Vector.h"
#include <ctime>    

using namespace std;

int a = 0;
int b = 1;
double eps = 0.0005;

double yTemp = 7.0 / 3;

Vect operator *(Matr A, Vect b)
{
	Vect c;
	/*c = InitVect(b.size);
	c = EnterZero(c);*/

	c = InitZeroVect(b.size);

	for (int i = 0; i < b.size; i++)
	{
		for (int j = 0; j < b.size; j++)
			c.V[i] += A.M[i][j] * b.V[j];
	}

	return c;
}


double y(double t)
{
	return t + yTemp;
}

double K(double t, double s, double x)
{
	return (t + x * x + 2) * (t + x * x + 2) / (t + s * s + 2);
}


double F(Vect t, Vect s, Vect x, Vect y, double h)
{
	double T = 0;

	double hHalf = h / 2;

	for (int i = 0; i < t.size; i++)
	{
		double temp = K(t.V[i], s.V[0], x.V[0]) * hHalf;
		for (int j = 1; j < t.size - 1; j++)
			temp += h * K(t.V[i], s.V[j], x.V[j]);
		temp += K(t.V[i], s.V[s.size - 1], x.V[x.size - 1]) * hHalf;

		double tmpForPow = temp - y.V[i];
		T += tmpForPow * tmpForPow;
	}
	return T;
}

double dK_dx(double t, double s, double x)
{
	return 4 * (t + x * x + 2) * x / (t + s * s + 2);
}

Vect grad(Vect t, Vect s, Vect x, Vect y, double h)
{
	//Vect T = InitVect(t.size);
	//T = EnterZero(T);

	double hHalf = h / 2;

	Vect T = InitZeroVect(t.size);

	for (int k = 0; k < T.size; k++)
	{
		double temp = 0;
		for (int i = 0; i < T.size; i++)
		{
			temp = K(t.V[i], s.V[0], x.V[0]) * hHalf;

			int tmpSizeIncr = T.size - 1;
			for (int j = 1; j < tmpSizeIncr; j++)
			{
				temp += h * K(t.V[i], s.V[j], x.V[j]);
			}
			temp += K(t.V[i], s.V[s.size - 1], x.V[x.size - 1]) * hHalf;
			temp = temp - y.V[i];
			
			temp = temp * dK_dx(t.V[i], s.V[k], x.V[k]) * h;

			if ((k != 0) && (k != tmpSizeIncr))
				temp = temp * 2;
		}
		T.V[k] = temp;
	}
	return T;
}


double mingrad(Vect t, Vect s, Vect x, Vect y, Vect g, double hh)
{
	double l = 0.0000001;
	double h = l;
	Vect x1 = InitVect(x.size);
	Vect x2 = InitVect(x.size);
	for (int i = 0; i < x.size; i++)
	{
		x1.V[i] = x.V[i];
		x2.V[i] = x.V[i] - h * g.V[i];
	}
	while (F(t, s, x1, y, hh) > F(t, s, x2, y, hh))
	{
		l += h;
		for (int i = 0; i < x.size; i++)
		{
			x1.V[i] = x2.V[i];
			x2.V[i] = x1.V[i] - h * g.V[i];
		}
	}
	return l;
}

int main()
{
	cout << "Hello world!" << endl;
	int N = 50;
	double h = ((double)(b - a)) / N;
	double lambda = 0.0001;
	N++;
	Vect t = InitVect(N);
	Vect s = InitVect(N);
	Vect f = InitVect(N);
	Vect x = InitVect(N);
	for (int i = 0; i < N; i++)
	{
		double ih = i*h;
		t.V[i] = ih;
		s.V[i] = ih;
		f.V[i] = y(ih);
		x.V[i] = y(ih);
	}
	int step = 0;
	Vect dx = InitVect(N);
	double nev = 0;
	//Print(t);
	do
	{
		//Print(x);
		//Print(f);
		lambda = mingrad(t, s, x, f, grad(t, s, x, f, h), h);
		cout << "Lambda ---- " << lambda << endl;
		dx = lambda * grad(t, s, x, f, h);
		x = x - dx;
		//Print(dx);
		//Print(x);
		nev = F(t, s, x, f, h);
		//nev = Norma(dx);
		cout << "!!!! Nevyazka - " << nev << "\n\n";
		step++;
		cout << "runtime = " << clock() / 1000.0 << endl;
		//_getch();
	} while (nev > eps);
	Print(x);
	cout << "runtime = " << clock() / 1000.0 << endl;
	cout << "Count of Iteration --- " << step;
	system("pause");
	return 0;
}
