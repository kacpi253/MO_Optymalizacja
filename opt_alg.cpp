#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
		for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
		Xopt.fit_fun(ff, ud1, ud2);							// obliczamy wartoœæ funkcji celu
		Xopt.ud = (Xopt.y);									// wartoœæ funckji celu zapamiêtujemy w macierzy ud
		if (Xopt.y < epsilon)								// sprawdzamy 1. kryterium stopu
		{
			Xopt.flag = 1;
			return Xopt;
		}
		solution X;
		while (true)
		{
			X = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
				X.x(i) = (ub(i) - lb(i)) * X.x(i) + lb(i);
			X.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (X.y < Xopt.y)
			{
				Xopt = X;
				if (Xopt.y < epsilon)							// sprawdzmy 1. kryterium stopu
				{
					Xopt.flag = 1;								// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
					break;
				}
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalnej liczby wywo³añ funkcji celu
				break;
			}
			Xopt.ud.add_row(Xopt.y);							// wartoœæ funckji dodajemy jako kolejny wiersz macierzy ud
		}
		Xopt.ud.add_row(Xopt.y);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)			//Ekspansja lab 1
{
	try
	{
		double* p = new double[2] { 0, 0 };
												// Owijamy liczby w obiekty klasy solution. 
												// Pozwala to na uniwersalne wywo³ywanie f. wielowymiarowych i liczenie f_calls
		solution X0(x0);
		solution X1(x0 + d);

		X0.fit_fun(ff, ud1, ud2);				//obliczenie wartoœci funkcji celu
		X1.fit_fun(ff, ud1, ud2);

		
		if (X0.y == X1.y)						// Krok 3-5 wg pseudokodu
		{
			p[0] = m2d(X0.x);					//matrix to double, dolna granica przedzialu
			p[1] = m2d(X1.x);					//górna granica przedzia³u
			return p;
		}

												// Krok 6-12 wg pseudokodu
		if (X1.y > X0.y)						//matrix jest przeciazona wiec mo¿enyt ak zrobiæ
		{
			d = -d;
			X1.x = x0 + d;
			X1.fit_fun(ff, ud1, ud2);

			if (X1.y >= X0.y)
			{
				p[0] = m2d(X1.x);
				p[1] = x0 - d;
				return p;
			}
		}

												// Zmienne do pêtli repeat-until
		solution X_prev = X0;
		solution X_curr = X1;
		solution X_next;

		int i = 0;

												// Pêtla wg kroku 13-19 z pseudokodu
		while (true)
		{
			if (solution::f_calls > Nmax)
				throw string("Przekroczono maksymalna liczbe wywolan funkcji celu w metodzie ekspansji.");

			i++;
			X_next.x = x0 + pow(alpha, i) * d;
			X_next.fit_fun(ff, ud1, ud2);

			if (X_curr.y <= X_next.y)
				break;

			X_prev = X_curr;
			X_curr = X_next;
		}

												// Krok 20-23 wg pseudokodu
		if (d > 0)
		{
			p[0] = m2d(X_prev.x);
			p[1] = m2d(X_next.x);
		}
		else
		{
			p[0] = m2d(X_next.x);
			p[1] = m2d(X_prev.x);
		}
		
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)		//Algorytm fibbonaciego
{
	try
	{
		solution Xopt;
		
		//inicjalizacja ciagu fibbonaciego, tablica 100 elemntow
		double fi[100];
		fi[0] = 1;
		fi[1] = 1;
		int k = 1;

		// Krok 1: znajdŸ najmniejsz¹ liczbê k spe³niaj¹c¹ nierównoœæ fi_k > (b-a)/epsilon
		while (fi[k] <= (b - a) / epsilon)
		{
			k++;
			fi[k] = fi[k - 1] + fi[k - 2];
			if (k >= 99) break; // Zabezpieczenie przed wyjœciem poza tablicê
		}

		// Krok 2: a_0 = a, b_0 = b
		double a_val = a;
		double b_val = b;

		// Krok 3, 4: oblicz c_0, d_0
		double c_val = b_val - (fi[k - 1] / fi[k]) * (b_val - a_val);
		double d_val = a_val + b_val - c_val;

		// Owijamy w obiekty solution i obliczamy po raz pierwszy funkcjê celu
		solution C(c_val), D(d_val);
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);

		Xopt.ud = matrix(b_val - a_val); // Inicjalizacja ud d³ugoœci¹ pocz¹tkow¹  // // //

		// Krok 5: pêtla
		for (int i = 0; i <= k - 3; ++i)
		{
			if (C.y < D.y) // Krok 6
			{
				// Krok 7, 8: a zostaje bez zmian, b = d
				b_val = d_val;

				// OPtymalizacja - stary punkt C staje sê=iê nowym punkntem D - nie musiy luczyæ funkcji 2 raz, przypisujemy obiekt
				d_val = c_val;
				D = C;

				// Krok 13: Liczymy tylko nowe C
				c_val = b_val - (fi[k - i - 2] / fi[k - i - 1]) * (b_val - a_val);
				C.x = c_val;
				C.fit_fun(ff, ud1, ud2);
			}
			else
			{
				// Krok 10, 11: b zostaje bez zmian, a = c
				a_val = c_val;

				// OPtymalizacja: Stary punkt D staje siê nowym C
				c_val = d_val;
				C = D;

				// Krok 14: Liczymy tylko nowe D
				d_val = a_val + b_val - c_val;
				D.x = d_val;
				D.fit_fun(ff, ud1, ud2);
			}

			Xopt.ud.add_row(b_val - a_val);  // // //
		}

		// Krok 16: Zwróæ c^(i+1)
		Xopt = C;
		
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)		//Algorytm Lagrange
{
	try
	{
		solution Xopt;
		// Nie ma c w nag³ówku wiêc przyjmujemy œrodek przedzia³u 
		double c = (a + b) / 2.0;

		// Owijamy w obiekty solution, ¿eby nie liczyæ wielokrotnie wartoœci funkcji
		solution A(a), B(b), C(c), D;
		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);

		double d_prev = c; // Do warunku stopu w pierwszej iteracji
		double d_val;


		Xopt.ud = matrix(b - a); // Inicjalizacja ud d³ugoœci¹ pocz¹tkow¹ // // //

		while (true)
		{
			// Wyci¹gamy wartoœci wspó³rzêdnych i funkcji z obiektów
			double a_val = m2d(A.x);
			double b_val = m2d(B.x);
			double c_val = m2d(C.x);

			double fa = m2d(A.y);
			double fb = m2d(B.y);
			double fc = m2d(C.y);

			// Krok 4 i 5 z pseudokodu: liczymy "l" oraz "m"
			// U¿ywamy mno¿enia (a_val * a_val) zamiast wbudowanej funkcji pow(), 
			// aby C++ nie mia³ w¹tpliwoœci jakich typów danych u¿ywamy.
			double l = fa * (b_val * b_val - c_val * c_val) +
				fb * (c_val * c_val - a_val * a_val) +
				fc * (a_val * a_val - b_val * b_val);

			double m = fa * (b_val - c_val) +
				fb * (c_val - a_val) +
				fc * (a_val - b_val);

			// Krok 6-8: zabezpieczenie przed m <= 0
			if (m <= 0)
				throw string("Mianownik m <= 0 w metodzie Lagrange'a");

			// Krok 9: obliczamy d
			d_val = 0.5 * l / m;

			D.x = d_val;
			D.fit_fun(ff, ud1, ud2); // Liczymy funkcjê tylko dla nowego punktu D
			double fd = m2d(D.y);

			// Krok 10-19: d znajduje siê w lewej po³ówce przedzia³u (miêdzy a i c)
			if (a_val < d_val && d_val < c_val)
			{
				if (fd < fc)
				{
					B = C; // Stare C staje siê nowym B (b = c)
					C = D; // Nowe D staje siê nowym C (c = d)
				}
				else
				{
					A = D; // Nowe D staje siê nowym A (a = d)
				}
			}
			// Krok 20-30: d znajduje siê w prawej po³ówce (miêdzy c i b)
			else if (c_val < d_val && d_val < b_val)
			{
				if (fd < fc)
				{
					A = C; // Stare C staje siê nowym A (a = c)
					C = D; // Nowe D staje siê nowym C (c = d)
				}
				else
				{
					B = D; // Nowe D staje siê nowym B (b = d)
				}
			}
			else
			{
				// Krok 31-33: D poza granicami
				throw string("Punkt d wyszedl poza przedzial poszukiwan");
			}

			// Krok 36-38: Sprawdzenie licznika wywo³añ
			if (solution::f_calls > Nmax)
				throw string("Przekroczono maksymalna liczbe wywolan funkcji celu (Nmax)");

			// Krok 39: Kryterium stopu. Pêtla koñczy siê, gdy przedzia³ [a, b] 
			// jest mniejszy od epsilon ALBO kolejne przybli¿enie 'd' przesta³o siê zmieniaæ.
			// Fabs to wyci¹gniêcie wartoœci bezwzglêdnej (modu³u) dla typu double.
			if (m2d(B.x) - m2d(A.x) < epsilon || fabs(d_val - d_prev) < gamma)
			{
				Xopt = D;

				return Xopt;

			}

			Xopt.ud.add_row(m2d(B.x) - m2d(A.x)); // // //
			d_prev = d_val;
		}
	}

	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
