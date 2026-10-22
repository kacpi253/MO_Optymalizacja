/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 15.09.2026
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	double epsilon;												// dok³adnoœæ
	int Nmax;													// maksymalna liczba wywo³añ funkcji celu
	matrix lb, ub;												// dolne oraz górne ograniczenie
	solution opt;												// rozwi¹zanie optymalne znalezione przez algorytm
	if (1) // Funkcja testowa
	{
		epsilon = 1e-2;
		Nmax = 10000;
		lb = matrix(2, 1, -5);
		ub = matrix(2, 1, 5);
		matrix a(2, 1);											// dok³adne rozwi¹zanie optymalne
		a(0) = -1;
		a(1) = 2;
		opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie liczników
		int* n = get_size(opt.ud);
		cout << n[0] << '\t' << n[1] << endl << endl;
	}
	if (1) // Wahad³o
	{
		Nmax = 1000;
		epsilon = 1e-2;
		lb = 0, ub = 5;
		double teta_opt = 1;									// maksymalne wychylenie wahad³a
		opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie liczników

		// Zapis zmian wartoœci funckji celu do pliku
		ofstream Sout("optymalicja_lab0.csv");					// definiujemy strumieñ do pliku .csv
		Sout << opt.ud;											// zapisujemy wyniki w pliku
		Sout.close();											// zamykamy strumieñ

		// Zapis symulacji do pliku csv
		matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
			MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
		matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
		Sout.open("symulacja_lab0.csv");						// otwieramy plik .csv
		Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
		Sout.close();											// zamykamy strumieñ
		Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
		Y[1].~matrix();
	}
}

void lab1()
{
	cout << "================== LAB 1: START ==================" << endl;

	// ---------------------------------------------------------------
	// CZÊŒÆ 1: FUNKCJA TESTOWA - 100 POWTÓRZEÑ DLA 3 WARTOŒCI ALFA
	// ---------------------------------------------------------------
	double alphas[3] = { 1.5, 2.0, 2.5 };
	double epsilon = 1e-5;
	double gamma = 1e-5;
	int Nmax = 1000;
	double d = 1.0;

	ofstream out_tab1("tabela1.csv");
	out_tab1 << "Alfa;Lp;x0;a;b;f_calls_exp;x_fib;y_fib;f_calls_fib;typ_fib;x_lag;y_lag;f_calls_lag;typ_lag\n";

	// Generowanie 100 powtarzalnych losowych punktów startowych x0 z przedzia³u [-100, 100]
	default_random_engine gen(12345);
	uniform_real_distribution<double> dist(-100.0, 100.0);

	double x0_arr[100];
	for (int i = 0; i < 100; ++i)
		x0_arr[i] = dist(gen);

	for (int a_idx = 0; a_idx < 3; ++a_idx)
	{
		double alpha = alphas[a_idx];
		cout << "\nLiczenie testowe dla alfa = " << alpha << "..." << endl;

		for (int i = 0; i < 100; ++i)
		{
			double x0 = x0_arr[i];

			// 1. Ekspansja
			solution::clear_calls();
			double* p = expansion(ff1T, x0, d, alpha, Nmax);
			int calls_exp = solution::f_calls;

			// 2. Metoda Fibonacciego w znalezionym przedziale
			solution::clear_calls();
			solution opt_fib = fib(ff1T, p[0], p[1], epsilon);
			int calls_fib = solution::f_calls;
			string typ_fib = (m2d(opt_fib.y) < -0.9) ? "Globalne" : "Lokalne";

			// 3. Metoda Lagrange'a w znalezionym przedziale
			solution::clear_calls();
			solution opt_lag;
			int calls_lag = 0;
			string typ_lag = "";
			try
			{
				opt_lag = lag(ff1T, p[0], p[1], epsilon, gamma, Nmax);
				calls_lag = solution::f_calls;
				typ_lag = (m2d(opt_lag.y) < -0.9) ? "Globalne" : "Lokalne";
			}
			catch (...)
			{
				typ_lag = "Brak zbieznosci";
			}

			out_tab1 << alpha << ";" << (i + 1) << ";" << x0 << ";"
				<< p[0] << ";" << p[1] << ";" << calls_exp << ";"
				<< m2d(opt_fib.x) << ";" << m2d(opt_fib.y) << ";" << calls_fib << ";" << typ_fib << ";"
				<< (typ_lag == "Brak zbieznosci" ? 0.0 : m2d(opt_lag.x)) << ";"
				<< (typ_lag == "Brak zbieznosci" ? 0.0 : m2d(opt_lag.y)) << ";"
				<< calls_lag << ";" << typ_lag << "\n";

			delete[] p;
		}
	}
	out_tab1.close();
	cout << "Wyniki zapisane do pliku 'tabela1.csv'." << endl;

	// ---------------------------------------------------------------
	// CZÊŒÆ 2: WYKRES D£UGOŒCI PRZEDZIA£U [a, b] BEZ EKSPANSJI
	// ---------------------------------------------------------------
	cout << "\nGenerowanie danych do wykresu (przedzial [-100, 100])..." << endl;
	solution::clear_calls();
	solution f_plot = fib(ff1T, -100.0, 100.0, epsilon);

	solution::clear_calls();
	solution l_plot = lag(ff1T, -100.0, 100.0, epsilon, gamma, Nmax);

	ofstream out_wykres("wykres.csv");
	out_wykres << "Iteracja;Fibonacci_ba;Lagrange_ba\n";
	int n_f = get_len(f_plot.ud);
	int n_l = get_len(l_plot.ud);
	int max_iter = max(n_f, n_l);

	for (int i = 0; i < max_iter; ++i)
	{
		out_wykres << (i + 1) << ";";
		if (i < n_f) out_wykres << f_plot.ud(i);
		out_wykres << ";";
		if (i < n_l) out_wykres << l_plot.ud(i);
		out_wykres << "\n";
	}
	out_wykres.close();
	cout << "Dane do wykresu zapisane w 'wykres.csv'." << endl;

	// ---------------------------------------------------------------
	// CZÊŒÆ 3: PROBLEM RZECZYWISTY (ZBIORNIKI)
	// ---------------------------------------------------------------
	cout << "\nOptymalizacja przekroju otworu zbiornika DA..." << endl;
	double DA_a = 1.0;
	double DA_b = 100.0;
	double eps_real = 1e-3;
	double gamma_real = 1e-3;

	// Optymalizacja Fibonaccim
	solution::clear_calls();
	solution opt_real_fib = fib(ff1R, DA_a, DA_b, eps_real);
	int calls_real_fib = solution::f_calls;

	// Optymalizacja Lagrange'em
	solution::clear_calls();
	solution opt_real_lag = lag(ff1R, DA_a, DA_b, eps_real, gamma_real, 200);
	int calls_real_lag = solution::f_calls;

	cout << "Fibonacci: DA* = " << m2d(opt_real_fib.x) << " cm^2, y* = " << m2d(opt_real_fib.y) << ", f_calls = " << calls_real_fib << endl;
	cout << "Lagrange:  DA* = " << m2d(opt_real_lag.x) << " cm^2, y* = " << m2d(opt_real_lag.y) << ", f_calls = " << calls_real_lag << endl;

	// Symulacja dla znalezionych parametrów DA* (do arkusza Symulacja)
	matrix Y0(3, 1);
	Y0(0) = 5.0; Y0(1) = 1.0; Y0(2) = 20.0;

	double DA_fib_m2 = m2d(opt_real_fib.x) * 1e-4;
	double DA_lag_m2 = m2d(opt_real_lag.x) * 1e-4;

	matrix* Y_fib = solve_ode(df1, 0, 1.0, 2000, Y0, NAN, DA_fib_m2);
	matrix* Y_lag = solve_ode(df1, 0, 1.0, 2000, Y0, NAN, DA_lag_m2);

	ofstream out_sym("symulacja.csv");
	out_sym << "t;VA_fib;VA_lag;VB_fib;VB_lag;TB_fib;TB_lag\n";
	int steps = get_len(Y_fib[0]);
	for (int i = 0; i < steps; ++i)
	{
		out_sym << Y_fib[0](i) << ";"
			<< Y_fib[1](i, 0) << ";" << Y_lag[1](i, 0) << ";"
			<< Y_fib[1](i, 1) << ";" << Y_lag[1](i, 1) << ";"
			<< Y_fib[1](i, 2) << ";" << Y_lag[1](i, 2) << "\n";
	}
	out_sym.close();
	cout << "Przebieg symulacji zapisany do pliku 'symulacja.csv'." << endl;

	Y_fib[0].~matrix(); Y_fib[1].~matrix();
	Y_lag[0].~matrix(); Y_lag[1].~matrix();

	cout << "================== LAB 1: KONIEC ==================" << endl;
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
