#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wspó³rzêdne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	int n = get_len(Y[0]);									// d³ugoœæ rozwi¹zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad³a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartoœæ funkcji celu (ud1 to za³o¿one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po³o¿enia to prêdkoœæ
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prêdkoœci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {			//Funkcja celu lab1 Testowa
	
	matrix y;											//bêdzie zawiera³a funkcjê celu

	double x_val = x(0);								// wartoœæ zmiennej x z macierzy 1x1

	//wzór funkcji zapisana matematycznym równaniem , M_PI - sta³a pi
	y = -cos(0.1 * x_val) * exp(-pow(0.1 * x_val - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x_val, 2);

	return y;
}


// Uk³ad równañ ró¿niczkowych opisuj¹cy bilans cieczy i temperatury w zbiornikach A i B
matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(3, 1); // dY(0) = dVA/dt, dY(1) = dVB/dt, dY(2) = dTB/dt

	// Parametry sta³e modelu z konspektu
	double PA = 1.5;             // Pole podstawy zbiornika A [m^2]
	double PB = 1.0;             // Pole podstawy zbiornika B [m^2]
	double DB = 36.5665e-4;      // Przekrój otworu odp³ywowego zbiornika B [m^2]
	double a = 0.98;             // Wspó³czynnik lepkoœci cieczy
	double b = 0.63;             // Wspó³czynnik zwê¿enia strumienia
	double g = 9.81;             // Przyspieszenie ziemskie [m/s^2]
	double F_in_B = 10.0 / 1000.0; // Dop³yw zewnêtrzny do B: 10 l/s -> 0.01 m^3/s
	double T_in_B = 20.0;        // Temperatura dop³ywu zewnêtrznego do B [°C]
	double TA_const = 90.0;      // Temperatura wody w zbiorniku A [°C] (sta³a, brak dop³ywu)

	double DA = m2d(ud2);        // Przekrój DA przekazany w m^2 przez ud2

	double VA = Y(0);
	double VB = Y(1);
	double TB = Y(2);

	// 1. Wyp³yw ze zbiornika A do B
	double FA_out = 0.0;
	if (VA > 0.0)
		FA_out = a * b * DA * sqrt(2.0 * g * (VA / PA));

	dY(0) = -FA_out;

	// 2. Wyp³yw ze zbiornika B
	double FB_out = 0.0;
	if (VB > 0.0)
		FB_out = a * b * DB * sqrt(2.0 * g * (VB / PB));

	dY(1) = FA_out + F_in_B - FB_out;

	// 3. Bilans energii i zmiana temperatury w B
	// dT/dt = (V_in / V) * (T_in - T) dla ka¿dego ze strumieni wp³ywu
	if (VB > 1e-6)
		dY(2) = (FA_out / VB) * (TA_const - TB) + (F_in_B / VB) * (T_in_B - TB);
	else
		dY(2) = 0.0;

	return dY;
}

// Funkcja celu dla problemu rzeczywistego (zbiorniki)
matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	double DA_cm2 = m2d(x);      // Zmienna decyzyjna w cm^2
	double DA_m2 = DA_cm2 * 1e-4; // Przeliczenie na m^2 do równañ fizycznych

	// Warunki pocz¹tkowe: VA(0) = 5 m^3, VB(0) = 1 m^3, TB(0) = 20 °C
	matrix Y0(3, 1);
	Y0(0) = 5.0;
	Y0(1) = 1.0;
	Y0(2) = 20.0;

	// Symulacja: czas od 0 do 2000 s z krokiem dt = 1 s
	matrix* Y = solve_ode(df1, 0, 1.0, 2000, Y0, NAN, DA_m2);

	int n = get_len(Y[0]);
	double TB_max = Y[1](0, 2); // Kolumna 2 w macierzy Y[1] to TB

	for (int i = 1; i < n; ++i)
	{
		if (Y[1](i, 2) > TB_max)
			TB_max = Y[1](i, 2);
	}

	// Cel: maksymalna temperatura TB ma wynosiæ 50 °C
	y = fabs(TB_max - 50.0);

	Y[0].~matrix();
	Y[1].~matrix();

	return y;
}