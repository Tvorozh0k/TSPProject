#include "NNSolver.h"
#include "ACOSolver.h"

#include <chrono>
using namespace chrono;

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << fixed << setprecision(6);

	Graph g("ft53.txt");

	NNSolver alg1;
	alg1.solve(g);

	double p = pow(0.05, 1.0 / g.n());

	map<string, any> params = {
		{"alpha", 1.0},
		{"beta", 2.0},
		{"rho", 0.2},
		{"n_ants", 25},
		{"tau0", 5.0 / alg1.len()},
		{"max_iter", 1500},
		{"a", (1 - p) / p / (g.n() / 2.0 - 1)},
		{"bst_tour_type", (string)"best-so-far"},
		{"local_search_type", (string)"2-opt"},
		{"local_search_tours", (string)"all"}
	};

	double len_sum = 0, time_sum = 0;

	for (int i = 0; i < 100; ++i)
	{
		ACOSolver alg2(MaxMinAS, params);

		system_clock::time_point start = system_clock::now();

		alg2.solve(g);

		system_clock::time_point end = system_clock::now();
		duration <double> delta = end - start;

		cout << "ACOSolver (затраченное время): " << delta.count() << "\n";

		alg2.print();

		len_sum += alg2.len();
		time_sum += delta.count();
	}

	cout << "AVG LEN: " << len_sum / 100 << "\n";
	cout << "AVG TIME: " << time_sum / 100 << " sec\n"; 

	return 0;
}