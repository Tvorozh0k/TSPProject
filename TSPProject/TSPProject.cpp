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

	// MMAS

	map<string, any> mmas_params = {
		{"alpha", 1.0},
		{"beta", 2.5},
		{"rho", 0.02},
		{"n_ants", g.n()},
		{"tau0", 5.0 / alg1.len()},
		{"max_iter", 1500},
		{"a", 2.0 * (1 - p) / (g.n() - 2) / p},
		{"local_search_type", (string)"None"},
		{"n_jobs", 1}
	};

	double len_sum = 0, time_sum = 0;

	for (int i = 0; i < 10; ++i)
	{
		ACOSolver alg2("MMAS", mmas_params);

		system_clock::time_point start = system_clock::now();

		alg2.solve(g);

		system_clock::time_point end = system_clock::now();
		duration <double> delta = end - start;

		cout << "ACOSolver (затраченное время): " << delta.count() << "\n";

		alg2.print();

		len_sum += alg2.len();
		time_sum += delta.count();
	}

	cout << "MMAS AVG LEN: " << len_sum / 10 << "\n";
	cout << "MMAS AVG TIME: " << time_sum / 10 << " sec\n";

	// MMAS + local search

	map<string, any> params = {
		{"alpha", 1.0},
		{"beta", 2.0},
		{"rho", 0.2},
		{"n_ants", 20},
		{"tau0", 5.0 / alg1.len()},
		{"max_iter", 1500},
		{"a", 2.0 * (1 - p) / (g.n() - 2) / p},
		{"local_search_type", (string)"None"},
		{"local_search_tours", (string)"k-random"},
		{"k", 0},
		{"n_jobs", 1}
	};

	vector<string> opt = { "2-opt", "2.5-opt", "3-opt" };
	vector<int> k = { 20, 10, 10 };

	for (int i = 0; i < opt.size(); ++i)
	{
		params["local_search_type"] = opt[i];
		params["k"] = k[i];

		double len_sum = 0, time_sum = 0;

		for (int i = 0; i < 10; ++i)
		{
			ACOSolver alg2("MMAS", params);

			system_clock::time_point start = system_clock::now();

			alg2.solve(g);

			system_clock::time_point end = system_clock::now();
			duration <double> delta = end - start;

			cout << "ACOSolver (затраченное время): " << delta.count() << "\n";

			alg2.print();

			len_sum += alg2.len();
			time_sum += delta.count();
		}

		cout << "MMAS + ls AVG LEN: " << len_sum / 10 << "\n";
		cout << "MMAS + ls AVG TIME: " << time_sum / 10 << " sec\n";
	}

	return 0;
}