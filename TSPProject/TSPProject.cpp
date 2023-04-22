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
	alg1.print();

	map<string, double> params = {
		{"alpha", 1},
		{"beta", 3},
		{"rho", 0.9},
		{"m", g.n()},
		{"tau0", 15 / 0.9 / alg1.len()},
		{"max_iter", 1000},
		{"w", 6}
	};

	double len_sum = 0, time_sum = 0;

	for (int i = 0; i < 10; ++i)
	{
		ACOSolver alg2(RankBasedAS, params);

		system_clock::time_point start = system_clock::now();

		alg2.solve(g);

		system_clock::time_point end = system_clock::now();
		duration <double> delta = end - start;

		cout << "ACOSolver (затраченное время): " << delta.count() << "\n";

		alg2.print();

		len_sum += alg2.len();
		time_sum += delta.count();
	}

	cout << "AVG LEN: " << len_sum / 10 << "\n";
	cout << "AVG TIME: " << time_sum / 10 << " sec\n";

	return 0;
}