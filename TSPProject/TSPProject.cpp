#include "NNSolver.h"
#include "ACOSolver.h"

#include <chrono>
#include <iomanip>
using namespace chrono;

/// <summary>
/// Поиск оптимальных значений параметров 
/// муравьиного алгоритма
/// </summary>
/// <param name="g"> граф </param>
/// <param name="alg"> муравьиный алгоритм </param>
void search_opt_params(Graph& g, string alg)
{
	NNSolver alg1;
	alg1.solve(g);

	map<string, any> params = {
		{"alpha", 1.0},
		{"beta", 5.0},
		{"rho", 0.9},
		{"n_ants", 24},
		{"a", 0.5 / g.n()},
		{"tau0",  5.0 / alg1.len()},
		{"max_iter", 1500},
		{"local_search_type", (string)"2-opt"},
		{"local_search_tours", (string)"k-random"},
		{"k", 24},
		{"n_jobs", 2}
	};

	cout << "+------------------------------+\n";
	cout << "|         Find Params          |\n";
	cout << "+------------------------------+\n\n";

	double best_beta = 0, best_rho = 0, best_time = 0, best_avg_len = 0;
	int best_iter = 0, best_len = INF;

	for (double beta = 1.1; beta < 4.1; beta += 0.4)
		for (double rho = 0.3; rho < 1; rho += 0.2)
			for (int max_iter = 1500; max_iter <= 1500; max_iter += 500)
			{
				cout << "beta = " << beta << ", rho = " << rho << ", max_iter = " << max_iter << " : \n\n";

				params["beta"] = beta;
				params["rho"] = rho;
				params["max_iter"] = max_iter;

				cout << "+-----+------------+-----------+\n";
				cout << "| No. |  Решение   |   Время   |\n";
				cout << "+-----+------------+-----------+\n";

				double len_sum = 0, time_sum = 0;
				int iters = 10, cur_best_len = INF;

				for (int i = 0; i < iters; ++i)
				{
					ACOSolver alg2(alg, params);

					system_clock::time_point start = system_clock::now();

					alg2.solve(g);

					system_clock::time_point end = system_clock::now();
					duration <double> delta = end - start;

					//cout << left << "| " << setw(2) << i + 1 << "  |  " << setw(8) << alg2.len() << "  |   " << setw(6) << delta.count() << "  |\n";

					len_sum += alg2.len();
					time_sum += delta.count();
					cur_best_len = min(best_len, alg2.len());
				}

				//cout << "+-----+------------+-----------+\n";
				cout << left << "| BST |  " << setw(8) << cur_best_len << "  |   " << setw(6) << "----" << "  |\n";
				cout << left << "| AVG |  " << setw(8) << len_sum / iters << "  |   " << setw(6) << time_sum / iters << "  |\n";
				cout << "+-----+------------+-----------+\n\n";

				if (cur_best_len < best_len)
				{
					best_beta = beta;
					best_rho = rho;
					best_iter = max_iter;

					best_len = cur_best_len;
					best_avg_len = len_sum / iters;
					best_time = time_sum / iters;
				}
				else if (len_sum < best_avg_len * iters)
				{
					best_beta = beta;
					best_rho = rho;
					best_iter = max_iter;

					best_len = cur_best_len;
					best_avg_len = len_sum / iters;
					best_time = time_sum / iters;
				}
				else if (time_sum < best_time * iters)
				{
					best_beta = beta;
					best_rho = rho;
					best_iter = max_iter;

					best_len = cur_best_len;
					best_avg_len = len_sum / iters;
					best_time = time_sum / iters;
				}
			} 

	cout << "BEST PARAMS: beta = " << best_beta << ", rho = " << best_rho << ", max_iter = " << best_iter << ", len = " << best_len << "\n";
	cout << "RESULTS: best_len = " << best_len << ", avg_len = " << best_avg_len << ", time = " << best_time << "\n";
}

/// <summary>
/// Запускаем AS
/// </summary>
/// <param name="g"> граф </param>
void check_as(Graph &g)
{
	NNSolver alg1;
	alg1.solve(g);

	double p = pow(0.05, 1.0 / g.n());

	map<string, any> params = {
		{"alpha", 1.0},
		{"beta", 3.5},
		{"rho", 0.2},
		{"n_ants", 24},
		{"a", 0.5 / g.n()},
		{"tau0",  5.0 / alg1.len()},
		{"max_iter", 200},
		{"local_search_type", (string)"3-opt"},
		{"local_search_tours", (string)"k-random"},
		{"k", 12},
		{"n_jobs", 4}
	};

	cout << "+------------------------------+\n";
	cout << "|          Ant System          |\n";
	cout << "+------------------------------+\n\n";

	cout << "+-----+------------+-----------+\n";
	cout << "| No. |  Решение   |   Время   |\n";
	cout << "+-----+------------+-----------+\n";

	double len_sum = 0, time_sum = 0;
	int iters = 25, best_len = INF;

	for (int i = 0; i < iters; ++i)
	{
		ACOSolver alg2("MMAS", params);

		system_clock::time_point start = system_clock::now();

		alg2.solve(g);

		system_clock::time_point end = system_clock::now();
		duration <double> delta = end - start;

		cout << left << "| " << setw(3) << i + 1 << " |  " << setw(8) << alg2.len() << "  |   " << setw(6) << delta.count() << "  |\n";

		len_sum += alg2.len();
		time_sum += delta.count();
		best_len = min(best_len, alg2.len());
	}

	cout << "+-----+------------+-----------+\n";
	cout << left << "| BST |  " << setw(8) << best_len << "  |   " << setw(6) << "----" << "  |\n";
	cout << left << "| AVG |  " << setw(8) << len_sum / iters << "  |   " << setw(6) << time_sum / iters << "  |\n";
	cout << "+-----+------------+-----------+\n\n"; 
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << fixed << setprecision(2);

	Graph g("ftv170.atsp");
	
	//search_opt_params(g, "MMAS");
	check_as(g);

	/*

	NNSolver alg1;
	alg1.solve(g);

	map<string, any> params = {
		{"alpha", 1.0},
		{"beta", 5.0},
		{"rho", 0.2},
		{"n_ants", 24},
		{"a", 0.5 / g.n()},
		{"tau0",  5.0 / alg1.len()},
		{"max_iter", 2000},
		{"local_search_type", (string)"2-opt"},
		{"local_search_tours", (string)"k-random"},
		{"k", 24},
		{"n_jobs", 1}
	};

	ACOSolver alg2("MMAS", params);
	alg2.solve(g);

	vector<int> sol = alg2.solution();

	int sum = 0;

	for (int i = 0; i < sol.size() - 1; ++i)
		sum += g[sol[i]][sol[i + 1]];

	cout << "sum = " << sum << ", ans = " << alg2.len() << "\n"; */


	/*
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
		{"n_jobs", 8}
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
	cout << "MMAS AVG TIME: " << time_sum / 10 << " sec\n"; */

	/*

	// MMAS + local search

	map<string, any> params = {
		{"alpha", 1.0},
		{"beta", 2.0},
		{"rho", 0.2},
		{"n_ants", 24},
		{"tau0", 5.0 / alg1.len()},
		{"max_iter", 1500},
		{"a", 2.0 * (1 - p) / (g.n() - 2) / p},
		{"local_search_type", (string)"None"},
		{"local_search_tours", (string)"k-random"},
		{"k", 0},
		{"n_jobs", 8}
	};

	vector<string> opt = { "2-opt", "2.5-opt", "3-opt" };
	vector<int> k = { 24, 16, 16 };

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
	} */

	return 0;
}
