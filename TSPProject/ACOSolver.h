#pragma once
#include <string>
#include <map>
#include <any>

#include "KOptSolver.h"
#include "Ant.h"

/// <summary>
/// Виды муравьиных алгоритмов
/// </summary>
enum ACOType { AS, ElitistAS, RankBasedAS, MaxMinAS };

/// <summary>
/// Муравьиный алгоритм решения 
/// задачи коммивояжера
/// </summary>
class ACOSolver : public TSPSolver
{
public:

	/// <summary>
	/// Вид алгоритма
	/// </summary>
	ACOType _type;

	/// <summary>
	/// Гиперпараметры:
	/// 
	/// alpha - значимость феромонов,
	/// beta - значимость привлекательности ребер,
	/// rho - коэффициент испарения феромонов,
	/// tau0 - изначальное кол-во феромонов на ребрах,
	/// tau_min - нижняя граница возможного числа феромонов,
	/// tau_max - верхняя граница возможного числа феромонов,
	/// причем tau_min = a * tau_max 
	/// </summary>
	double _alpha, _beta, _rho, _tau0, _tau_min, _tau_max, _a;

	/// <summary>
	/// Гиперпараметры:
	/// 
	/// n_ants - число муравьев,
	/// max_iter - число поколений,
	/// e - коэффициент значимости лучшего тура
	/// </summary>
	int _n_ants, _max_iter, _e;

	/// <summary>
	/// Гиперпараметры:
	/// 
	/// bst_tour_type - тип лучшего тура,
	/// обновляющего феромоны на ребрах
	/// local_search_type - алгоритм локальной 
	/// оптимизации решения
	/// local_search_tours - какие туры оптимизируются
	/// </summary>
	string _bst_tour_type, _local_search_type, _local_search_tours;

	/// <summary>
	/// Конструктор
	/// </summary>
	/// <param name="g"> граф </param>
	/// <param name="alpha"> альфа </param>
	/// <param name="beta"> бета </param>
	/// <param name="rho"> ро </param>
	/// <param name="m"> число мурвьев </param>
	/// <param name="tau0"> начальное число феромонов </param>
	/// <param name="max_iter"> число поколений </param>
	ACOSolver(ACOType type, map<string, any> params) : _type(type)
	{
		// Повсеместные параметры
		_alpha = any_cast<double>(params["alpha"]);
		_beta = any_cast<double>(params["beta"]);
		_rho = any_cast<double>(params["rho"]);
		_tau0 = any_cast<double>(params["tau0"]);

		_n_ants = any_cast<int>(params["n_ants"]);
		_max_iter = any_cast<int>(params["max_iter"]);

		_local_search_type = any_cast<string>(params["local_search_type"]);

		// В случае наличия локального поиска оптимального решения
		if (_local_search_type != "None") _local_search_tours = any_cast<string>(params["local_search_tours"]);
		else  _local_search_tours = "None";

		// В случае элитной муравьиной системы
		if (type == ElitistAS)
		{
			_e = any_cast<int>(params["e"]);
			_bst_tour_type = "best-so-far";
		}

		// В случае ранговой муравьиной системы
		if (type == RankBasedAS)
		{
			_e = any_cast<int>(params["w"]);
			_bst_tour_type = "best-so-far";
		}

		// В случае максминной муравьиной системы
		if (type == MaxMinAS)
		{
			_e = 1;
			_a = any_cast<double>(params["a"]);
			_bst_tour_type = any_cast<string>(params["bst_tour_type"]);
		}
	}

	/// <summary>
	/// Решаем задачу
	/// </summary>
	/// <param name="g"> граф </param>
	void solve(Graph& g)
	{
		_len = INF;
		_n_cities = g.n();

		// В случае максминной муравьиной системы
		if (_type == MaxMinAS)
		{
			_tau_max = _tau0;
			_tau_min = _a * _tau_max;
		}
		else
		{
			_tau_min = 0;
			_tau_max = INF;
		}

		// матрица концентрации феромонов на ребрах графа
		vector<vector<double>> tau(_n_cities, vector<double>(_n_cities, _tau0));

		// матрица привлекательности ребер графа 
		// (все элементы возведены в степень -_beta)
		vector<vector<double>> eta_beta(_n_cities, vector<double>(_n_cities));

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			for (int j = 0; j < _n_cities; ++j)
				eta_beta[i][j] = pow(g[i][j], -_beta);

		// матрица весов (значимостей) каждого из ребер
		// (см. формула вероятности попадания муравья в тот или иной город)
		vector<vector<double>> weights(_n_cities, vector<double>(_n_cities));

		// все вершины графа
		vector<int> vertices(_n_cities);

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			vertices[i] = i;

		// генерируем муравьев 
		vector<Ant> ants(_n_ants);

		// visited[i] = 0 <=> вершина i НЕ посещена
		// visited[i] = 1 <=> вершина i посещена
		vector<vector<int>> visited(GridThreadsNum, vector<int>(_n_cities));

		// не посещенные вершины
		vector<vector<int>> choices(GridThreadsNum, vector<int>(_n_cities));

		for (int it = 0; it < _max_iter; ++it)
		{
			#pragma omp parallel for
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					weights[i][j] = pow(tau[i][j], _alpha) * eta_beta[i][j];

			// перемешиваем вершины графа
			shuffle(vertices.begin(), vertices.end(), gen);

			// лучшее за итерацию решение
			vector<vector<int>> iter_best_sol(GridThreadsNum);
			vector<int> iter_best_len(GridThreadsNum, INF);

			// муравей ищет решение	
			#pragma omp parallel for schedule(dynamic) 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				ants[i].solve(g, vertices[i], visited[thread_num], choices[thread_num], weights);

				for (int j = 0; j < _n_cities; ++j)
					visited[thread_num][j] = 0;

				// В случае наличия локального поиска оптимального решения
				if (_local_search_tours == "all")
				{
					KOptSolver ls_ant;
					ls_ant.solve(g, _local_search_type, ants[i].solution(), ants[i].len());

					if (ls_ant.len() < iter_best_len[thread_num])
					{
						iter_best_sol[thread_num] = ls_ant.solution();
						iter_best_len[thread_num] = ls_ant.len();
					}
				}
				else if (ants[i].len() < iter_best_len[thread_num])
				{
					iter_best_sol[thread_num] = ants[i].solution();
					iter_best_len[thread_num] = ants[i].len();
				}
			}

			// лучшее решение по потокам 
			int pos_min = min_element(iter_best_len.begin(), iter_best_len.end()) - iter_best_len.begin();

			// В случае наличия локального поиска оптимального решения
			if (_local_search_tours == "best_sol")
			{
				KOptSolver ls_best;
				ls_best.solve(g, _local_search_type, iter_best_sol[pos_min], iter_best_len[pos_min]);

				iter_best_sol[pos_min] = ls_best.solution();
				iter_best_len[pos_min] = ls_best.len();
			}

			// сохраняем решение
			if (iter_best_len[pos_min] < _len)
			{
				_solution = iter_best_sol[pos_min];
				_len = iter_best_len[pos_min];

				// В случае максминной муравьиной системы
				if (_type == MaxMinAS)
				{
					_tau_max = 1.0 / _rho / _len;
					_tau_min = _a * _tau_max;
				}
			}

			// феромоны испаряются 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					tau[i][j] = max(_tau_min, tau[i][j] * (1 - _rho));	

			// В случае ранговой муравьиной системы
			if (_type == RankBasedAS) sort(ants.begin(), ants.end());

			int iter;

			if (_type == RankBasedAS || _type == MaxMinAS) iter = _e - 1;
			else iter = _n_ants;

			// оставляем след феромонов на ребрах 
			for (int i = 0; i < iter; ++i)
			{
				vector<int> solution = ants[i].solution();
				double w = ((_type == RankBasedAS) ? (_e - i - 1.0) : 1.0) / ants[i].len();

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					tau[solution[j]][solution[j + 1]] += w;
			}

			// В случае НЕ классической муравьиной системы
			if (_type != AS)
			{
				vector<int> solution;
				double w;

				if (_bst_tour_type == "best-so-far")
				{
					solution = _solution;
					w = (double)_e / _len;
				}
				else 
				{
					solution = iter_best_sol[pos_min];
					w = (double)_e / iter_best_len[pos_min];
				}

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					tau[solution[j]][solution[j + 1]] = min(_tau_max, tau[solution[j]][solution[j + 1]] + w);
			}
		}
	}

	/// <summary>
	/// Выводим решение на экран
	/// </summary>
	void print()
	{
		cout << "ACOSolver (найденное решение): ";

		for (int i = 0; i < _n_cities; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", длина пути = " << _len << "\n";
	}
};