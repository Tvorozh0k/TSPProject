#pragma once
#include <string>
#include <map>
#include <any>

#include <chrono>
using namespace chrono;

#include "KOptSolver.h"
#include "Ant.h"

ofstream fout("iter_time.txt");

/// <summary>
/// Муравьиный алгоритм решения 
/// задачи коммивояжера
/// </summary>
class ACOSolver : public TSPSolver
{
	/// <summary>
	/// Вид алгоритма:
	/// <para> - AS (Ant System) </para> 
	/// <para> - EAS (Elitist Ant System) </para>
	/// <para> - ASRank (Rank-based Ant System) </para>
	/// <para> - MMAS (Max-Min Ant System) </para>
	/// </summary>
	string _type;

	/// <summary>
	/// Повсеместные гиперпараметры:
	/// <para> - alpha (значимость феромонов) </para>
	/// <para> - beta (значимость привлекательности ребер) </para>
	/// <para> - rho (коэффициент испарения феромонов) </para>
	/// <para> - tau0 (изначальное кол-во феромонов на ребрах) </para>
	/// </summary>
	double _alpha, _beta, _rho, _tau0;

	/// <summary>
	/// Повсеместные гиперпараметры:
	/// <para> - n_ants (число муравьев) </para>
	/// <para> - max_iter (число поколений) </para>
	/// <para> - n_jobs (число используемых потоков) </para>
	/// </summary>
	int _n_ants, _max_iter, _n_jobs;

	/// <summary>
	/// Гиперпараметр w:
	/// <para> - EAS (значимость лучшего тура) </para>
	/// <para> - ASRank (число лучших муравьев, учитываемых
	/// при перераспределении феромонов на ребрах) </para>
	/// </summary>
	int _w;

	/// <summary>
	/// Максминная муравьиная система (MMAS):
	/// <para> - tau_min (нижняя граница возможного числа феромонов) </para>
	/// <para> - tau_max (верхняя граница возможного числа феромонов) </para>
	/// <para> - a (отношение tau_max / tau_min) </para>
	/// </summary>
	double _tau_min, _tau_max, _a;

	/// <summary>
	/// Локальная оптимизация решений, найденных муравьями:
	/// <para> - None (без локальной оптимизации) </para>
	/// <para> - 2-opt </para>
	/// <para> - 2.5-opt </para>
	/// <para> - 3-opt </para>
	/// </summary>
	string _local_search_type;

	/// <summary>
	/// Локальная оптимизация решений, найденных муравьями:
	/// <para> - k-best (k лучших муравьев) </para>
	/// <para> - k-random (k случайных муравьев) </para>
	/// </summary>
	string _local_search_tours;

	/// <summary>
	/// Гиперпараметр k:
	/// Число муравьев, решение которых будет локально
	/// оптимизировано
	/// </summary>
	int _k;

	/// <summary>
	/// Вспомогательные матрицы для работы алгоритма:
	/// <para> - tau (матрица концентрации феромонов на ребрах) </para>
	/// <para> - eta_beta (матрица привлекательности ребер,
	/// возведенных в степень -beta) </para>
	/// <para> - weights (матрица значимостей ребер, феромоны + 
	/// привлекательность) </para>
	/// </summary>
	vector<vector<double>> _tau, _eta_beta, _weights;

	/// <summary>
	/// Вспомогательный массив для работы алгоритма:
	/// <para> - vertices (массив вершин графа) </para> 
	/// </summary>
	vector<int> _vertices;

	/// <summary>
	/// Вспомогательный массив для работы алгоритма:
	/// <para> - ants (массив муравьев) </para>
	/// </summary>
	vector<Ant> _ants;
	
	/// <summary>
	/// Вспомогательные массивы для работы алгоритма:
	/// <para> - choices (массив непосещенных муравьем городов) </para>
	/// <para> - visited (массив с флагами посещен / не посещен ли муравьем город) </para>
	/// </summary>
	vector<vector<int>> _choices, _visited;

	/// <summary>
	/// (Elitist) Ant System
	/// </summary>
	/// <param name="g"> граф </param>
	void elitist_ant_system(Graph& g)
	{
		for (int it = 0; it < _max_iter; ++it)
		{
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_weights[i][j] = pow(_tau[i][j], _alpha) * _eta_beta[i][j];

			// перемешиваем вершины графа
			shuffle(_vertices.begin(), _vertices.end(), gen);

			// лучшее за итерацию решение
			vector<int> pos(_n_jobs);

			// муравей ищет решение	
			#pragma omp parallel for 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				_ants[i].solve(g, _vertices[i], _visited[thread_num], _choices[thread_num], _weights);

				for (int j = 0; j < _n_cities; ++j)
					_visited[thread_num][j] = 0;

				if (_ants[i] < _ants[pos[thread_num]]) pos[thread_num] = i;
			}

			for (int i = 0; i < _n_jobs; ++i)
				if (_ants[pos[i]].len() < _len)
				{
					_solution = _ants[pos[i]].solution();
					_len = _ants[pos[i]].len();
				}

			// феромоны испаряются 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_tau[i][j] *= (1 - _rho); 

			// оставляем след феромонов на ребрах 
			for (int i = 0; i < _n_ants; ++i)
			{
				vector<int> solution = _ants[i].solution();
				double w = 1.0 / _ants[i].len();

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					_tau[solution[j]][solution[j + 1]] += w;
			}

			// В случае элитной муравьиной системы
			if (_type == "EAS")
			{
				double w = (double)_w / _len;

				#pragma omp parallel for
				for (int j = 0; j < _n_cities; ++j)
					_tau[_solution[j]][_solution[j + 1]] += w;
			}
		}
	}

	/// <summary>
	/// Rank-Based Ant System
	/// </summary>
	/// <param name="g"> граф </param>
	void rank_based_ant_system(Graph& g)
	{
		for (int it = 0; it < _max_iter; ++it)
		{
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_weights[i][j] = pow(_tau[i][j], _alpha) * _eta_beta[i][j];

			// перемешиваем вершины графа
			shuffle(_vertices.begin(), _vertices.end(), gen);

			// лучшее за итерацию решение
			vector<int> pos(_n_jobs);

			// муравей ищет решение	
			#pragma omp parallel for 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				_ants[i].solve(g, _vertices[i], _visited[thread_num], _choices[thread_num], _weights);

				for (int j = 0; j < _n_cities; ++j)
					_visited[thread_num][j] = 0;

				if (_ants[i] < _ants[pos[thread_num]]) pos[thread_num] = i;
			}

			for (int i = 0; i < _n_jobs; ++i)
				if (_ants[pos[i]].len() < _len)
				{
					_solution = _ants[pos[i]].solution();
					_len = _ants[pos[i]].len();
				}

			// феромоны испаряются 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_tau[i][j] *= (1 - _rho);

			// сортируем муравьев по длине решения
			sort(_ants.begin(), _ants.end());

			// оставляем след феромонов на ребрах 
			for (int i = 0; i < _w - 1; ++i)
			{
				vector<int> solution = _ants[i].solution();
				double w = (_w - i - 1.0) / _ants[i].len();

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					_tau[solution[j]][solution[j + 1]] += w;
			}

			// Учитываем лучшее решение
			double w = (double)_w / _len;

			#pragma omp parallel for
			for (int j = 0; j < _n_cities; ++j)
				_tau[_solution[j]][_solution[j + 1]] += w;
		}
	}

	/// <summary>
	/// Max-Min Ant System
	/// </summary>
	/// <param name="g"> граф </param>
	void max_min_ant_system(Graph& g)
	{
		_tau_max = _tau0;
		_tau_min = _a * _tau_max;

		for (int it = 0; it < _max_iter; ++it)
		{
			system_clock::time_point start = system_clock::now();

			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_weights[i][j] = pow(_tau[i][j], _alpha) * _eta_beta[i][j];

			// перемешиваем вершины графа
			shuffle(_vertices.begin(), _vertices.end(), gen);

			// локальные оптимизации решений
			vector<KOptSolver> ls_ants(_k);

			// лучшее за итерацию решение
			vector<int> pos(_n_jobs);

			// муравей ищет решение	
			#pragma omp parallel for 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				_ants[i].solve(g, _vertices[i], _visited[thread_num], _choices[thread_num], _weights);

				for (int j = 0; j < _n_cities; ++j)
					_visited[thread_num][j] = 0;

				if (_local_search_type == "None" && _ants[i] < _ants[pos[thread_num]]) 
					pos[thread_num] = i;
			}

			// Если есть локальная оптимизация решений
			if (_local_search_type != "None")
			{
				if (_local_search_tours == "k-best") sort(_ants.begin(), _ants.end());
				else shuffle(_ants.begin(), _ants.end(), gen);

				#pragma omp parallel for 
				for (int i = 0; i < _k; ++i)
				{
					int thread_num = omp_get_thread_num();

					ls_ants[i].solve(g, _local_search_type, _ants[i].solution(), _ants[i].len());

					if (ls_ants[i] < ls_ants[pos[thread_num]]) pos[thread_num] = i;
				}

				// Обновляем лучшее решение
				for (int i = 0; i < _n_jobs; ++i)
					if (ls_ants[pos[i]].len() < _len)
					{
						_solution = ls_ants[pos[i]].solution();
						_len = ls_ants[pos[i]].len();
					}
			}
			else
			{
				// Обновляем лучшее решение
				for (int i = 0; i < _n_jobs; ++i)
					if (_ants[pos[i]].len() < _len)
					{
						_solution = _ants[pos[i]].solution();
						_len = _ants[pos[i]].len();
					}
			}
			
			_tau_max = 1.0 / _rho / _len;
			_tau_min = _a * _tau_max;

			// феромоны испаряются 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_tau[i][j] = max(_tau_min, _tau[i][j] * (1 - _rho));

			// Учитываем лучшее решение
			double w = 1.0 / _len;

			#pragma omp parallel for 
			for (int j = 0; j < _n_cities; ++j)
				_tau[_solution[j]][_solution[j + 1]] = min(_tau_max, _tau[_solution[j]][_solution[j + 1]] + w);

			system_clock::time_point end = system_clock::now();
			duration <double> delta = end - start;

			fout << delta.count() << " ";
		}

		fout << "\n";
	}

public:

	/// <summary>
	/// Конструктор
	/// </summary>
	/// <param name="type"> вид применяемого алгоритма </param>
	/// <param name="params"> словарь параметров алгоритма </param>
	ACOSolver(string type, map<string, any> params) : _type(type)
	{
		// Повсеместные параметры
		_alpha = any_cast<double>(params["alpha"]);
		_beta = any_cast<double>(params["beta"]);
		_rho = any_cast<double>(params["rho"]);
		_tau0 = any_cast<double>(params["tau0"]);
		
		// Повсеместные параметры
		_n_ants = any_cast<int>(params["n_ants"]);
		_max_iter = any_cast<int>(params["max_iter"]);
		_n_jobs = any_cast<int>(params["n_jobs"]);

		// В случае элитной муравьиной системы
		if (type == "EAS" || type == "ASRank") _w = any_cast<int>(params["w"]);

		// В случае максминной муравьиной системы
		if (type == "MMAS")
		{
			_a = any_cast<double>(params["a"]);

			_local_search_type = any_cast<string>(params["local_search_type"]);

			// В случае наличия локального поиска оптимального решения
			if (_local_search_type != "None")
			{
				_local_search_tours = any_cast<string>(params["local_search_tours"]);
				_k = any_cast<int>(params["k"]);
			}
		}
	}

	/// <summary>
	/// Решаем задачу
	/// </summary>
	/// <param name="g"> граф </param>
	void solve(Graph& g)
	{
		// Этап инициализации:

		_len = INF;
		_n_cities = g.n();
		omp_set_num_threads(_n_jobs);

		_tau = vector<vector<double>>(_n_cities, vector<double>(_n_cities, _tau0));
		_eta_beta = vector<vector<double>>(_n_cities, vector<double>(_n_cities));
		_weights = vector<vector<double>>(_n_cities, vector<double>(_n_cities));

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			for (int j = 0; j < _n_cities; ++j)
				_eta_beta[i][j] = pow(g[i][j], -_beta);

		_vertices = vector<int>(_n_cities);

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			_vertices[i] = i;

		_ants = vector<Ant>(_n_ants);

		_choices = vector<vector<int>>(_n_jobs, vector<int>(_n_cities));
		_visited = vector<vector<int>>(_n_jobs, vector<int>(_n_cities));

		if (_type == "AS" || _type == "EAS") elitist_ant_system(g);
		else if (_type == "ASRank") rank_based_ant_system(g);
		else if (_type == "MMAS") max_min_ant_system(g);
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