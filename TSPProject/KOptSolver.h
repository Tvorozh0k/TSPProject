#pragma once
#include "TSPSolver.h"

/// <summary>
/// Жадный алгоритм решения задачи
/// коммивояжера: k-opt алгоритм 
/// </summary>
class KOptSolver : public TSPSolver
{
	/// <summary>
	/// Вид алгоритма
	/// </summary>
	string _type;

	/// <summary>
	/// Алгоритм 2-opt
	/// </summary>
	/// <param name="g"> граф </param>
	void two_opt(Graph& g)
	{
		// префикс длины пути первоначального решения
		vector<int> pref(_n_cities + 1);

		for (int i = 1; i <= _n_cities; ++i)
			pref[i] = pref[i - 1] + g[_solution[i - 1]][_solution[i]];

		// префикс длины перевернутого пути первоначального решения
		vector<int> pref_rev(_n_cities + 1);

		for (int i = 1; i <= _n_cities; ++i)
			pref_rev[i] = pref_rev[i - 1] + g[_solution[_n_cities - i + 1]][_solution[_n_cities - i]];
		
		for (int i = 0; i < _n_cities - 2; ++i)
			for (int j = i + 2; j < _n_cities; ++j)
			{
				// удаляем дуги (i, i + 1) и (j, j + 1)

				// первая часть пути: 0 -> 1 -> ... -> i
				int new_len = pref[i];

				// новая дуга: i -> j
				new_len += g[_solution[i]][_solution[j]]; 

				// вторая часть пути: j -> j - 1 -> ... -> i + 1
				new_len += pref_rev[_n_cities - i - 1] - pref_rev[_n_cities - j];

				// новая дуга: i + 1 -> j + 1
				new_len += g[_solution[i + 1]][_solution[j + 1]]; 

				// третья часть пути: j + 1 -> j + 2 -> n - 1 -> 0 (n)
				new_len += _len - pref[j + 1]; 

				if (new_len < _len)
				{
					reverse(_solution.begin() + i + 1, _solution.begin() + j + 1);
					_len = new_len;

					two_opt(g);
					return;
				}
			}
	}

	/// <summary>
	/// Алгоритм 2.5-opt
	/// </summary>
	/// <param name="g"> граф </param>
	void two_half_opt(Graph& g)
	{
		for (int i = 1; i < _n_cities - 2; ++i)
		{
			// длина нового пути
			int new_len = _len;

			// (n - 1) -> 0 -> 1 -> ... -> i -> i + 1 -> ...
			new_len -= g[_solution[_n_cities - 1]][_solution[0]];
			new_len -= g[_solution[0]][_solution[1]];
			new_len -= g[_solution[i]][_solution[i + 1]];

			// i -> 0 -> i + 1 -> ... -> (n - 1) -> 1 -> ...
			new_len += g[_solution[i]][_solution[0]];
			new_len += g[_solution[0]][_solution[i + 1]];
			new_len += g[_solution[_n_cities - 1]][_solution[1]];

			if (new_len < _len)
			{
				rotate(_solution.begin() + 1, _solution.begin() + i + 1, _solution.begin() + _n_cities);
				_len = new_len;

				two_half_opt(g);
				return;
			}
		}

		for (int i = 1; i < _n_cities; ++i)
			for (int j = 1; j < _n_cities - 2; ++j)
			{
				// длина нового пути
				int new_len = _len;

				// (i - 1) -> i -> (i + 1) -> ... -> (i + j) -> (i + j + 1) -> ...
				new_len -= g[_solution[i - 1]][_solution[i]];
				new_len -= g[_solution[i]][_solution[i + 1]];
				new_len -= g[_solution[(i + j) % _n_cities]][_solution[(i + j + 1) % _n_cities]];

				// (i + j) -> i -> (i + j + 1) -> ... -> (i - 1) -> (i + 1) -> ...
				new_len += g[_solution[(i + j) % _n_cities]][_solution[i]];
				new_len += g[_solution[i]][_solution[(i + j + 1) % _n_cities]];
				new_len += g[_solution[i - 1]][_solution[i + 1]];

				if (new_len < _len)
				{
					int pos1 = min(i, (i + j) % _n_cities + (i + j) / _n_cities);
					int pos2 = max(i, (i + j) % _n_cities + (i + j) / _n_cities);

					int pos3 = i + 1 - (i + j) / _n_cities;

					rotate(_solution.begin() + pos1, _solution.begin() + pos3, _solution.begin() + pos2 + 1);
					_len = new_len;

					two_half_opt(g);
					return;
				}
			}
	}

	/// <summary>
	/// Алгоритм 3-opt
	/// </summary>
	/// <param name="g"> граф </param>
	void three_opt(Graph& g)
	{
		for (int i = 0; i < _n_cities - 2; ++i)
			for (int j = i + 1; j < _n_cities - 1; ++j)
				for (int k = j + 1; k < _n_cities; ++k)
				{
					// длина нового пути
					int new_len = _len;

					// ... -> x -> (x + 1) -> ...
					new_len -= g[_solution[i]][_solution[i + 1]];
					new_len -= g[_solution[j]][_solution[j + 1]];
					new_len -= g[_solution[k]][_solution[k + 1]];

					// i -> j + 1 -> ... -> k -> i + 1 -> ... -> j -> k + 1 -> ... 
					new_len += g[_solution[i]][_solution[j + 1]];
					new_len += g[_solution[k]][_solution[i + 1]];
					new_len += g[_solution[j]][_solution[k + 1]];

					if (new_len < _len)
					{
						rotate(_solution.begin() + i + 1, _solution.begin() + j + 1, _solution.begin() + k + 1);
						_len = new_len;

						three_opt(g);
						return;
					}
				}
	}

public:

	/// <summary>
	/// Конструктор
	/// </summary>
	KOptSolver() {}

	/// <summary>
	/// Решаем задачу
	/// </summary>
	/// <param name="g"> граф </param>
	/// <param name="type"> вид применяемого алгоритма </param>
	/// <param name="init_sol"> первоначальное решение </param>
	/// <param name="init_len"> длина первоначального решения </param>
	void solve(Graph& g, string type, vector<int> init_sol, int init_len)
	{
		_type = type;

		_solution = init_sol;
		_len = init_len;

		_n_cities = init_sol.size() - 1;

		if (_type == "2-opt") two_opt(g);
		else if (_type == "2.5-opt") two_half_opt(g);
		else if (_type == "3-opt") three_opt(g);
	}

	/// <summary>
	/// Выводим решение на экран
	/// </summary>
	void print()
	{
		cout << "TwoOptSolver (найденное решение): ";

		for (int i = 0; i < _n_cities; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", длина пути = " << _len << "\n";
	}
};