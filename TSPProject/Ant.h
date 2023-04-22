#pragma once
#include "TSPSolver.h"

/// <summary>
/// Класс муравей
/// </summary>
class Ant : public TSPSolver
{
	/// <summary>
	/// Муравей выбирает следуюшую вершину
	/// </summary>
	/// <param name="visited"> массив с флагом посещена/не посещена для каждой вершмны </param>
	/// <param name="tau_alpha"> концентрация феромонов, уже возведенная в степень alpha </param>
	/// <param name="eta"> привлекательность ребер в степени beta </param>
	/// <returns> следующая вершина </returns>
	int next(vector<int>& visited, vector<int>& choices, vector<double>& weights)
	{
		int pos = 0;
		double sum = 0;

		for (int i = 0; i < _size; ++i)
			if (!visited[i])
			{
				choices[pos++] = i;
				sum += weights[i];
			}

		// генерируем случайную величину 
		uniform_real_distribution <double> interval_double(0, sum);
		double rnd = interval_double(gen);

		for (int i = 0; i < pos; ++i)
			if (rnd > weights[choices[i]]) rnd -= weights[choices[i]];
			else return choices[i];
	}

public:

	Ant() {}

	/// <summary>
	/// Муравей прокладывает свой путь
	/// </summary>
	/// <param name="tau_alpha"> матрица концентраций феромонов (в степени alpha) </param>
	/// <param name="eta"> матрица привлекательности ребер (в степени beta) </param>
	void solve(Graph& g, int s, vector<int>& visited, vector<int>& choices, vector<vector<double>>& weights)
	{
		_size = g.n();
		_solution.resize(_size + 1);

		// путь: s -> ...
		_solution[0] = s; _len = 0;
		visited[s] = 1;

		for (int i = 1; i < _size; ++i)
		{
			int from = _solution[i - 1];

			// выбираем следующую вершину
			int to = next(visited, choices, weights[from]);

			// путь: ... -> from -> to -> ...
			_solution[i] = to; _len += g[from][to];
			visited[to] = 1;
		}

		// путь: s -> ... -> t -> s
		int t = _solution[_size - 1];
		_solution[_size] = s; _len += g[t][s];
	}

	/// <summary>
	/// Выводим решение на экран
	/// </summary>
	void print()
	{
		cout << "Ant (найденное решение): ";

		for (int i = 0; i < _size; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", длина пути = " << _len << "\n";
	}

	/// <summary>
	/// Сравнение муравьев посредством длин 
	/// найденных решений
	/// </summary>
	/// <param name="ant"> муравей </param>
	/// <returns> значение _len МЕНЬШЕ ant.len() </returns>
	bool operator < (Ant& ant)
	{
		return _len < ant.len();
	}
};