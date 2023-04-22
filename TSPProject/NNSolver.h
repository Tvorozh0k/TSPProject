#pragma once
#include "TSPSolver.h"

/// <summary>
/// Жадный алгоритм решения задачи
/// коммивояжера: метод ближайшего соседа
/// </summary>
class NNSolver : public TSPSolver
{
public:

	/// <summary>
	/// Решаем задачу
	/// </summary>
	void solve(Graph& g)
	{
		_size = g.n();

		// visited[i] = 0 <=> вершина i НЕ посещена
		// visited[i] = 1 <=> вершина i посещена
		vector <int> visited(_size);

		// генерируем первую вершину случайно
		uniform_int_distribution <int> interval_int(0, _size - 1);
		int cur = interval_int(gen);

		_solution.resize(_size + 1);

		// путь: cur -> ... 
		_solution[0] = cur; _len = 0;
		visited[cur] = 1;

		// достраиваем путь до конца
		for (int i = 1; i < _size; ++i)
		{
			// ищем вершину next
			int next = cur, dist = INF;

			// проходимся по всем НЕ посещенным вершинам j 
			// и ещем наиболее близкую к вершине cur
			for (int j = 0; j < _size; ++j)
				if (!visited[j] && g[cur][j] < dist)
				{
					next = j;
					dist = g[cur][j];
				}

			// путь: ... -> cur -> next -> ...
			_solution[i] = next; _len += dist;
			visited[next] = 1;

			cur = next;
		}

		// путь: s -> ... -> t -> s
		int s = _solution[0], t = _solution[_size - 1];
		_solution[_size] = s; _len += g[t][s];
	}

	/// <summary>
	/// Выводим решение на экран
	/// </summary>
	void print()
	{
		cout << "NNSolver (найденное решение): ";

		for (int i = 0; i < _size; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", длина пути = " << _len << "\n";
	}
};