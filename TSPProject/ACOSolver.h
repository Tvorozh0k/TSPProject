#pragma once
#include <algorithm>
#include <map>
#include "TSPSolver.h"
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
	/// Гиперпараметры
	/// </summary>
	double _alpha, _beta, _rho, _m, _tau0, _max_iter, _e;

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
	ACOSolver(ACOType type, map<string, double> params) : _type(type)
	{
		// Повсеместные параметры
		_alpha = params["alpha"];
		_beta = params["beta"];
		_rho = params["rho"];
		_m = params["m"];
		_tau0 = params["tau0"];
		_max_iter = params["max_iter"];

		// В случае элитной муравьиной системы
		if (type == ElitistAS) _e = params["e"];

		// В случае ранговой муравьиной системы
		if (type == RankBasedAS) _e = params["w"];
	}

	/// <summary>
	/// Решаем задачу
	/// </summary>
	/// <param name="g"> граф </param>
	void solve(Graph& g)
	{
		_len = INF;
		_size = g.n();

		// матрица концентрации феромонов на ребрах графа
		vector<vector<double>> tau(_size, vector<double>(_size, _tau0));

		// матрица привлекательности ребер графа 
		// (все элементы возведены в степень -_beta)
		vector<vector<double>> eta_beta(_size, vector<double>(_size));

		for (int i = 0; i < _size; ++i)
			for (int j = 0; j < _size; ++j)
				eta_beta[i][j] = pow(g[i][j], -_beta);

		// матрица весов (значимостей) каждого из ребер
		// (см. формула вероятности попадания муравья в тот или иной город)
		vector<vector<double>> weights(_size, vector<double>(_size));

		// все вершины графа
		vector<int> vertices(_size);

		for (int i = 0; i < _size; ++i)
			vertices[i] = i;

		// генерируем муравьев 
		vector<Ant> ants(_m);

		// visited[i] = 0 <=> вершина i НЕ посещена
		// visited[i] = 1 <=> вершина i посещена
		vector<int> visited(_size);

		// не посещенные вершины
		vector<int> choices(_size);

		for (int it = 0; it < _max_iter; ++it)
		{
			for (int i = 0; i < _size; ++i)
				for (int j = 0; j < _size; ++j)
					weights[i][j] = pow(tau[i][j], _alpha) * eta_beta[i][j];

			// перемешиваем вершины графа
			shuffle(vertices.begin(), vertices.end(), gen);

			// позиция лучшего решения итерации
			int pos = 0;

			// муравей ищет решение
			for (int i = 0; i < _m; ++i)
			{
				ants[i].solve(g, vertices[i], visited, choices, weights);

				for (int j = 0; j < _size; ++j)
					visited[j] = 0;

				// сохраняем решение
				if (ants[i].len() < ants[pos].len()) pos = i;
			}

			// сохраняем решение
			if (ants[pos].len() < _len)
			{
				_solution = ants[pos].solution();
				_len = ants[pos].len();
			}

			// феромоны испаряются 
			for (int i = 0; i < _size; ++i)
				for (int j = 0; j < _size; ++j)
					tau[i][j] *= _rho;

			// В случае ранговой муравьиной системы
			if (_type == RankBasedAS) sort(ants.begin(), ants.end());

			int iter = (_type == RankBasedAS) ? _e - 1 : _m;

			// оставляем след феромонов на ребрах 
			for (int i = 0; i < iter; ++i)
			{
				vector<int> solution = ants[i].solution();
				double w = ((_type == RankBasedAS) ? (_e - i - 1) : 1.0) / ants[i].len();

				for (int j = 0; j < _size; ++j)
					tau[solution[j]][solution[j + 1]] += w;
			}

			// В случае элитной или ранговой муравьиной систем
			if (_type == ElitistAS || _type == RankBasedAS)
			{
				double w = _e / _len;

				for (int j = 0; j < _size; ++j)
					tau[_solution[j]][_solution[j + 1]] += w;
			}
		}
	}

	void print()
	{
		cout << "ACOSolver (найденное решение): ";

		for (int i = 0; i < _size; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", длина пути = " << _len << "\n";
	}
};