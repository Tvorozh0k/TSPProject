#pragma once
#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include "Graph.h"

/// <summary>
/// Формальная бесконечность
/// </summary>
const int INF = 1e9 + 7;

/// <summary>
/// Генератор псевдослучайных чисел
/// </summary>
mt19937 gen(time(NULL));

/// <summary>
/// Абстрактный класс - алгоритм,
/// разрешающий задачу коммивояжера
/// </summary>
class TSPSolver
{
protected:

	/// <summary>
	/// Решение
	/// </summary>
	vector<int> _solution;

	/// <summary>
	/// Длина (стоимость) решения
	/// </summary>
	int _len;

	/// <summary>
	/// Число вершин в графе
	/// </summary>
	int _size;

public:

	/// <summary>
	/// Выводим решение на экран
	/// </summary>
	virtual void print() = 0;

	/// <summary>
	/// Само решение
	/// </summary>
	/// <returns> само решение (путь) </returns>
	vector<int> solution() { return _solution; }

	/// <summary>
	/// Длина (стоимость) решения
	/// </summary>
	/// <returns> длина (стоимость) решения </returns>
	int len() { return _len; }
};