#pragma once
#include <fstream>
#include <cassert>
#include <vector>
using namespace std;

/// <summary>
/// Класс граф
/// </summary>
class Graph
{
	/// <summary>
	/// Число вершин в графе
	/// </summary>
	int _N;

	/// <summary>
	/// Матрица смежности графа
	/// </summary>
	vector<vector<int>> _adjMat;

public:

	Graph() {}

	/// <summary>
	/// Конструктор: считываем матрицу смежности
	/// из файла
	/// </summary>
	/// <param name="filePath"> путь к файлу </param>
	Graph(string filePath)
	{
		ifstream fin(filePath);
		assert(fin.is_open());

		fin >> _N;

		_adjMat.resize(_N);

		for (int i = 0; i < _N; ++i)
			_adjMat[i].resize(_N);

		for (int i = 0; i < _N; ++i)
			for (int j = 0; j < _N; ++j)
				fin >> _adjMat[i][j];

		fin.close();
	}

	/// <summary>
	/// Число вершин в графе
	/// </summary>
	/// <returns> число вершин в графе </returns>
	int n() { return _N; }

	/// <summary>
	/// Перегрузка операции адресации
	/// </summary>
	/// <param name="i"> индекс </param>
	/// <returns> строка матрицы смежности, соответствующая вершине с индексом i </returns>
	vector<int>& operator [] (int i) { return _adjMat[i]; }
};