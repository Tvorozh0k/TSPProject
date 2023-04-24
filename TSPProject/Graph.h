#pragma once
#include <fstream>
#include <cassert>
#include <vector>
using namespace std;

/// <summary>
/// ����� ����
/// </summary>
class Graph
{
	/// <summary>
	/// ����� ������ � �����
	/// </summary>
	int _N;

	/// <summary>
	/// ������� ��������� �����
	/// </summary>
	vector<vector<int>> _adjMat;

public:

	Graph() {}

	/// <summary>
	/// �����������: ��������� ������� ���������
	/// �� �����
	/// </summary>
	/// <param name="filePath"> ���� � ����� </param>
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
	/// ����� ������ � �����
	/// </summary>
	/// <returns> ����� ������ � ����� </returns>
	int n() { return _N; }

	/// <summary>
	/// ���������� �������� ���������
	/// </summary>
	/// <param name="i"> ������ </param>
	/// <returns> ������ ������� ���������, ��������������� ������� � �������� i </returns>
	vector<int>& operator [] (int i) { return _adjMat[i]; }
};