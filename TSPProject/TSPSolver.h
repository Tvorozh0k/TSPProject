#pragma once
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>

#include "Graph.h"

/// <summary>
/// ���������� �������������
/// </summary>
const int INF = 1e9 + 7;

/// <summary>
/// ��������� ��������������� �����
/// </summary>
mt19937 gen(time(NULL));

/// <summary>
/// ����������� ����� - ��������,
/// ����������� ������ ������������
/// </summary>
class TSPSolver
{
protected:

	/// <summary>
	/// �������
	/// </summary>
	vector<int> _solution;

	/// <summary>
	/// ����� (���������) �������
	/// </summary>
	int _len;

	/// <summary>
	/// ����� ������ � �����
	/// </summary>
	int _n_cities;

public:

	/// <summary>
	/// ������� ������� �� �����
	/// </summary>
	virtual void print() = 0;

	/// <summary>
	/// ���� �������
	/// </summary>
	/// <returns> ���� ������� (����) </returns>
	vector<int> solution() { return _solution; }

	/// <summary>
	/// ����� (���������) �������
	/// </summary>
	/// <returns> ����� (���������) ������� </returns>
	int len() { return _len; }

	/// <summary>
	/// ��������� ���������� ����������� ���� 
	/// ��������� �������
	/// </summary>
	/// <param name="solver"> �������� </param>
	/// <returns> �������� _len ������ solver.len()  </returns>
	bool operator < (TSPSolver& solver)
	{
		return _len < solver.len();
	}
};