#pragma once
#include "TSPSolver.h"

/// <summary>
/// ����� �������
/// </summary>
class Ant : public TSPSolver
{
	/// <summary>
	/// ������� �������� ��������� �������
	/// </summary>
	/// <param name="visited"> ������ � ������ ��������/�� �������� ��� ������ ������� </param>
	/// <param name="tau_alpha"> ������������ ���������, ��� ����������� � ������� alpha </param>
	/// <param name="eta"> ����������������� ����� � ������� beta </param>
	/// <returns> ��������� ������� </returns>
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

		// ���������� ��������� �������� 
		uniform_real_distribution <double> interval_double(0, sum);
		double rnd = interval_double(gen);

		for (int i = 0; i < pos; ++i)
			if (rnd > weights[choices[i]]) rnd -= weights[choices[i]];
			else return choices[i];
	}

public:

	Ant() {}

	/// <summary>
	/// ������� ������������ ���� ����
	/// </summary>
	/// <param name="tau_alpha"> ������� ������������ ��������� (� ������� alpha) </param>
	/// <param name="eta"> ������� ����������������� ����� (� ������� beta) </param>
	void solve(Graph& g, int s, vector<int>& visited, vector<int>& choices, vector<vector<double>>& weights)
	{
		_size = g.n();
		_solution.resize(_size + 1);

		// ����: s -> ...
		_solution[0] = s; _len = 0;
		visited[s] = 1;

		for (int i = 1; i < _size; ++i)
		{
			int from = _solution[i - 1];

			// �������� ��������� �������
			int to = next(visited, choices, weights[from]);

			// ����: ... -> from -> to -> ...
			_solution[i] = to; _len += g[from][to];
			visited[to] = 1;
		}

		// ����: s -> ... -> t -> s
		int t = _solution[_size - 1];
		_solution[_size] = s; _len += g[t][s];
	}

	/// <summary>
	/// ������� ������� �� �����
	/// </summary>
	void print()
	{
		cout << "Ant (��������� �������): ";

		for (int i = 0; i < _size; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", ����� ���� = " << _len << "\n";
	}

	/// <summary>
	/// ��������� �������� ����������� ���� 
	/// ��������� �������
	/// </summary>
	/// <param name="ant"> ������� </param>
	/// <returns> �������� _len ������ ant.len() </returns>
	bool operator < (Ant& ant)
	{
		return _len < ant.len();
	}
};