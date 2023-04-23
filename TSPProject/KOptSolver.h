#pragma once
#include "TSPSolver.h"

/// <summary>
/// ������ �������� ������� ������
/// ������������: k-opt �������� 
/// </summary>
class KOptSolver : public TSPSolver
{
	/// <summary>
	/// ��� ���������
	/// </summary>
	string _type;

	/// <summary>
	/// ������ ������ ���������� 2-opt
	/// </summary>
	void two_opt(Graph& g)
	{
		// ������� ����� ���� ��������������� �������
		vector<int> pref(_n_cities + 1);

		for (int i = 1; i <= _n_cities; ++i)
			pref[i] = pref[i - 1] + g[_solution[i - 1]][_solution[i]];

		// ������� ����� ������������� ���� ��������������� �������
		vector<int> pref_rev(_n_cities + 1);

		for (int i = 1; i <= _n_cities; ++i)
			pref_rev[i] = pref_rev[i - 1] + g[_solution[_n_cities - i + 1]][_solution[_n_cities - i]];

		for (int i = 0; i < _n_cities - 2; ++i)
			for (int j = i + 2; j < _n_cities; ++j)
			{
				// ������� ���� (i, i + 1) � (j, j + 1)

				// ������ ����� ����: 0 -> 1 -> ... -> i
				int new_len = pref[i];

				// ����� ����: i -> j
				new_len += g[_solution[i]][_solution[j]]; 

				// ������ ����� ����: j -> j - 1 -> ... -> i + 1
				new_len += pref_rev[_n_cities - i - 1] - pref_rev[_n_cities - j];

				// ����� ����: i + 1 -> j + 1
				new_len += g[_solution[i + 1]][_solution[j + 1]]; 

				// ������ ����� ����: j + 1 -> j + 2 -> n - 1 -> 0 (n)
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

public:

	/// <summary>
	/// �����������
	/// </summary>
	KOptSolver() {}

	/// <summary>
	/// ������ ������
	/// </summary>
	/// <param name="g"> ���� </param>
	/// <param name="type"> ��� ������������ ��������� </param>
	/// <param name="init_sol"> �������������� ������� </param>
	/// <param name="init_len"> ����� ��������������� ������� </param>
	void solve(Graph& g, string type, vector<int> init_sol, int init_len)
	{
		_type = type;

		_solution = init_sol;
		_len = init_len;

		_n_cities = init_sol.size() - 1;

		if (_type == "2-opt") two_opt(g);
	}

	/// <summary>
	/// ������� ������� �� �����
	/// </summary>
	void print()
	{
		cout << "TwoOptSolver (��������� �������): ";

		for (int i = 0; i < _n_cities; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", ����� ���� = " << _len << "\n";
	}
};