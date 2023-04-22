#pragma once
#include <algorithm>
#include <map>
#include "TSPSolver.h"
#include "Ant.h"

/// <summary>
/// ���� ���������� ����������
/// </summary>
enum ACOType { AS, ElitistAS, RankBasedAS, MaxMinAS };

/// <summary>
/// ���������� �������� ������� 
/// ������ ������������
/// </summary>
class ACOSolver : public TSPSolver
{
public:

	/// <summary>
	/// ��� ���������
	/// </summary>
	ACOType _type;

	/// <summary>
	/// ��������������
	/// </summary>
	double _alpha, _beta, _rho, _m, _tau0, _max_iter, _e;

	/// <summary>
	/// �����������
	/// </summary>
	/// <param name="g"> ���� </param>
	/// <param name="alpha"> ����� </param>
	/// <param name="beta"> ���� </param>
	/// <param name="rho"> �� </param>
	/// <param name="m"> ����� ������� </param>
	/// <param name="tau0"> ��������� ����� ��������� </param>
	/// <param name="max_iter"> ����� ��������� </param>
	ACOSolver(ACOType type, map<string, double> params) : _type(type)
	{
		// ������������ ���������
		_alpha = params["alpha"];
		_beta = params["beta"];
		_rho = params["rho"];
		_m = params["m"];
		_tau0 = params["tau0"];
		_max_iter = params["max_iter"];

		// � ������ ������� ���������� �������
		if (type == ElitistAS) _e = params["e"];

		// � ������ �������� ���������� �������
		if (type == RankBasedAS) _e = params["w"];
	}

	/// <summary>
	/// ������ ������
	/// </summary>
	/// <param name="g"> ���� </param>
	void solve(Graph& g)
	{
		_len = INF;
		_size = g.n();

		// ������� ������������ ��������� �� ������ �����
		vector<vector<double>> tau(_size, vector<double>(_size, _tau0));

		// ������� ����������������� ����� ����� 
		// (��� �������� ��������� � ������� -_beta)
		vector<vector<double>> eta_beta(_size, vector<double>(_size));

		for (int i = 0; i < _size; ++i)
			for (int j = 0; j < _size; ++j)
				eta_beta[i][j] = pow(g[i][j], -_beta);

		// ������� ����� (�����������) ������� �� �����
		// (��. ������� ����������� ��������� ������� � ��� ��� ���� �����)
		vector<vector<double>> weights(_size, vector<double>(_size));

		// ��� ������� �����
		vector<int> vertices(_size);

		for (int i = 0; i < _size; ++i)
			vertices[i] = i;

		// ���������� �������� 
		vector<Ant> ants(_m);

		// visited[i] = 0 <=> ������� i �� ��������
		// visited[i] = 1 <=> ������� i ��������
		vector<int> visited(_size);

		// �� ���������� �������
		vector<int> choices(_size);

		for (int it = 0; it < _max_iter; ++it)
		{
			for (int i = 0; i < _size; ++i)
				for (int j = 0; j < _size; ++j)
					weights[i][j] = pow(tau[i][j], _alpha) * eta_beta[i][j];

			// ������������ ������� �����
			shuffle(vertices.begin(), vertices.end(), gen);

			// ������� ������� ������� ��������
			int pos = 0;

			// ������� ���� �������
			for (int i = 0; i < _m; ++i)
			{
				ants[i].solve(g, vertices[i], visited, choices, weights);

				for (int j = 0; j < _size; ++j)
					visited[j] = 0;

				// ��������� �������
				if (ants[i].len() < ants[pos].len()) pos = i;
			}

			// ��������� �������
			if (ants[pos].len() < _len)
			{
				_solution = ants[pos].solution();
				_len = ants[pos].len();
			}

			// �������� ���������� 
			for (int i = 0; i < _size; ++i)
				for (int j = 0; j < _size; ++j)
					tau[i][j] *= _rho;

			// � ������ �������� ���������� �������
			if (_type == RankBasedAS) sort(ants.begin(), ants.end());

			int iter = (_type == RankBasedAS) ? _e - 1 : _m;

			// ��������� ���� ��������� �� ������ 
			for (int i = 0; i < iter; ++i)
			{
				vector<int> solution = ants[i].solution();
				double w = ((_type == RankBasedAS) ? (_e - i - 1) : 1.0) / ants[i].len();

				for (int j = 0; j < _size; ++j)
					tau[solution[j]][solution[j + 1]] += w;
			}

			// � ������ ������� ��� �������� ���������� ������
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
		cout << "ACOSolver (��������� �������): ";

		for (int i = 0; i < _size; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", ����� ���� = " << _len << "\n";
	}
};