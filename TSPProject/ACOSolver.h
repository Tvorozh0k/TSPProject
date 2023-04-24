#pragma once
#include <string>
#include <map>
#include <any>

#include "KOptSolver.h"
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
	/// ��������������:
	/// 
	/// alpha - ���������� ���������,
	/// beta - ���������� ����������������� �����,
	/// rho - ����������� ��������� ���������,
	/// tau0 - ����������� ���-�� ��������� �� ������,
	/// tau_min - ������ ������� ���������� ����� ���������,
	/// tau_max - ������� ������� ���������� ����� ���������,
	/// ������ tau_min = a * tau_max 
	/// </summary>
	double _alpha, _beta, _rho, _tau0, _tau_min, _tau_max, _a;

	/// <summary>
	/// ��������������:
	/// 
	/// n_ants - ����� ��������,
	/// max_iter - ����� ���������,
	/// e - ����������� ���������� ������� ����
	/// </summary>
	int _n_ants, _max_iter, _e;

	/// <summary>
	/// ��������������:
	/// 
	/// bst_tour_type - ��� ������� ����,
	/// ������������ �������� �� ������
	/// local_search_type - �������� ��������� 
	/// ����������� �������
	/// local_search_tours - ����� ���� ��������������
	/// </summary>
	string _bst_tour_type, _local_search_type, _local_search_tours;

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
	ACOSolver(ACOType type, map<string, any> params) : _type(type)
	{
		// ������������ ���������
		_alpha = any_cast<double>(params["alpha"]);
		_beta = any_cast<double>(params["beta"]);
		_rho = any_cast<double>(params["rho"]);
		_tau0 = any_cast<double>(params["tau0"]);

		_n_ants = any_cast<int>(params["n_ants"]);
		_max_iter = any_cast<int>(params["max_iter"]);

		_local_search_type = any_cast<string>(params["local_search_type"]);

		// � ������ ������� ���������� ������ ������������ �������
		if (_local_search_type != "None") _local_search_tours = any_cast<string>(params["local_search_tours"]);
		else  _local_search_tours = "None";

		// � ������ ������� ���������� �������
		if (type == ElitistAS)
		{
			_e = any_cast<int>(params["e"]);
			_bst_tour_type = "best-so-far";
		}

		// � ������ �������� ���������� �������
		if (type == RankBasedAS)
		{
			_e = any_cast<int>(params["w"]);
			_bst_tour_type = "best-so-far";
		}

		// � ������ ���������� ���������� �������
		if (type == MaxMinAS)
		{
			_e = 1;
			_a = any_cast<double>(params["a"]);
			_bst_tour_type = any_cast<string>(params["bst_tour_type"]);
		}
	}

	/// <summary>
	/// ������ ������
	/// </summary>
	/// <param name="g"> ���� </param>
	void solve(Graph& g)
	{
		_len = INF;
		_n_cities = g.n();

		// � ������ ���������� ���������� �������
		if (_type == MaxMinAS)
		{
			_tau_max = _tau0;
			_tau_min = _a * _tau_max;
		}
		else
		{
			_tau_min = 0;
			_tau_max = INF;
		}

		// ������� ������������ ��������� �� ������ �����
		vector<vector<double>> tau(_n_cities, vector<double>(_n_cities, _tau0));

		// ������� ����������������� ����� ����� 
		// (��� �������� ��������� � ������� -_beta)
		vector<vector<double>> eta_beta(_n_cities, vector<double>(_n_cities));

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			for (int j = 0; j < _n_cities; ++j)
				eta_beta[i][j] = pow(g[i][j], -_beta);

		// ������� ����� (�����������) ������� �� �����
		// (��. ������� ����������� ��������� ������� � ��� ��� ���� �����)
		vector<vector<double>> weights(_n_cities, vector<double>(_n_cities));

		// ��� ������� �����
		vector<int> vertices(_n_cities);

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			vertices[i] = i;

		// ���������� �������� 
		vector<Ant> ants(_n_ants);

		// visited[i] = 0 <=> ������� i �� ��������
		// visited[i] = 1 <=> ������� i ��������
		vector<vector<int>> visited(GridThreadsNum, vector<int>(_n_cities));

		// �� ���������� �������
		vector<vector<int>> choices(GridThreadsNum, vector<int>(_n_cities));

		for (int it = 0; it < _max_iter; ++it)
		{
			#pragma omp parallel for
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					weights[i][j] = pow(tau[i][j], _alpha) * eta_beta[i][j];

			// ������������ ������� �����
			shuffle(vertices.begin(), vertices.end(), gen);

			// ������ �� �������� �������
			vector<vector<int>> iter_best_sol(GridThreadsNum);
			vector<int> iter_best_len(GridThreadsNum, INF);

			// ������� ���� �������	
			#pragma omp parallel for schedule(dynamic) 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				ants[i].solve(g, vertices[i], visited[thread_num], choices[thread_num], weights);

				for (int j = 0; j < _n_cities; ++j)
					visited[thread_num][j] = 0;

				// � ������ ������� ���������� ������ ������������ �������
				if (_local_search_tours == "all")
				{
					KOptSolver ls_ant;
					ls_ant.solve(g, _local_search_type, ants[i].solution(), ants[i].len());

					if (ls_ant.len() < iter_best_len[thread_num])
					{
						iter_best_sol[thread_num] = ls_ant.solution();
						iter_best_len[thread_num] = ls_ant.len();
					}
				}
				else if (ants[i].len() < iter_best_len[thread_num])
				{
					iter_best_sol[thread_num] = ants[i].solution();
					iter_best_len[thread_num] = ants[i].len();
				}
			}

			// ������ ������� �� ������� 
			int pos_min = min_element(iter_best_len.begin(), iter_best_len.end()) - iter_best_len.begin();

			// � ������ ������� ���������� ������ ������������ �������
			if (_local_search_tours == "best_sol")
			{
				KOptSolver ls_best;
				ls_best.solve(g, _local_search_type, iter_best_sol[pos_min], iter_best_len[pos_min]);

				iter_best_sol[pos_min] = ls_best.solution();
				iter_best_len[pos_min] = ls_best.len();
			}

			// ��������� �������
			if (iter_best_len[pos_min] < _len)
			{
				_solution = iter_best_sol[pos_min];
				_len = iter_best_len[pos_min];

				// � ������ ���������� ���������� �������
				if (_type == MaxMinAS)
				{
					_tau_max = 1.0 / _rho / _len;
					_tau_min = _a * _tau_max;
				}
			}

			// �������� ���������� 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					tau[i][j] = max(_tau_min, tau[i][j] * (1 - _rho));	

			// � ������ �������� ���������� �������
			if (_type == RankBasedAS) sort(ants.begin(), ants.end());

			int iter;

			if (_type == RankBasedAS || _type == MaxMinAS) iter = _e - 1;
			else iter = _n_ants;

			// ��������� ���� ��������� �� ������ 
			for (int i = 0; i < iter; ++i)
			{
				vector<int> solution = ants[i].solution();
				double w = ((_type == RankBasedAS) ? (_e - i - 1.0) : 1.0) / ants[i].len();

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					tau[solution[j]][solution[j + 1]] += w;
			}

			// � ������ �� ������������ ���������� �������
			if (_type != AS)
			{
				vector<int> solution;
				double w;

				if (_bst_tour_type == "best-so-far")
				{
					solution = _solution;
					w = (double)_e / _len;
				}
				else 
				{
					solution = iter_best_sol[pos_min];
					w = (double)_e / iter_best_len[pos_min];
				}

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					tau[solution[j]][solution[j + 1]] = min(_tau_max, tau[solution[j]][solution[j + 1]] + w);
			}
		}
	}

	/// <summary>
	/// ������� ������� �� �����
	/// </summary>
	void print()
	{
		cout << "ACOSolver (��������� �������): ";

		for (int i = 0; i < _n_cities; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", ����� ���� = " << _len << "\n";
	}
};