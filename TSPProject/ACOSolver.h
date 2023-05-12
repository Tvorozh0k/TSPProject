#pragma once
#include <string>
#include <map>
#include <any>

#include <chrono>
using namespace chrono;

#include "KOptSolver.h"
#include "Ant.h"

ofstream fout("iter_time.txt");

/// <summary>
/// ���������� �������� ������� 
/// ������ ������������
/// </summary>
class ACOSolver : public TSPSolver
{
	/// <summary>
	/// ��� ���������:
	/// <para> - AS (Ant System) </para> 
	/// <para> - EAS (Elitist Ant System) </para>
	/// <para> - ASRank (Rank-based Ant System) </para>
	/// <para> - MMAS (Max-Min Ant System) </para>
	/// </summary>
	string _type;

	/// <summary>
	/// ������������ ��������������:
	/// <para> - alpha (���������� ���������) </para>
	/// <para> - beta (���������� ����������������� �����) </para>
	/// <para> - rho (����������� ��������� ���������) </para>
	/// <para> - tau0 (����������� ���-�� ��������� �� ������) </para>
	/// </summary>
	double _alpha, _beta, _rho, _tau0;

	/// <summary>
	/// ������������ ��������������:
	/// <para> - n_ants (����� ��������) </para>
	/// <para> - max_iter (����� ���������) </para>
	/// <para> - n_jobs (����� ������������ �������) </para>
	/// </summary>
	int _n_ants, _max_iter, _n_jobs;

	/// <summary>
	/// ������������� w:
	/// <para> - EAS (���������� ������� ����) </para>
	/// <para> - ASRank (����� ������ ��������, �����������
	/// ��� ����������������� ��������� �� ������) </para>
	/// </summary>
	int _w;

	/// <summary>
	/// ���������� ���������� ������� (MMAS):
	/// <para> - tau_min (������ ������� ���������� ����� ���������) </para>
	/// <para> - tau_max (������� ������� ���������� ����� ���������) </para>
	/// <para> - a (��������� tau_max / tau_min) </para>
	/// </summary>
	double _tau_min, _tau_max, _a;

	/// <summary>
	/// ��������� ����������� �������, ��������� ���������:
	/// <para> - None (��� ��������� �����������) </para>
	/// <para> - 2-opt </para>
	/// <para> - 2.5-opt </para>
	/// <para> - 3-opt </para>
	/// </summary>
	string _local_search_type;

	/// <summary>
	/// ��������� ����������� �������, ��������� ���������:
	/// <para> - k-best (k ������ ��������) </para>
	/// <para> - k-random (k ��������� ��������) </para>
	/// </summary>
	string _local_search_tours;

	/// <summary>
	/// ������������� k:
	/// ����� ��������, ������� ������� ����� ��������
	/// ��������������
	/// </summary>
	int _k;

	/// <summary>
	/// ��������������� ������� ��� ������ ���������:
	/// <para> - tau (������� ������������ ��������� �� ������) </para>
	/// <para> - eta_beta (������� ����������������� �����,
	/// ����������� � ������� -beta) </para>
	/// <para> - weights (������� ����������� �����, �������� + 
	/// �����������������) </para>
	/// </summary>
	vector<vector<double>> _tau, _eta_beta, _weights;

	/// <summary>
	/// ��������������� ������ ��� ������ ���������:
	/// <para> - vertices (������ ������ �����) </para> 
	/// </summary>
	vector<int> _vertices;

	/// <summary>
	/// ��������������� ������ ��� ������ ���������:
	/// <para> - ants (������ ��������) </para>
	/// </summary>
	vector<Ant> _ants;
	
	/// <summary>
	/// ��������������� ������� ��� ������ ���������:
	/// <para> - choices (������ ������������ �������� �������) </para>
	/// <para> - visited (������ � ������� ������� / �� ������� �� �������� �����) </para>
	/// </summary>
	vector<vector<int>> _choices, _visited;

	/// <summary>
	/// (Elitist) Ant System
	/// </summary>
	/// <param name="g"> ���� </param>
	void elitist_ant_system(Graph& g)
	{
		for (int it = 0; it < _max_iter; ++it)
		{
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_weights[i][j] = pow(_tau[i][j], _alpha) * _eta_beta[i][j];

			// ������������ ������� �����
			shuffle(_vertices.begin(), _vertices.end(), gen);

			// ������ �� �������� �������
			vector<int> pos(_n_jobs);

			// ������� ���� �������	
			#pragma omp parallel for 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				_ants[i].solve(g, _vertices[i], _visited[thread_num], _choices[thread_num], _weights);

				for (int j = 0; j < _n_cities; ++j)
					_visited[thread_num][j] = 0;

				if (_ants[i] < _ants[pos[thread_num]]) pos[thread_num] = i;
			}

			for (int i = 0; i < _n_jobs; ++i)
				if (_ants[pos[i]].len() < _len)
				{
					_solution = _ants[pos[i]].solution();
					_len = _ants[pos[i]].len();
				}

			// �������� ���������� 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_tau[i][j] *= (1 - _rho); 

			// ��������� ���� ��������� �� ������ 
			for (int i = 0; i < _n_ants; ++i)
			{
				vector<int> solution = _ants[i].solution();
				double w = 1.0 / _ants[i].len();

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					_tau[solution[j]][solution[j + 1]] += w;
			}

			// � ������ ������� ���������� �������
			if (_type == "EAS")
			{
				double w = (double)_w / _len;

				#pragma omp parallel for
				for (int j = 0; j < _n_cities; ++j)
					_tau[_solution[j]][_solution[j + 1]] += w;
			}
		}
	}

	/// <summary>
	/// Rank-Based Ant System
	/// </summary>
	/// <param name="g"> ���� </param>
	void rank_based_ant_system(Graph& g)
	{
		for (int it = 0; it < _max_iter; ++it)
		{
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_weights[i][j] = pow(_tau[i][j], _alpha) * _eta_beta[i][j];

			// ������������ ������� �����
			shuffle(_vertices.begin(), _vertices.end(), gen);

			// ������ �� �������� �������
			vector<int> pos(_n_jobs);

			// ������� ���� �������	
			#pragma omp parallel for 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				_ants[i].solve(g, _vertices[i], _visited[thread_num], _choices[thread_num], _weights);

				for (int j = 0; j < _n_cities; ++j)
					_visited[thread_num][j] = 0;

				if (_ants[i] < _ants[pos[thread_num]]) pos[thread_num] = i;
			}

			for (int i = 0; i < _n_jobs; ++i)
				if (_ants[pos[i]].len() < _len)
				{
					_solution = _ants[pos[i]].solution();
					_len = _ants[pos[i]].len();
				}

			// �������� ���������� 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_tau[i][j] *= (1 - _rho);

			// ��������� �������� �� ����� �������
			sort(_ants.begin(), _ants.end());

			// ��������� ���� ��������� �� ������ 
			for (int i = 0; i < _w - 1; ++i)
			{
				vector<int> solution = _ants[i].solution();
				double w = (_w - i - 1.0) / _ants[i].len();

				#pragma omp parallel for 
				for (int j = 0; j < _n_cities; ++j)
					_tau[solution[j]][solution[j + 1]] += w;
			}

			// ��������� ������ �������
			double w = (double)_w / _len;

			#pragma omp parallel for
			for (int j = 0; j < _n_cities; ++j)
				_tau[_solution[j]][_solution[j + 1]] += w;
		}
	}

	/// <summary>
	/// Max-Min Ant System
	/// </summary>
	/// <param name="g"> ���� </param>
	void max_min_ant_system(Graph& g)
	{
		_tau_max = _tau0;
		_tau_min = _a * _tau_max;

		for (int it = 0; it < _max_iter; ++it)
		{
			system_clock::time_point start = system_clock::now();

			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_weights[i][j] = pow(_tau[i][j], _alpha) * _eta_beta[i][j];

			// ������������ ������� �����
			shuffle(_vertices.begin(), _vertices.end(), gen);

			// ��������� ����������� �������
			vector<KOptSolver> ls_ants(_k);

			// ������ �� �������� �������
			vector<int> pos(_n_jobs);

			// ������� ���� �������	
			#pragma omp parallel for 
			for (int i = 0; i < _n_ants; ++i)
			{
				int thread_num = omp_get_thread_num();

				_ants[i].solve(g, _vertices[i], _visited[thread_num], _choices[thread_num], _weights);

				for (int j = 0; j < _n_cities; ++j)
					_visited[thread_num][j] = 0;

				if (_local_search_type == "None" && _ants[i] < _ants[pos[thread_num]]) 
					pos[thread_num] = i;
			}

			// ���� ���� ��������� ����������� �������
			if (_local_search_type != "None")
			{
				if (_local_search_tours == "k-best") sort(_ants.begin(), _ants.end());
				else shuffle(_ants.begin(), _ants.end(), gen);

				#pragma omp parallel for 
				for (int i = 0; i < _k; ++i)
				{
					int thread_num = omp_get_thread_num();

					ls_ants[i].solve(g, _local_search_type, _ants[i].solution(), _ants[i].len());

					if (ls_ants[i] < ls_ants[pos[thread_num]]) pos[thread_num] = i;
				}

				// ��������� ������ �������
				for (int i = 0; i < _n_jobs; ++i)
					if (ls_ants[pos[i]].len() < _len)
					{
						_solution = ls_ants[pos[i]].solution();
						_len = ls_ants[pos[i]].len();
					}
			}
			else
			{
				// ��������� ������ �������
				for (int i = 0; i < _n_jobs; ++i)
					if (_ants[pos[i]].len() < _len)
					{
						_solution = _ants[pos[i]].solution();
						_len = _ants[pos[i]].len();
					}
			}
			
			_tau_max = 1.0 / _rho / _len;
			_tau_min = _a * _tau_max;

			// �������� ���������� 
			#pragma omp parallel for 
			for (int i = 0; i < _n_cities; ++i)
				for (int j = 0; j < _n_cities; ++j)
					_tau[i][j] = max(_tau_min, _tau[i][j] * (1 - _rho));

			// ��������� ������ �������
			double w = 1.0 / _len;

			#pragma omp parallel for 
			for (int j = 0; j < _n_cities; ++j)
				_tau[_solution[j]][_solution[j + 1]] = min(_tau_max, _tau[_solution[j]][_solution[j + 1]] + w);

			system_clock::time_point end = system_clock::now();
			duration <double> delta = end - start;

			fout << delta.count() << " ";
		}

		fout << "\n";
	}

public:

	/// <summary>
	/// �����������
	/// </summary>
	/// <param name="type"> ��� ������������ ��������� </param>
	/// <param name="params"> ������� ���������� ��������� </param>
	ACOSolver(string type, map<string, any> params) : _type(type)
	{
		// ������������ ���������
		_alpha = any_cast<double>(params["alpha"]);
		_beta = any_cast<double>(params["beta"]);
		_rho = any_cast<double>(params["rho"]);
		_tau0 = any_cast<double>(params["tau0"]);
		
		// ������������ ���������
		_n_ants = any_cast<int>(params["n_ants"]);
		_max_iter = any_cast<int>(params["max_iter"]);
		_n_jobs = any_cast<int>(params["n_jobs"]);

		// � ������ ������� ���������� �������
		if (type == "EAS" || type == "ASRank") _w = any_cast<int>(params["w"]);

		// � ������ ���������� ���������� �������
		if (type == "MMAS")
		{
			_a = any_cast<double>(params["a"]);

			_local_search_type = any_cast<string>(params["local_search_type"]);

			// � ������ ������� ���������� ������ ������������ �������
			if (_local_search_type != "None")
			{
				_local_search_tours = any_cast<string>(params["local_search_tours"]);
				_k = any_cast<int>(params["k"]);
			}
		}
	}

	/// <summary>
	/// ������ ������
	/// </summary>
	/// <param name="g"> ���� </param>
	void solve(Graph& g)
	{
		// ���� �������������:

		_len = INF;
		_n_cities = g.n();
		omp_set_num_threads(_n_jobs);

		_tau = vector<vector<double>>(_n_cities, vector<double>(_n_cities, _tau0));
		_eta_beta = vector<vector<double>>(_n_cities, vector<double>(_n_cities));
		_weights = vector<vector<double>>(_n_cities, vector<double>(_n_cities));

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			for (int j = 0; j < _n_cities; ++j)
				_eta_beta[i][j] = pow(g[i][j], -_beta);

		_vertices = vector<int>(_n_cities);

		#pragma omp parallel for
		for (int i = 0; i < _n_cities; ++i)
			_vertices[i] = i;

		_ants = vector<Ant>(_n_ants);

		_choices = vector<vector<int>>(_n_jobs, vector<int>(_n_cities));
		_visited = vector<vector<int>>(_n_jobs, vector<int>(_n_cities));

		if (_type == "AS" || _type == "EAS") elitist_ant_system(g);
		else if (_type == "ASRank") rank_based_ant_system(g);
		else if (_type == "MMAS") max_min_ant_system(g);
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