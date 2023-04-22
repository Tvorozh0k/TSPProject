#pragma once
#include "TSPSolver.h"

/// <summary>
/// ������ �������� ������� ������
/// ������������: ����� ���������� ������
/// </summary>
class NNSolver : public TSPSolver
{
public:

	/// <summary>
	/// ������ ������
	/// </summary>
	void solve(Graph& g)
	{
		_size = g.n();

		// visited[i] = 0 <=> ������� i �� ��������
		// visited[i] = 1 <=> ������� i ��������
		vector <int> visited(_size);

		// ���������� ������ ������� ��������
		uniform_int_distribution <int> interval_int(0, _size - 1);
		int cur = interval_int(gen);

		_solution.resize(_size + 1);

		// ����: cur -> ... 
		_solution[0] = cur; _len = 0;
		visited[cur] = 1;

		// ����������� ���� �� �����
		for (int i = 1; i < _size; ++i)
		{
			// ���� ������� next
			int next = cur, dist = INF;

			// ���������� �� ���� �� ���������� �������� j 
			// � ���� �������� ������� � ������� cur
			for (int j = 0; j < _size; ++j)
				if (!visited[j] && g[cur][j] < dist)
				{
					next = j;
					dist = g[cur][j];
				}

			// ����: ... -> cur -> next -> ...
			_solution[i] = next; _len += dist;
			visited[next] = 1;

			cur = next;
		}

		// ����: s -> ... -> t -> s
		int s = _solution[0], t = _solution[_size - 1];
		_solution[_size] = s; _len += g[t][s];
	}

	/// <summary>
	/// ������� ������� �� �����
	/// </summary>
	void print()
	{
		cout << "NNSolver (��������� �������): ";

		for (int i = 0; i < _size; ++i)
			cout << _solution[i] << " - ";

		cout << _solution.back() << ", ����� ���� = " << _len << "\n";
	}
};