#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>;
#include <windows.h>;

using namespace std;

double eps = 0.0000001;
const int N = 4;

void swap_rows(double A[N][N], int n, int x, int y)
{
	for (int i = 0; i < n; i++) {
		swap(A[x][i], A[y][i]);
	}
}

void swap_columns(double A[N][N], int n, int x, int y)
{
	for (int i = 0; i < n; i++) {
		swap(A[i][x], A[i][y]);
	}
}


bool LU_full(double A[N][N], double L[N][N],
	double U[N][N], int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			U[i][j] = A[i][j];
		}
	}
	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++) {
			if (fabs(U[i][i]) <= eps) {
				int ansx = i;
				int ansy = i;
				for (int j = i; j < n; j++) {
					for (int k = i; k < n; k++) {
						if (fabs(U[j][k]) > fabs(U[ansx][ansy])) {
							ansx = j;
							ansy = k;
							break;
						}
					}
				}
				if (i != ansx) {
					swap_rows(U, n, i, ansx);
					swap_rows(L, n, i, ansx);
				}
				else {
					swap_columns(U, n, i, ansy);
					swap_columns(L, n, i, ansy);
				}
			}
			if (fabs(U[i][i]) <= eps) {
				return false;
			}
			for (int j = i; j < n; j++) {
				L[j][i] = U[j][i] / U[i][i];
			}
		}
		for (int i = k; i < n; i++) {
			for (int j = k - 1; j < n; j++) {
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
			}
		}
	}
	return true;
}

void proisv(double A[N][N], double B[N][N],
	double R[N][N], int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
}

void show(double A[N][N], int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%8.5lf ", A[i][j]);
		}
		cout << endl;
	}
}

double det(double a[N][N], int n)
{
	double res = 1;
	for (int i = 0; i < n; i++) {
		res *= a[i][i];
	}
	return res;
}

int main()
{
	// ��������� �� ������� ������� 
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	ifstream fin("input.txt");

	double A[N][N] = { 0 };
	double L[N][N] = { 0 };
	double R[N][N] = { 0 };
	double U[N][N] = { 0 };
	double temp[N][N] = { 0 };
	int n = N;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			fin >> A[i][j];
			temp[i][j] = A[i][j];
		}
	}
	cout << "�������� ����� ����" << endl;
	cout << "1 - LU ����������" << endl;
	cout << "2 - �����" << endl;
	int key;
	cin >> key;
	bool res;
	if (key == 1) {
		res = LU_full(A, L, U, n);
	}
	if (key == 2) {
		exit(1);
	}
	if (!res) {
		cout << "No solve" << endl;
		system("pause");
		return 0;
	}

	cout << "����������� �������" << endl << det(U, n) << endl;
	cout << "���� �������" << endl;
	show(A, n);
	cout << "U �������" << endl;
	show(U, n);
	cout << "L �������" << endl;
	show(L, n);
	proisv(L, U, R, n);
	cout << "L*U �������" << endl;
	show(R, n);
	system("pause");
	return 0;
	
}
