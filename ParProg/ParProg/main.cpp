#include <iostream>
#include "omp.h"


using namespace std;

double CannonAlgorithm(double * AM, double * BM, double * CM, int N, int B)
{
	double start = 0, finish = 0;
	
	start = omp_get_wtime();
	for (int block = 0; block < B*B; block++)
	{
		int ib = block / B;
		int jb = block - ib * B;

		int jj = (ib + jb) % B;

		for (int k = 0; k < B; k++)
		{
			for (int mi = ib * (N / B); mi < (ib + 1)*(N / B); mi++)
				for (int mj = jb * (N / B); mj < (jb + 1)*(N / B); mj++)
				{
					double SUM = 0.0;
					for (int mk = jj * (N / B); mk < (jj + 1)*(N / B); mk++)
						SUM += AM[mi*N + mk] * BM[mk*N + mj];
					CM[mi*N + mj] += SUM;
				}

			jj = (jj + 1) % B;
		}
	}
	finish = omp_get_wtime();

	return (finish - start);
}


void main()
{
	int B = 3;			//grid size B*B
	int N = B * 2;		//matrix size: N*N, block size 2*2

	double *AM = new double[N*N];
	double *BM = new double[N*N];
	double *CM = new double[N*N];

	for (int i = 0; i < N*N; i++)
	{
		AM[i] = rand() % 5 + 1;
		BM[i] = rand() % 5 + 1;
		CM[i] = 0.0;
	}

	double TIME = 0.0;
	TIME = CannonAlgorithm(AM, BM, CM, N, B);
	cout << "TIME: " << TIME << endl;

	bool visual = true;
	if (visual) 
	{
		cout << "A:" << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout << AM[i*N + j] << "\t";
			}
			cout << endl;
		}
		cout << endl;
		cout << "B:" << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout << BM[i*N + j] << "\t";
			}
			cout << endl;
		}
		cout << endl;
		cout << "C:" << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout << CM[i*N + j] << "\t";
			}
			cout << endl;
		}
	}
	
	cin.get();
}