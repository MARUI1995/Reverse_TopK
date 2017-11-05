#include  "engine.h"
#include <string>
#include <stdio.h>
#include "dos.h"
#include <iostream>
using namespace std;
#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libeng.lib")
#define doubleSize(count) count * sizeof(double)

Engine* openEng()
{
	Engine* pEng  = NULL;
	if (!(pEng = engOpen(NULL)))
		printf("Open matlab enging fail!");
	return pEng;
}

/*
*Description: ���þ���H, �Խ����ϵ���Ϊ2�� ������Ϊ0
*Parameters:
* H: float*
* column: int
*/
void setH(double* H, int column)
{
	for ( int i = 0; i < column; i++ )
		for( int j = 0; j < column; j++ )
		{
			H[i * column + j] = 0.0;
			H[i * column + i] = 2.0;
		}
}

/*
*Description: ����f��f = (-2)*f
*Parameters:
* f: float*
* dimension: ������ά��
*/
void setf(double* f, int dimension)
{
	for( int i = 0; i < dimension; i++ )
		f[ i ] = ( -2 ) * f[ i ];
}

/*
*Description: �����������ĸ�ά���
*Parameters:
* a: double* ��һ������
* b: doulbe* �ڶ�������
* dimension: int ������ά��
* Return: double* ������������Ľ��
*/
double* vector_substraction(const double a[], const double b[], int count)
{
	if( a == NULL )
		return NULL;
	double* result = new double[count];
	int i = -1;
	while( ++i < count )
		result[i] = a[i] - b[i];
	return result;
}

/*
*Description: ������������ˣ�������һά����ķ�ʽ���д洢
*Parameters:
* a: float* ��һ������
* b: float* �ڶ�������
* count: int �����������˼�������
* dimension: int ������ά��
* Return: float ����������˵Ľ��
*/
void vector_multiply(double* result, const double a[], const double b[], int count, int dimension)
{
	int i = -1;
	int cnt = 0;
	while( ++i < dimension* count )
	{
		double tmp = 0.0;
		int j = -1;
		while( ++j < dimension )
		{
			tmp += a[i] * b[i];
			i++;
		}
		i--;
		result[cnt] = tmp;
		cnt++;
	}
}

/*
*Description: ����һ����ʾ�����һά���飬��������ת�ɰ����к��д�ŵ���ʽ
*Parameters:
* w: ��һ����row�У�column����ɵĶ�ά���������к��еķ�ʽ���д洢��һ��һά����
* row:
* column: 
*/
void reverseVectorToCR(double* rcVector, int row, int column)
{
	double* tmp = new double[row * column];
	for( int i = 0; i < row; i++ )
		for( int j = 0; j < column; j++ )
		{
			int c_pos = i * column + j ;
			int m_pos = j * row + i;
			tmp[ m_pos ] = rcVector[ c_pos ];
		}
		memcpy(rcVector, tmp, doubleSize( row * column ));
		delete [] tmp;
}

/*
*Description: ����һ����ʾ�����һά���飬��������ת�ɰ��������д�ŵ���ʽ
*Parameters:
* arr: ��һ����row�У�column����ɵĶ�ά���������к��еķ�ʽ���д洢��һ��һά����
* row:
* column: 
*/
void reverseVectorToRC( int* arr, int row, int column )
{
	int* tmp = new int[ row * column ];
	for( int i = 0; i < row; i++ )
	{
		for( int j = 0; j < column; j++ )
		{
			int c_pos = i * column + j ;
			int m_pos = j * row + i;
			tmp[ c_pos ] = arr[ m_pos ];
		}
	}
	memcpy((void*)arr, (void*)tmp, row*column*sizeof(float));
	delete[] tmp;
}

/*
*Description:�޲��յ�ʱ���޸Ĳ�ѯ��Pʱ���ã�ʹ�� quadprog(H,f,A,b,lb,ub)
*Parameter:
* q: ԭ���Ĳ�ѯ��
* w: ��ʾ��ʧ��ʸ��
* wp: ��Ӧһ��w�����ڵ�top-k���ĵ�
* count: w �д�ŵ������м���
* dimension: w�д�ŵ������Ǽ�ά��
*Return: �޸ĺ�Ĳ�ѯ�� q
*/
float* q_Quadprog(float* q, float* w, float* wp, int count, int dimension)
{
	Engine* pEng = openEng();
	int num = count * dimension;

	double *_q, *_w, *_wp; 
	_q = new double[dimension];
	_w = new double[count * dimension];
	_wp = new double[count * dimension];

	 int i = 0;
	//��float��ָ��ָ�������ת�浽double��ָ����
	for(; i < count*dimension; i++)
	{
		if( i < dimension )
			_q[i] = q[i];
		_w[i] = (double)w[i];
		_wp[i] = (double)wp[i];
	}

	if ( pEng == NULL )
		return NULL;

	//����H
	mxArray *H = NULL;
	int _H_size = dimension * dimension;
	double* _H = new double[_H_size];

	setH(_H, dimension);
	H = mxCreateDoubleMatrix( dimension, dimension, mxREAL );
	memcpy((void*)mxGetPr(H), (void*)_H, doubleSize(_H_size));
	engPutVariable(pEng, "H", H);

	// ����f
	double* _f = new double[dimension];
	memcpy(_f, _q, doubleSize(dimension));
	setf(_f, dimension);
	mxArray *f = NULL;
	f = mxCreateDoubleMatrix(dimension, 1, mxREAL);
	memcpy((void*)mxGetPr(f), (void*)_f, doubleSize(dimension));
	engPutVariable(pEng, "f", f);

	// ���� A
	mxArray *A = NULL;
	double* _A = new double[ num ];
	memcpy(_A, _w, doubleSize(num));
	reverseVectorToCR(_A, count, dimension );
	A = mxCreateDoubleMatrix( count, dimension, mxREAL);
	memcpy((void*)mxGetPr(A),(void*)_A, doubleSize(num));
	engPutVariable(pEng, "A", A);

	// ����b
	mxArray *b = NULL;
	double* _b = new double[count];
	vector_multiply(_b, _wp, _w, count, dimension);
	b = mxCreateDoubleMatrix( count, 1, mxREAL);
	memcpy((void*)mxGetPr(b), (void*)_b,  doubleSize(count));
	engPutVariable(pEng, "b", b);

	//���� ub
	mxArray* ub = NULL;
	ub = mxCreateDoubleMatrix( dimension, 1, mxREAL );
	memcpy((void*)mxGetPr(ub), (void*)_q, doubleSize(dimension));
	engPutVariable(pEng, "ub", ub);

	//����lb
	mxArray* lb = NULL;
	double* _lb = new double[dimension];
	memset(_lb, 0, doubleSize(dimension));
	lb = mxCreateDoubleMatrix ( dimension, 1, mxREAL );
	memcpy((void*)mxGetPr(lb), (void*)_lb, doubleSize(dimension));
	engPutVariable(pEng, "lb", lb);

	//ִ��
	char out[3000];
	engOutputBuffer(pEng, out, 3000);
	engEvalString(pEng, "modified_q = quadprog (H,f,A,b,[],[],lb,ub);");
	mxArray* result;
	result = engGetVariable(pEng, "modified_q");

	double* _modified_q = mxGetPr(result);
	float* modified_q = new float[dimension];
	i = -1;
	while( ++i < dimension )
	{
		modified_q[i] = (float)_modified_q[i];
		if( modified_q[i] < 0 )
		{
			cout<<modified_q[i]<<endl;
			modified_q[i] = 0;
		}
	}
//	engClose(pEng);
	delete[] _q;
	delete[] _w;
	delete[] _wp;
	delete[] _H;
	delete[] _A;
	delete[] _f;
	delete[] _b;
	delete[] _lb;
	return modified_q;
}

///*
//* Description: �ڱ߽���ȡ��ʱ����һ����������
//* Parameters:
//*  factors: �������̵�ϵ������һά�����ʾ
//*  maxValue: δ֪����ȡֵ���Ͻ�
//*  dimension: ʸ����ά�ȣ���δ֪���ĸ���
//*  count: Ҫȡ���ĸ���
//*/
//int* calIndefEqual(int* factors, int maxValue, int dimension, int count)
//{
//	char symbols[5] = {'x', 'y', 'z', 'm', 'n'};  // ���ű���
//	cout<<symbols[0];
//	char commStr[100];//����matlab�������ַ���
//
//	// �� matlab engine
//	Engine* pEng = openEng();
//
//	memset(commStr,'\0',sizeof(commStr));
//   //------------------------------------
//   //��matlab�ж������
//   //------------------------------------
//	strcpy(commStr,"[");
//	for( int i =0 ; i < dimension; i++ )
//			sprintf(commStr,"%s%c,",commStr,symbols[i]);
//    commStr[strlen(commStr)-1]='\0';
//	sprintf(commStr, "%s]=ndgrid(0:%d);",commStr, maxValue);
//	engEvalString(pEng,commStr);
//
//	memset(commStr,'\0',sizeof(commStr));
//	//-----------------------------------
//	// ִ�� find ��������������Ϊ�� ax + by + ...  = 0
//   //------------------------------------
//	strcpy(commStr,"result = find(");
//	for( int i =0 ; i < dimension; i++ )
//		sprintf(commStr,"%s%d*%c+",commStr,factors[i],symbols[i]);
//	commStr[strlen(commStr)-1]='\0';
//	strcat(commStr, "==0);");
//	engEvalString(pEng,commStr);
//
//	memset(commStr,'\0',sizeof(commStr));
//	//------------------------------------
//	//  ��ִ�н��д����� results
//   //-------------------------------------
//	strcpy(commStr,"results=[");
//	for( int i =0 ; i < dimension; i++ )
//		sprintf(commStr,"%s%c(result),",commStr,symbols[i]);
//	commStr[strlen(commStr)-1]='\0';
//	strcat(commStr, "];");
//	engEvalString(pEng,commStr);
//
//	memset(commStr,'\0',sizeof(commStr));
//	//------------------------------------
//	//���ҽ�������˳��
//   //-------------------------------------
//	sprintf(commStr,"n=randperm(size(results))");
//	engEvalString(pEng,commStr);
//
//	memset(commStr,'\0',sizeof(commStr));
//	//------------------------------------
//	// �Ӵ���˳���ľ���ȡ��ǰ������
//   //-------------------------------------
//	sprintf(commStr," results=results(n(1:%d),:)",count);
//	engEvalString(pEng,commStr);
//
//
//	int* results = new int[count * dimension];
//	results = (int*)(mxGetPr(engGetVariable(pEng, "results")));
//	reverseVectorToRC(results, count, dimension);
//		engClose(pEng);
//	return results;
//}