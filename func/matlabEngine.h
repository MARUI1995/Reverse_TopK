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
*Description: 设置矩阵H, 对角线上的数为2， 其它的为0
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
*Description: 设置f，f = (-2)*f
*Parameters:
* f: float*
* dimension: 向量的维数
*/
void setf(double* f, int dimension)
{
	for( int i = 0; i < dimension; i++ )
		f[ i ] = ( -2 ) * f[ i ];
}

/*
*Description: 将两个向量的各维相减
*Parameters:
* a: double* 第一个向量
* b: doulbe* 第二个向量
* dimension: int 向量的维度
* Return: double* 两个向量相减的结果
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
*Description: 将两个矩阵相乘，矩阵以一维数组的方式进行存储
*Parameters:
* a: float* 第一个向量
* b: float* 第二个向量
* count: int 数组里面存放了几个向量
* dimension: int 向量的维度
* Return: float 两个向量相乘的结果
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
*Description: 输入一个表示矩阵的一维数组，将该数组转成按先列后行存放的形式
*Parameters:
* w: 将一个用row行，column列组成的多维数组以先列后行的方式进行存储的一个一维数组
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
*Description: 输入一个表示矩阵的一维数组，将该数组转成按先行再列存放的形式
*Parameters:
* arr: 将一个用row行，column列组成的多维数组以先行后列的方式进行存储的一个一维数组
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
*Description:无参照点时，修改查询点P时调用，使用 quadprog(H,f,A,b,lb,ub)
*Parameter:
* q: 原来的查询点
* w: 表示丢失的矢量
* wp: 对应一个w，处于第top-k个的点
* count: w 中存放的向量有几个
* dimension: w中存放的向量是几维的
*Return: 修改后的查询点 q
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
	//将float型指针指向的数据转存到double型指针中
	for(; i < count*dimension; i++)
	{
		if( i < dimension )
			_q[i] = q[i];
		_w[i] = (double)w[i];
		_wp[i] = (double)wp[i];
	}

	if ( pEng == NULL )
		return NULL;

	//设置H
	mxArray *H = NULL;
	int _H_size = dimension * dimension;
	double* _H = new double[_H_size];

	setH(_H, dimension);
	H = mxCreateDoubleMatrix( dimension, dimension, mxREAL );
	memcpy((void*)mxGetPr(H), (void*)_H, doubleSize(_H_size));
	engPutVariable(pEng, "H", H);

	// 设置f
	double* _f = new double[dimension];
	memcpy(_f, _q, doubleSize(dimension));
	setf(_f, dimension);
	mxArray *f = NULL;
	f = mxCreateDoubleMatrix(dimension, 1, mxREAL);
	memcpy((void*)mxGetPr(f), (void*)_f, doubleSize(dimension));
	engPutVariable(pEng, "f", f);

	// 设置 A
	mxArray *A = NULL;
	double* _A = new double[ num ];
	memcpy(_A, _w, doubleSize(num));
	reverseVectorToCR(_A, count, dimension );
	A = mxCreateDoubleMatrix( count, dimension, mxREAL);
	memcpy((void*)mxGetPr(A),(void*)_A, doubleSize(num));
	engPutVariable(pEng, "A", A);

	// 设置b
	mxArray *b = NULL;
	double* _b = new double[count];
	vector_multiply(_b, _wp, _w, count, dimension);
	b = mxCreateDoubleMatrix( count, 1, mxREAL);
	memcpy((void*)mxGetPr(b), (void*)_b,  doubleSize(count));
	engPutVariable(pEng, "b", b);

	//设置 ub
	mxArray* ub = NULL;
	ub = mxCreateDoubleMatrix( dimension, 1, mxREAL );
	memcpy((void*)mxGetPr(ub), (void*)_q, doubleSize(dimension));
	engPutVariable(pEng, "ub", ub);

	//设置lb
	mxArray* lb = NULL;
	double* _lb = new double[dimension];
	memset(_lb, 0, doubleSize(dimension));
	lb = mxCreateDoubleMatrix ( dimension, 1, mxREAL );
	memcpy((void*)mxGetPr(lb), (void*)_lb, doubleSize(dimension));
	engPutVariable(pEng, "lb", lb);

	//执行
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
//* Description: 在边界上取点时，解一个不定方程
//* Parameters:
//*  factors: 不定方程的系数，以一维数组表示
//*  maxValue: 未知数的取值的上界
//*  dimension: 矢量的维度，即未知数的个数
//*  count: 要取样的个数
//*/
//int* calIndefEqual(int* factors, int maxValue, int dimension, int count)
//{
//	char symbols[5] = {'x', 'y', 'z', 'm', 'n'};  // 符号变量
//	cout<<symbols[0];
//	char commStr[100];//传入matlab的命令字符串
//
//	// 打开 matlab engine
//	Engine* pEng = openEng();
//
//	memset(commStr,'\0',sizeof(commStr));
//   //------------------------------------
//   //在matlab中定义变量
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
//	// 执行 find 函数，不定方程为： ax + by + ...  = 0
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
//	//  把执行结果写入矩阵 results
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
//	//打乱结果矩阵的顺序
//   //-------------------------------------
//	sprintf(commStr,"n=randperm(size(results))");
//	engEvalString(pEng,commStr);
//
//	memset(commStr,'\0',sizeof(commStr));
//	//------------------------------------
//	// 从打乱顺序后的矩阵取出前几数组
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