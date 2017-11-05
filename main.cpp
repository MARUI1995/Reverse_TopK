#include ".//rtree/rtree.h"
#include ".//blockfile/cache.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <io.h>
#include<vector>
#include<windows.h>
#include<unordered_map>
// #define MAX_SKYLINE_RESULT	102400
using namespace std;
char setFileName[100];
string rslt_fname("");
const int MAX_FILE_NAME = 512;
RTree* buildtree(char *_fname, int _dimension, Cache* c);
template <class T>
void writeSet(T*, int, int, char*);
template <class T>
void readSet(T* set, char* filePath, int dimension);
int readSet(float* queryPoints, float* weightSet, char* filePath, int dimension, int missCnt, int * missIdSet);
ofstream outfile;
class qPoint
{
public:
	float *m_qry_pt;
	float *m_reWeight;
	float *m_notReWeight;
	//float *m_minParameter;
	float *m_missID;
	float *m_tableID;
	float *m_atrriContent;

	//float m_spjTime;
	int m_dimension;
	int m_reWeightcnt;
	int m_notReWeightcnt;
	int m_missCnt;
	int m_tableCnt; //连接表的个数作为实验参数给定
	int m_attriCnt; //属性个数作为参数给定
	
	qPoint(float *_qry_pt, int dimension, int missCnt, int attriCnt, int tableCnt)
		:m_dimension(dimension),
		m_reWeightcnt(0),
		m_notReWeightcnt(0),
		m_missCnt(missCnt),
		m_tableCnt(tableCnt),
		m_attriCnt(attriCnt)
		//m_spjTime(0)
	{
		m_qry_pt = new float[dimension];
		memcpy(m_qry_pt, _qry_pt, sizeof(float) * dimension);
		m_reWeight = new float[dimension * 100];
		m_notReWeight = new float[dimension * 100];
		m_missID = new float[missCnt];
		m_tableID = new float[tableCnt];
		m_atrriContent = new float[m_attriCnt * 2];
		//m_minParameter = new float[100];

	}
};
class attribute
{
public:
	int m_dimension;
	float * m_atrribute;
	int m_id;
	attribute(int dimension,float *attribute,int id)
		:m_dimension(dimension),
		m_id(id)
	{
		m_atrribute = new float[dimension];
		memcpy(m_atrribute, attribute, sizeof(float) * dimension);
	}
};
void ifReverse(vector<qPoint> &qvector, int qryCnt, float* queryPoints,float* weightSet, int missWeightCnt, int dimension, Cache *c, RTree *rtree, Heap *hp, int k, int *missIdSet)
{
	//------------------------------begin-------------------------------------------------------------------
	// 执行查看查询点是否是反Top-k结果
	//---------------------------------------------------------------------	
	float* missWeightSet = new float[missWeightCnt * dimension];//对当前数据集的单个查询点的一组丢失权重向量，可能是一个，可能是多个，根据数据集情况，每一行是一个权重
	float *_qry_pt = new float[dimension];//对当前数据集的一个查询点
	qPoint *qpoint;//一个查询点
	int * missId;//一个查询点的缺失向量所有ID
	//float *parameter;
	clock_t begin_tick, end_tick;
	for (int cnt = 0; cnt < qryCnt; cnt++) //按照查询点个数进行循环
	{
		//parameter = new float[100];
		//outfile << "第\t" << cnt << "\t个查询点:\n";
		missId = new int[missWeightCnt];
		memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//赋值给查询点  _qry_pt
		memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //赋值给缺失权重
		//memcpy(parameter, &missIdSet[cnt*1*missWeightCnt], sizeof(float) * 1 * missWeightCnt);//赋值给权重属性
		memcpy(missId, &missIdSet[cnt * 1 * missWeightCnt], sizeof(int) * 1 * missWeightCnt);//赋值给权重ID
		c->page_faults = 0;//cache c
		qpoint = new qPoint(_qry_pt, dimension, missWeightCnt);
		qpoint->m_missID = missId;
		//------------------------------begin-------------------------------------------------------------------
		// 赋值给查询点缺失权重的属性
		//---------------------------------------------------------------------	
		/*begin_tick = clock();
		float *currentPa = new float[1];
		float *minPa = new float[1];
		memcpy(minPa, parameter, sizeof(float) * dimension);//初始化
		for (int missweightcnt = 0; missweightcnt < missWeightCnt; missweightcnt++)
		{
			memcpy(currentPa, &parameter[missweightcnt *1], sizeof(float) * dimension);//赋值给查询点  _qry_pt
			for (int j = 0; j < 1; j++)
			{
				if (currentPa[j] < minPa[j])
					minPa[j] = currentPa[j];
			}
		}
		end_tick = clock();*/
		//memcpy(qpoint->m_minParameter, minPa, sizeof(float) * dimension);//初始化
		//qpoint->m_spjTime = (end_tick - begin_tick) / (float)1000;


		cout << "赋值给查询点时：";
		for (int i = 0; i<dimension; i++)
		{
			cout << _qry_pt[i] << " ";

		}
		cout << endl;
		cout << "建立相应qpoint：";
		for (int i = 0; i<dimension; i++)
		{
			cout << qpoint->m_qry_pt[i] << " ";

		}

		
		cout << endl;

		/*计算当前查询点的主导点和不可比点*/
		int ranking, rank;
		int bedomnum, incomnum;  //incomnum是incomparable点的个数，bedomnum是dominate查询点q的点的个数，
		float *incomparable = new float[3000000];
		rtree->incomparable_points_reuse(hp, incomparable, incomnum, _qry_pt, bedomnum);   //求incomparable点的集合，放在数组incomparable中，并求出incomnum和 bedomnum的值
		//outfile << "主导点有" << bedomnum - 1 << "个" << endl;//因为主导点包含了查询点自己
		//outfile << "不可比点有" << incomnum << "个" << endl;
		for (int j = 0; j < incomnum; j++)
		{


			for (int i = 0; i < dimension; i++)
			{
				//outfile << incomparable[j*dimension + i] << "　";

			}
			//outfile << endl;
		}

		float * currentWeight = new float[dimension];

		for (int misweightCnt = 0; misweightCnt < missWeightCnt; misweightCnt++)
		{
			memcpy(currentWeight, &missWeightSet[misweightCnt * dimension], sizeof(float) * dimension);//赋值给查询点  _qry_pt
			rtree->bbs_ranking_querypoint(incomparable, incomnum, _qry_pt, currentWeight, ranking);//此处传入currentWeight必须是单个权重
			rank = ranking + bedomnum;
			//outfile << "当前rank=" << rank << "　" << "当前Topk值=" << k << endl;

			if (rank <= k)
			{
				//outfile << "==========是反Top-k结果，不用改正==========" << "rank=" << rank << endl;
				memcpy(qpoint->m_reWeight + qpoint->m_reWeightcnt*dimension, currentWeight, sizeof(float) * dimension);
				qpoint->m_reWeightcnt++;

			}
			else
			{
				//outfile << "==========不是反Top-k结果，需要改正==========" << "rank=" << rank << endl;
				memcpy(qpoint->m_notReWeight + qpoint->m_notReWeightcnt*dimension, currentWeight, sizeof(float) * dimension);
				qpoint->m_notReWeightcnt++;
			}
		}
		qvector.push_back(*qpoint);

		cout << "加入一个qpoint" << endl;
		cout << "此查询点有" << qpoint->m_missCnt << "个缺失结果" << "    " << "缺失结果的ID分别为";
		for (int i = 0; i < qpoint->m_missCnt; i++)
		{
			cout << qpoint->m_missID[i] << " " << endl;
		}
		cout << "当前qvector大小为" << qvector.size() << endl;
		cout << endl;
	}
}
void modifyKW(vector<qPoint> &qvector, char *fname, int qryCnt, int missWeightCnt, float* queryPoints, float* weightSet, Cache *c, RTree *rtree, int dimension, Heap *hp, long samplesize, float  quality_of_answer, float probability_guarantee,int k, int &modifyk)
{
	float* missWeightSet = new float[missWeightCnt * dimension];//对当前数据集的单个查询点的一组丢失权重向量，可能是一个，可能是多个，根据数据集情况，每一行是一个权重
	float *_qry_pt = new float[dimension];//对当前数据集的一个查询点
	float *modifyweight = new float[dimension * missWeightCnt];//修改后的权重向量，每一列是一个权重
	outfile.open(setFileName, ios::out | ios::app);
	cout << "进入修改KW" << endl;
	outfile << "========== 修改 k w 的近似算法 ==========\n";
	outfile << "共有\t" << qryCnt << "\t个查询点\n";
	outfile << "数据集：" << fname << endl;
	outfile << endl;
	clock_t begin_tick, end_tick;
	float pruning = 0.0, prunings = 0.0;
	float penalty = 0.0, penaltys = 0.0;
	float ticksCnt = 0, page_faults = 0;
	int totalMissCnt = missWeightCnt;
	//float * minPa;
	//float spjTime = 0;
	//if(qvector.size()==qryCnt)
	//cout<<"qvector大小与查询点数目相同"<<endl;

	for (int cnt = 0; cnt < qryCnt; cnt++) //按照查询点个数进行循环
	{
		outfile << "第\t" << cnt << "\t个查询点:\n";
		//memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//赋值给查询点  _qry_pt
		//memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //赋值给缺失权重，此处被赋值的权重必须是非反Top-k结果。
		_qry_pt = qvector[cnt].m_qry_pt;
		missWeightSet = qvector[cnt].m_notReWeight;
		missWeightCnt = qvector[cnt].m_notReWeightcnt;
		//minPa = qvector[cnt].m_minParameter;
		//spjTime = qvector[cnt].m_spjTime;
		c->page_faults = 0;//cache c
		begin_tick = clock();
		rtree->modifyWandK_approximate_reuse(hp, modifyweight, modifyk, penalty, _qry_pt, missWeightSet, missWeightCnt, k, quality_of_answer, probability_guarantee, pruning);//修改kw
		end_tick = clock();

		//将结果写入文件
		outfile << "---------- 参数设置----------\n";
		outfile << "查询点q：( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << _qry_pt[i] << " , ";
		outfile << _qry_pt[dimension - 1] << " )" << endl;   //输出当前查询点的各个维度的值

		outfile << "反top-k查询的k值： " << k << endl;  //原k值
		outfile << "维度d = " << dimension << endl;
		outfile << "丢失的权重向量集合大小=" << totalMissCnt << endl;
		outfile << "丢失的非反Top-k权重集合大小|W| = " << missWeightCnt << ",需要修正" << endl; //缺失权重个数
		outfile << "需要修正的向量集合:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			//outfile << "在why-not向量下的rank = " << rank << " : ( ";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << missWeightSet[i * dimension + j] << " , ";
			outfile << missWeightSet[j] << " )" << endl;
		}

		outfile << "Pr = " << probability_guarantee << endl;
		outfile << "T = " << quality_of_answer << endl;

		outfile << "\n----------执行结果----------" << endl;
		outfile << "修正后向量集合:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			outfile << i << ". (";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << modifyweight[i*dimension + j] << " , ";
			outfile << modifyweight[j] << " )" << endl;
		}
		outfile << "k = " << modifyk << endl;
		//outfile << "IO代价 = " << c->page_faults / (float)100 << endl;
		//outfile << "CPU代价 =  " << (end_tick - begin_tick) / (float)1000 << endl;

		//删除掉SPJ原始修改
		/*outfile << "修正后spj部分属性：" << endl;
		for (int i = 0; i < parameterCnt; i++)
		{
			outfile << minPa[i] << "	";
		}
		*/
		outfile << endl;
		outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
		outfile << "penalty = " << penalty << endl;
		//outfile << "剪枝效果 = " << pruning << endl;
		//outfile << "剪枝个数 = " << (int)(samplesize * pruning) << endl;
		ticksCnt += (end_tick - begin_tick) ;
		page_faults += (c->page_faults);
		prunings += pruning;
		penaltys += penalty;



		outfile << "--------------------\n";
	}
	ticksCnt = ticksCnt / (float)1000 / (float)qryCnt;
	page_faults = page_faults / (float)100 / (float)qryCnt;
	penaltys = penaltys / (float)qryCnt;
	prunings = prunings / (float)qryCnt;
	outfile << "==========KW平均结果==========\n";
	//outfile << "IO代价 = " << page_faults << endl;
	//outfile << "CPU代价 =  " << ticksCnt << endl;
	outfile << "TotalTime = " << page_faults + ticksCnt << endl;
	outfile << "penalty = " << penaltys << endl;
	//outfile << "剪枝效果 = " << prunings << endl;
	outfile << "取样个数 = " << samplesize << endl;
	//outfile << "剪枝个数 = " << (int)(samplesize * prunings) << endl;
	outfile << "========================================\n\n";
	outfile.flush();
}

void modifyQ(queue< float* > &qmodifiedQueue, vector<qPoint> &qvector, int dimension, Cache *c, RTree *rtree, char *fname, int qryCnt, int missWeightCnt, float* queryPoints, float* weightSet, Heap *hp, int k)
{
	// outfile.open(setFileName, ios::out | ios::app);
	cout << "进入修改q" << endl;
	float* missWeightSet = new float[missWeightCnt * dimension];//对当前数据集的单个查询点的一组丢失权重向量，可能是一个，可能是多个，根据数据集情况，每一行是一个权重
	float *_qry_pt = new float[dimension];//对当前数据集的一个查询点
	float* qmodified = new float[dimension];//修改后的q放在qmodified
	float qabsolute;
	float ticksCnt = 0.0, page_faults = 0.0;
	float penalty = 0.0, penaltys = 0.0;
	clock_t begin_tick, end_tick;
	//float * minPa;
	//float spjTime = 0;
	//memset(setFileName, '\0', 100);
	if (outfile.is_open())
		cout << "open";
	outfile << "========== 修改 q 的精确算法 ==========\n---------- 参数设置----------\n";
	outfile << "共有\t" << qryCnt << "\t个查询点\n";
	outfile << "数据集：" << fname << endl;

	for (int cnt = 0; cnt < qryCnt; cnt++)
	{
		outfile << "第\t" << cnt << "\t个查询点:\n";
		//memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//赋值给查询点
		//memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt);//赋值给缺失向量
		_qry_pt = qvector[cnt].m_qry_pt;
		missWeightSet = qvector[cnt].m_notReWeight;
		missWeightCnt = qvector[cnt].m_notReWeightcnt;
		//minPa = qvector[cnt].m_minParameter;
		//spjTime = qvector[cnt].m_spjTime;
		c->page_faults = 0;

		if (missWeightCnt>0)
		{
			begin_tick = clock();
			rtree->modifyQ_accurate(hp, _qry_pt, missWeightSet, missWeightCnt, k, qmodified);//精确修改Q
			end_tick = clock();
		}
		else
		{
			begin_tick = clock();
			memcpy(qmodified, _qry_pt, sizeof(float) * dimension);
			end_tick = clock();
		}

		qmodifiedQueue.push(qmodified);

		outfile << "查询点q：( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << _qry_pt[i] << " , ";
		outfile << _qry_pt[dimension - 1] << " )" << endl;

		outfile << "反top-k查询的k值： " << k << endl;
		outfile << "维度d = " << dimension << endl;
		outfile << "丢失的非反Top-k向量集合大小|W| = " << missWeightCnt << endl;
		outfile << "丢失的向量集合:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			//outfile << "在why-not向量下的rank = " << rank << " : ( ";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << missWeightSet[i * dimension + j] << " , ";
			outfile << missWeightSet[j] << " )" << endl;
		}

		outfile << "\n----------执行结果----------" << endl;
		outfile << "查询点修改后: (";

		qabsolute = penalty = 0.0;
		for (int i = 0; i < dimension - 1; i++)
		{
			outfile << qmodified[i] << " , ";
			penalty += (qmodified[i] - _qry_pt[i]) * (qmodified[i] - _qry_pt[i]);
			qabsolute += _qry_pt[i] * _qry_pt[i];
		}
		outfile << qmodified[dimension - 1] << " )" << endl;
		penalty += (qmodified[dimension - 1] - _qry_pt[dimension - 1]) * (qmodified[dimension - 1] - _qry_pt[dimension - 1]);
		qabsolute += _qry_pt[dimension - 1] * _qry_pt[dimension - 1];
		penalty = penalty / qabsolute;
		
		//删除掉原始SPJ修正
		/*outfile << "修正后spj部分属性：" << endl;
		for (int i = 0; i < parameterCnt; i++)
		{
			outfile << minPa[i] << "	";
		}
		outfile << endl;
		*/

		//outfile << "IO代价 = " << c->page_faults / (float)100 << endl;
		//outfile << "CPU代价 =  " << (end_tick - begin_tick) / (float)1000 << endl;
		outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
		outfile << "penalty = " << penalty << endl;

		ticksCnt += (end_tick - begin_tick) ;
		page_faults += (c->page_faults);
		penaltys += penalty;
		outfile << "--------------------\n";


	}
	ticksCnt = ticksCnt / (float)1000 / (float)qryCnt;
	page_faults = page_faults / (float)100 / (float)qryCnt;
	penaltys = penaltys / (float)qryCnt;
	outfile << "==========AccurateQ平均结果==========\n";
	//outfile<<"IO代价 = "<< page_faults<<endl;
	//outfile<<"CPU代价 =  "<< ticksCnt<<endl;		
	outfile << "TotalTime = " << page_faults + ticksCnt << endl;
	outfile << "penalty = " << penaltys << endl;
	outfile << "========================================\n\n";
}

void modifyKWQ(queue< float* > &qmodifiedQueue, vector<qPoint> &qvector,  int dimension, Cache *c, RTree *rtree, char *fname, int qryCnt, int missWeightCnt, float* queryPoints, float* weightSet, Heap *hp, int k, long samplesize, float *modifyweight, int modifyk, float  quality_of_answer, float probability_guarantee)
{
	float* missWeightSet = new float[missWeightCnt * dimension];//对当前数据集的单个查询点的一组丢失权重向量，可能是一个，可能是多个，根据数据集情况，每一行是一个权重
	float *_qry_pt = new float[dimension];//对当前数据集的一个查询点
	float* qmodified_kwq = new float[dimension];
	float* qmodified = new float[dimension];
	float ticksCnt = 0.0, page_faults = 0.0;
	float pruning = 0.0, prunings = 0.0;
	float penalty = 0.0, penaltys = 0.0;
	clock_t begin_tick, end_tick;
	//float * minPa;
	//float spjTime = 0;
	int totalMissCnt = missWeightCnt;
	modifyk = k;
	//memset(setFileName, '\0', 100);
	cout << "进入修改KWQ" << endl;
	if (outfile.is_open())
		cout << "open";
	outfile << "========== 修改 q k w 的算法 ==========\n---------- 参数设置----------\n";
	outfile << "共有\t" << qryCnt << "\t个查询点\n";
	outfile << "数据集：" << fname << endl;
	int hpUsed = 0;

	for (int cnt = 0; cnt < qryCnt; cnt++)
	{
		outfile << "第\t" << cnt << "\t个查询点:\n";
		//memcpy( _qry_pt,  &queryPoints[cnt * dimension], sizeof( float ) * dimension );
		//memcpy( missWeightSet,  &weightSet[cnt * dimension * missWeightCnt], sizeof( float ) * dimension * missWeightCnt);
		_qry_pt = qvector[cnt].m_qry_pt;
		missWeightSet = qvector[cnt].m_notReWeight;
		missWeightCnt = qvector[cnt].m_notReWeightcnt;
		//minPa = qvector[cnt].m_minParameter;
		//spjTime = qvector[cnt].m_spjTime;

		pruning = 0.0;
		qmodified = qmodifiedQueue.front();
		qmodifiedQueue.pop();
		c->page_faults = 0;
		hp->used = 0;
		if (missWeightCnt>0)
		{
			begin_tick = clock();
			rtree->modify_k_W_Q_reuse(hp, _qry_pt, qmodified, missWeightSet, missWeightCnt, k, qmodified_kwq, modifyweight, modifyk, quality_of_answer, probability_guarantee, penalty, pruning);
			end_tick = clock();
		}
		else
		{
			begin_tick = clock();
			memcpy(qmodified_kwq, qmodified, sizeof(float) * dimension);
			end_tick = clock();
		}

		//将结果写入文件
		outfile << "查询点q：( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << _qry_pt[i] << " , ";
		outfile << _qry_pt[dimension - 1] << " )" << endl;
		hpUsed += hp->used;
		outfile << "修改Q方案得到的查询点: ( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << qmodified[i] << " , ";
		outfile << qmodified[dimension - 1] << " )" << endl;

		outfile << "反top-k查询的k值： " << k << endl;
		outfile << "维度d = " << dimension << endl;
		outfile << "丢失的权重向量集合大小=" << totalMissCnt << endl;
		outfile << "丢失的非反Top-k权重集合大小|W| = " << missWeightCnt << ",需要修正" << endl; //缺失权重个数
		outfile << "需要修正的向量集合:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			//outfile<<"在why-not向量下的rank = "<<*(qry_ranks_p - i )<<" : ( ";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << missWeightSet[i*dimension + j] << " , ";
			outfile << missWeightSet[j] << " )" << endl;
		}

		outfile << "T = " << quality_of_answer << endl;
		outfile << "Pr = " << probability_guarantee << endl;

		outfile << "\n----------执行结果----------" << endl;
		outfile << "修正后向量集合:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			outfile << i << ". (";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << modifyweight[i*dimension + j] << " , ";
			outfile << modifyweight[j] << " )" << endl;
		}
		outfile << "k = " << modifyk << endl;
		outfile << "修改查询点q：( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << qmodified_kwq[i] << " , ";
		outfile << qmodified_kwq[dimension - 1] << " )" << endl;
		
		//删除原始SPJ修改
		/*outfile << "修正后spj部分属性：" << endl;
		for (int i = 0; i < parameterCnt; i++)
		{
			outfile << minPa[i] << "	";
		}
		outfile << endl;
		*/

		//outfile<<"IO代价 = "<< c->page_faults / ( float )100<<endl;
		//outfile<<"CPU代价 =  "<< ( end_tick - begin_tick ) / ( float )1000<<endl;		
		outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
		outfile << "penalty = " << penalty << endl;
		//outfile<<"剪枝效果 = "<<pruning<<endl;
		//outfile<<"hp resued = "<<hp->used<<endl;
		//outfile<<"剪枝个数 = "<<(int)(samplesize * pruning)<<endl;
		ticksCnt += (end_tick - begin_tick) ;
		page_faults += (c->page_faults);
		prunings += pruning;
		penaltys += penalty;
		outfile << "--------------------\n";
	}
	ticksCnt = ticksCnt / (float)1000 / (float)qryCnt;
	page_faults = page_faults / (float)100 / (float)qryCnt;
	penaltys = penaltys / (float)qryCnt;
	prunings = prunings / (float)qryCnt;
	outfile << "==========KWQ平均结果==========\n";
	//outfile<<"IO代价 = "<< page_faults<<endl;
	//outfile<<"CPU代价 =  "<< ticksCnt<<endl;		
	outfile << "TotalTime = " << page_faults + ticksCnt << endl;
	outfile << "penalty = " << penaltys << endl;
	//outfile<<"剪枝效果 = "<<prunings<<endl;
	//outfile<<"hp resued = "<<hpUsed / qryCnt<<endl;
	outfile << "取样个数 = " << samplesize << endl;
	//outfile<<"剪枝个数 = "<<(int)( samplesize * prunings )<<endl;
	outfile << "========================================\n\n";
}

// 随机生成查询点，范围为min,max,并写入文件, min = 0 
void randQ(int dimension, int count, int min, int max)
{
	float* queryPoints = new float[dimension * count];

	for (int i = 0; i < dimension * count; i++)
	{
		int int_r = rand();
		long base = RAND_MAX - 1;
		float f_r = ((float)int_r) / base;
		queryPoints[i] = (max - min) * f_r + min;
	}
	char fname[100];
	memset(fname, '\0', 100);
	sprintf(fname, ".\\QueryPoints_%dd.txt", dimension);
	writeSet(queryPoints, dimension, count, fname);
}
void GetReady(vector<string> &datesets)
{ 
	//把权重文件名放入datasets中。
	ifstream infile("datasets_list.txt");
	if (!infile)
	{
		cout << "Open datasets_list.txt failed!" << endl;
		exit(1);
	}
	string file_name;
	string prefix(".\\dataset\\");
	while (getline(infile, file_name))
	{
		file_name = prefix + file_name;
		datesets.push_back(file_name);
	}
	for (int i = 0; i < datesets.size(); i++)
	{
		cout << '[' << i + 1 << ']' << '\t' << datesets[i] << endl;
	}
	cout << "Please select which dataset to search('0' to quit):";
    
	
}

//void Test()
//{
//	Cache *c = NULL;
//	RTree *rtree = NULL;
//	int dimension = 2;
//	char *fname = new char[MAX_FILE_NAME];
//	try{
//
//		memset(fname, 0, MAX_FILE_NAME);
//		sprintf(fname, ".\\dataset\\Test.txt");
//		c = new Cache(0, 1024);//分配内存
//		rtree = buildtree(fname, dimension, c);//对当前数据集建RTree，利用fname路径下的文件
//		cout << rtree->dimension << endl;
//		clock_t begin_tick, end_tick;
//		double duration = 0.0;//算法持续时间
//		Heap *hp = new Heap();
//		hp->init(dimension);//堆初始化
//
//		//---------------------begin--------------------------------------------
//		// 寻找查询点
//		//-----------------------------------------------------------------
//		float* queryPoints = new float[dimension * 1000];//查询点，每一列是一个查询点，dimension维
//		float* weightSet = new float[dimension * 1000];//缺失权重向量，每一列是一个权重向量，dimension维
//		int qryCnt = 1; // 满足条件的查询点个数 
//		int missWeightCnt = 1;// 丢失的向量个数
//		char setFileName[100];
//
//		memset(setFileName, '\0', 100);
//		sprintf(setFileName, ".\\querypoints\\4_1_Test.txt");//给查询点文件名赋值，查询点排名_缺失向量数_当前使用数据集名称
//
//		qryCnt = readSet(queryPoints, weightSet, setFileName, dimension, missWeightCnt);//读取querypoints的一个文件，把缺失向量放在weightSet里，把查询点放在queryPoints里
//		cout << "共有" << qryCnt << "个查询点" << endl;
//
//		float* missWeightSet = new float[missWeightCnt * dimension];//对当前数据集的单个查询点的一组丢失权重向量，可能是一个，可能是多个，根据数据集情况，每一行是一个权重
//		float *_qry_pt = new float[dimension];//对当前数据集的一个查询点
//
//		float *modifyweight = new float[dimension * missWeightCnt];//修改后的权重向量，每一列是一个权重
//		float  quality_of_answer = 0.5, probability_guarantee = 0.5;//计算抽样的空间大小S
//		int  k = 3, modifyk=k;
//		float pruning = 0.0, prunings = 0.0;
//		float penalty = 0.0, penaltys = 0.0;
//		//ofstream outfile;
//		float ticksCnt = 0, page_faults = 0;
//		long samplesize = (long)(log(1 - probability_guarantee) / log(1 - quality_of_answer / 100)) + 1;//计算抽样的空间大小S
//		cout << "抽样空间大小为" << samplesize << endl;
//		memset(setFileName, '\0', 100);
//		sprintf(setFileName, ".\\results\\%s", &fname[10]);//把输出结果文件名赋值给setFilename
//		outfile.open(setFileName, ios::out | ios::app);//如果没有文件，那么生成空文件；如果有文件，那么在文件尾追加。
//
//		//------------------------------begin-------------------------------------------------------------------
//		// 执行查看查询点是否是反Top-k结果
//		//---------------------------------------------------------------------	
//		//for (int cnt = 0; cnt < qryCnt; cnt++) //按照查询点个数进行循环
//		//{
//		//	outfile << "第\t" << cnt << "\t个查询点:\n";
//		//	memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//赋值给查询点  _qry_pt
//		//	memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //赋值给缺失权重
//		//	c->page_faults = 0;//cache c
//		//	/*计算当前查询点的主导点和不可比点*/
//		//	int ranking, rank;
//		//	int bedomnum, incomnum;  //incomnum是incomparable点的个数，bedomnum是dominate查询点q的点的个数，
//		//	float *incomparable = new float[3000000];
//		//	rtree->incomparable_points_reuse(hp, incomparable, incomnum, _qry_pt, bedomnum);   //求incomparable点的集合，放在数组incomparable中，并求出incomnum和 bedomnum的值
//		//	outfile << "主导点有" << bedomnum - 1 << "个" << endl;//因为主导点包含了查询点自己
//		//	outfile << "不可比点有" << incomnum << "个" << endl;
//		//	for (int j = 0; j < incomnum; j++)
//		//	{
//		//		for (int i = 0; i < dimension; i++)
//		//		{
//		//			outfile << incomparable[j*dimension + i] << "　";
//		//		}
//		//		outfile << endl;
//		//	}
//		//	float * currentWeight = new float[dimension];
//		//	for (int misweightCnt = 0; misweightCnt < missWeightCnt; misweightCnt++)
//		//	{
//		//		memcpy(currentWeight, &missWeightSet[misweightCnt * dimension], sizeof(float) * dimension);//赋值给查询点  _qry_pt
//		//		rtree->bbs_ranking_querypoint(incomparable, incomnum, _qry_pt, currentWeight, ranking);//此处传入currentWeight必须是单个权重
//		//		rank = ranking + bedomnum;
//		//		outfile << "当前rank=" << rank << "　" << "当前Topk值=" << k << endl;
//		//		if (rank <= k)
//		//		{
//		//			outfile << "==========是反Top-k结果，不用改正==========" << "rank=" << rank << endl;
//		//		}
//		//		else outfile << "==========不是反Top-k结果，需要改正==========" << "rank=" << rank << endl;
//		//	}
//		//}
//	
//		vector<qPoint> qvector;
//		cout<<"建好了vector"<<endl;
//		ifReverse(qvector, qryCnt, _qry_pt, queryPoints, missWeightSet, weightSet, missWeightCnt, dimension, c, rtree, hp, k);
//		cout<<"完成了反Top-k查询"<<endl;
//		cout<<"当前qvector大小为"<<qvector.size()<<endl;
//		/*for(int i=0;i<qryCnt;i++)
//		{
//			int missCnt=qvector[i].m_missCnt;
//			cout<<"当前为第"<<i<<"个查询点:	";
//			for(int j=0;j<dimension;j++)
//				cout<<qvector[i].m_qry_pt[j]<<"	";
//			cout<<endl<<"共有"<<missCnt<<"个缺失权重	"<<endl;
//			cout<<qvector[i].m_notReWeightcnt<<"个非反Top-k结果的权重:"<<endl;
//			for(int j=0;j<qvector[i].m_notReWeightcnt;j++)
//			{
//				for(int z=0;z<dimension;z++)
//				{
//					cout<<qvector[i].m_notReWeight[j*dimension+z]<<",";
//				}
//				cout<<endl;
//			}
//			cout<<qvector[i].m_reWeightcnt<<"个是反Top-k结果的权重:"<<endl;
//			for(int j=0;j<qvector[i].m_reWeightcnt;j++)
//			{
//				for(int z=0;z<dimension;z++)
//				{
//					cout<<qvector[i].m_reWeight[j*dimension+z]<<",";
//				}
//				cout<<endl;
//			}
//		}
//		cout<<"查看是否是反Top-k结果"<<endl;*/
//		
//		modifyKW(qvector,_qry_pt,fname, qryCnt, missWeightCnt, queryPoints, missWeightSet, weightSet, c, rtree, dimension, hp, modifyweight, samplesize, quality_of_answer, probability_guarantee, k, modifyk);
//		cout<<"kw修改完成"<<endl;
//		queue< float* > qmodifiedQueue;
//		modifyQ(qmodifiedQueue,qvector,_qry_pt,dimension,c,rtree,fname,qryCnt,missWeightCnt,queryPoints,missWeightSet,weightSet,hp,k);
//		cout<<"q修改完成"<<endl;
//		modifyKWQ(qmodifiedQueue,qvector,_qry_pt,dimension,c,rtree,fname,qryCnt,missWeightCnt,queryPoints,missWeightSet,weightSet,hp,k,samplesize,modifyweight,modifyk,quality_of_answer,probability_guarantee);
//		cout<<"kwq修改完成"<<endl;
//
//		//Sleep(100000);
//		//------------------------------begin-------------------------------------------------------------------
//		// 执行近似修改k w 的算法
//		//---------------------------------------------------------------------	
//		//outfile << "========== 修改 k w 的近似算法 ==========\n";
//		//outfile << "共有\t" << qryCnt << "\t个查询点\n";
//		//outfile << "数据集：" << fname << endl;
//		//outfile << endl;
//		//for (int cnt = 0; cnt < qryCnt; cnt++) //按照查询点个数进行循环
//		//{
//		//	outfile << "第\t" << cnt << "\t个查询点:\n";
//		//	memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//赋值给查询点  _qry_pt
//		//	memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //赋值给缺失权重
//		//	c->page_faults = 0;//cache c			
//		//			begin_tick = clock();
//		//			rtree->modifyWandK_approximate_reuse(hp, modifyweight, modifyk, penalty, _qry_pt, missWeightSet, missWeightCnt, k, quality_of_answer, probability_guarantee, pruning);//修改kw
//		//			end_tick = clock();
//		//			//将结果写入文件
//		//			outfile << "---------- 参数设置----------\n";
//		//			outfile << "查询点q：( ";
//		//			for (int i = 0; i < dimension - 1; i++)
//		//				outfile << _qry_pt[i] << " , ";
//		//			outfile << _qry_pt[dimension - 1] << " )" << endl;   //输出当前查询点的各个维度的值
//		//			outfile << "反top-k查询的k值： " << k << endl;  //原k值
//		//			outfile << "维度d = " << dimension << endl;
//		//			outfile << "丢失的向量集合大小|W| = " << missWeightCnt << endl; //缺失权重个数
//		//			outfile << "丢失的向量集合:" << endl;
//		//			for (int i = 0; i < missWeightCnt; i++)
//		//			{
//		//				//outfile << "在why-not向量下的rank = " << rank << " : ( ";
//		//				int j;
//		//				for (j = 0; j < dimension - 1; j++)
//		//					outfile << missWeightSet[i * dimension + j] << " , ";
//		//				outfile << missWeightSet[j] << " )" << endl;
//		//			}
//		//			outfile << "Pr = " << probability_guarantee << endl;
//		//			outfile << "T = " << quality_of_answer << endl;
//		//			outfile << "\n----------执行结果----------" << endl;
//		//			outfile << "向量集合:" << endl;
//		//			for (int i = 0; i < missWeightCnt; i++)
//		//			{
//		//				outfile << i << ". (";
//		//				int j;
//		//				for (j = 0; j < dimension - 1; j++)
//		//					outfile << modifyweight[i*dimension + j] << " , ";
//		//				outfile << modifyweight[j] << " )" << endl;
//		//			}
//		//			outfile << "k = " << modifyk << endl;
//		//			outfile << "IO代价 = " << c->page_faults / (float)100 << endl;
//		//			outfile << "CPU代价 =  " << (end_tick - begin_tick) / (float)1000 << endl;
//		//			outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
//		//			outfile << "penalty = " << penalty << endl;
//		//			outfile << "剪枝效果 = " << pruning << endl;
//		//			outfile << "剪枝个数 = " << (int)(samplesize * pruning) << endl;
//		//			ticksCnt += (end_tick - begin_tick);
//		//			page_faults += (c->page_faults);
//		//			prunings += pruning;
//		//			penaltys += penalty;
//		//		
//		//	
//		//	outfile << "--------------------\n";
//		//}
//		//ticksCnt = ticksCnt / (float)1000 / (float)qryCnt;
//		//page_faults = page_faults / (float)100 / (float)qryCnt;
//		//penaltys = penaltys / (float)qryCnt;
//		//prunings = prunings / (float)qryCnt;
//		//outfile << "==========KW平均结果==========\n";
//		//outfile << "IO代价 = " << page_faults << endl;
//		//outfile << "CPU代价 =  " << ticksCnt << endl;
//		//outfile << "TotalTime = " << page_faults + ticksCnt << endl;
//		//outfile << "penalty = " << penaltys << endl;
//		//outfile << "剪枝效果 = " << prunings << endl;
//		//outfile << "取样个数 = " << samplesize << endl;
//		//outfile << "剪枝个数 = " << (int)(samplesize * prunings) << endl;
//		//outfile << "========================================\n\n";
//		//Sleep(100000);
//		//-----------------------------------end-----------------------------------------------------------------------
//
//		//------------------------------begin-------------------------------------------------------------------
//		// 执行精确修改Q 的算法
//		//---------------------------------------------------------------------
//		//float* qmodified = new float[dimension];//修改后的q放在qmodified
//		////queue< float* > qmodifiedQueue;
//		//float qabsolute;
//		//ticksCnt = 0.0, page_faults = 0.0;
//		//penalty = 0.0, penaltys = 0.0;
//		//memset(setFileName, '\0', 100);
//		//if(outfile.is_open())
//		//	cout<<"open";
//		//outfile<<"========== 修改 q 的精确算法 ==========\n---------- 参数设置----------\n";
//		//outfile<<"共有\t"<<qryCnt<<"\t个查询点\n";
//		//outfile<<"数据集："<<fname<<endl;
//		//for (int cnt = 0; cnt < qryCnt; cnt++)
//		//{
//		//	outfile << "第\t" << cnt << "\t个查询点:\n";
//		//	memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//赋值给查询点
//		//	memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt);//赋值给缺失向量
//		//	c->page_faults = 0;			
//		//	begin_tick = clock();
//		//	rtree->modifyQ_accurate(hp, _qry_pt, missWeightSet, missWeightCnt, k, qmodified);//精确修改Q
//		//	end_tick = clock();
//		//	qmodifiedQueue.push(qmodified);
//		//	outfile << "查询点q：( ";
//		//	for (int i = 0; i < dimension - 1; i++)
//		//		outfile << _qry_pt[i] << " , ";
//		//	outfile << _qry_pt[dimension - 1] << " )" << endl;
//		//	outfile << "反top-k查询的k值： " << k << endl;
//		//	outfile << "维度d = " << dimension << endl;
//		//	outfile << "丢失的向量集合大小|W| = " << missWeightCnt << endl;
//		//	outfile << "丢失的向量集合:" << endl;
//		//	for (int i = 0; i < missWeightCnt; i++)
//		//	{
//		//		//outfile << "在why-not向量下的rank = " << rank << " : ( ";
//		//		int j;
//		//		for (j = 0; j < dimension - 1; j++)
//		//			outfile << missWeightSet[i * dimension + j] << " , ";
//		//		outfile << missWeightSet[j] << " )" << endl;
//		//	}
//		//	outfile << "\n----------执行结果----------" << endl;
//		//	outfile << "查询点修改: (";
//		//	qabsolute = penalty = 0.0;
//		//	for (int i = 0; i < dimension - 1; i++)
//		//	{
//		//		outfile << qmodified[i] << " , ";
//		//		penalty += (qmodified[i] - _qry_pt[i]) * (qmodified[i] - _qry_pt[i]);
//		//		qabsolute += _qry_pt[i] * _qry_pt[i];
//		//	}
//		//	outfile << qmodified[dimension - 1] << " )" << endl;
//		//	penalty += (qmodified[dimension - 1] - _qry_pt[dimension - 1]) * (qmodified[dimension - 1] - _qry_pt[dimension - 1]);
//		//	qabsolute += _qry_pt[dimension - 1] * _qry_pt[dimension - 1];
//		//	penalty = penalty / qabsolute;
//		//	outfile << "IO代价 = " << c->page_faults / (float)100 << endl;
//		//	outfile << "CPU代价 =  " << (end_tick - begin_tick) / (float)1000 << endl;
//		//	outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
//		//	outfile << "penalty = " << penalty << endl;
//		//	ticksCnt += (end_tick - begin_tick);
//		//	page_faults += (c->page_faults);
//		//	penaltys += penalty;
//		//	outfile << "--------------------\n";			
//		//}
//		//ticksCnt  =  ticksCnt / ( float) 1000 / ( float ) qryCnt;
//		//page_faults = page_faults / ( float ) 100 / ( float ) qryCnt;
//		//penaltys = penaltys / ( float ) qryCnt;
//		//outfile<<"==========AccurateQ平均结果==========\n";
//		//outfile<<"IO代价 = "<< page_faults<<endl;
//		//outfile<<"CPU代价 =  "<< ticksCnt<<endl;		
//		//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
//		//outfile<<"penalty = "<<penaltys <<endl;
//		//outfile<<"========================================\n\n";
//		//-----------------------------------end-----------------------------------------------------------------------
//
//
////------------------------------begin-------------------------------------------------------------------
//// 执行修改K W Q 的算法
////---------------------------------------------------------------------
//			//float* qmodified_kwq = new float[dimension];
//			//ticksCnt = 0.0, page_faults = 0.0;
//			//pruning = 0.0, prunings = 0.0;
//			//penalty = 0.0, penaltys = 0.0;
//			//memset(setFileName, '\0', 100);
//			//if(outfile.is_open())
//			//	cout<<"open";
//			//outfile<<"========== 修改 q k w 的算法 ==========\n---------- 参数设置----------\n";
//			//outfile<<"共有\t"<<qryCnt<<"\t个查询点\n";
//			//outfile<<"数据集："<<fname<<endl;
//			//int hpUsed = 0;
//			//for( int cnt = 0; cnt < qryCnt; cnt++ )
//			//{
//			//	outfile<<"第\t"<<cnt<<"\t个查询点:\n";
//			//	memcpy( _qry_pt,  &queryPoints[cnt * dimension], sizeof( float ) * dimension );
//			//	memcpy( missWeightSet,  &weightSet[cnt * dimension * missWeightCnt], sizeof( float ) * dimension * missWeightCnt);
//			//	pruning = 0.0;
//			//	qmodified = qmodifiedQueue.front();
//			//	qmodifiedQueue.pop();
//			//	c->page_faults = 0;
//			//	hp->used = 0;
//			//	begin_tick = clock();
//			//	cout << "xiugaiqian" << endl;
//			//	rtree->modify_k_W_Q_reuse(hp, _qry_pt,  qmodified, missWeightSet, missWeightCnt, k, qmodified_kwq, modifyweight, modifyk, quality_of_answer, probability_guarantee, penalty, pruning);
//			//	end_tick = clock();
//			//	cout << "paozhene" << endl;
//			//	 //将结果写入文件
//			//	outfile<<"查询点q：( ";
//			//	for( int i = 0; i < dimension - 1; i++ )
//			//		outfile<<_qry_pt[i]<<" , ";
//			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;
//			//	hpUsed += hp->used;
//			//	outfile<<"修改Q方案得到的查询点: ( ";
//			//	for( int i = 0; i < dimension - 1; i++ )
//			//		outfile<<qmodified[i]<<" , ";
//			//	outfile<<qmodified[dimension - 1]<<" )"<<endl;
//			//	outfile<<"反top-k查询的k值： "<<k<<endl;
//			//	outfile<<"维度d = "<<dimension<<endl;
//			//	outfile<<"丢失的向量集合大小|W| = "<<missWeightCnt<<endl;
//			//	outfile<<"丢失的向量集合:"<<endl;
//			//	for( int i = 0; i < missWeightCnt; i++ )
//			//	{
//			//		//outfile<<"在why-not向量下的rank = "<<*(qry_ranks_p - i )<<" : ( ";
//			//		int j;
//			//		for( j = 0; j < dimension - 1; j++ )
//			//			outfile<<missWeightSet[i*dimension + j]<<" , ";
//			//		outfile<<missWeightSet[j]<<" )"<<endl;
//			//	}
//			//	outfile<<"T = "<<quality_of_answer<<endl;
//			//	outfile<<"Pr = "<<probability_guarantee<<endl;
//			//	outfile<<"\n----------执行结果----------"<<endl;
//			//	outfile<<"向量集合:"<<endl;
//			//	for( int i = 0; i < missWeightCnt; i++ )
//			//	{
//			//		outfile<<i<<". (";
//			//		int j;
//			//		for( j = 0; j < dimension - 1; j++ )
//			//			outfile<<modifyweight[i*dimension + j]<<" , ";
//			//		outfile<<modifyweight[j]<<" )"<<endl;
//			//	}
//			//	outfile<<"k = "<<modifyk<<endl;
//			//	outfile<<"修改查询点q：( ";
//			//	for( int i = 0; i < dimension - 1; i++ )
//			//		outfile<<qmodified_kwq[i]<<" , ";
//			//	outfile<<qmodified_kwq[dimension - 1]<<" )"<<endl;
//			//	outfile<<"IO代价 = "<< c->page_faults / ( float )100<<endl;
//			//	outfile<<"CPU代价 =  "<< ( end_tick - begin_tick ) / ( float )1000<<endl;		
//			//	outfile<<"TotalTime = "<< c -> page_faults / ( float )100 + ( end_tick - begin_tick ) / ( float )1000<<endl;
//			//	outfile<<"penalty = "<<penalty<<endl;
//			//	outfile<<"剪枝效果 = "<<pruning<<endl;
//			//	outfile<<"hp resued = "<<hp->used<<endl;
//			//	outfile<<"剪枝个数 = "<<(int)(samplesize * pruning)<<endl;
//			//	ticksCnt += ( end_tick - begin_tick );
//			//	page_faults += ( c->page_faults );
//			//	prunings += pruning;
//			//	penaltys += penalty;
//			//	outfile<<"--------------------\n";
//			//}
//			//ticksCnt  =  ticksCnt / ( float) 1000 / ( float ) qryCnt;
//			//page_faults = page_faults / ( float ) 100 / ( float ) qryCnt;
//			//penaltys = penaltys / ( float ) qryCnt;
//			//prunings = prunings / ( float ) qryCnt;
//			//outfile<<"==========KWQ平均结果==========\n";
//			//outfile<<"IO代价 = "<< page_faults<<endl;
//			//outfile<<"CPU代价 =  "<< ticksCnt<<endl;		
//			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
//			//outfile<<"penalty = "<<penaltys <<endl;
//			//outfile<<"剪枝效果 = "<<prunings<<endl;
//			//outfile<<"hp resued = "<<hpUsed / qryCnt<<endl;
//			//outfile<<"取样个数 = "<<samplesize<<endl;
//			//outfile<<"剪枝个数 = "<<(int)( samplesize * prunings )<<endl;
//			//outfile<<"========================================\n\n";
//			//-----------------------------------end-----------------------------------------------------------------------		
//
//
//
//		outfile.close();
//		//---------------------------------------------
//
//		//---------------------------------------------
//		if (rtree)
//		{
//			delete rtree;
//			rtree = NULL;
//		}
//	}
//	catch (exception e)
//	{
//		if (rtree)
//		{
//			delete rtree;
//			rtree = NULL;
//		}
//	}
//	delete c;
//	delete[] fname;
//	
//}
int main()
{

	//Test();
	// do_test();
	RTree *rtree = NULL;
	Cache *c = NULL;
	char *fname = new char[MAX_FILE_NAME];
	int query_number = 0;
	int dateset_no = 0;
	int dimension = 2;
	int k = 10;
	int queryRank = 101;
	int missWeightCnt = 1;
	float  quality_of_answer = 0.6;
	float  probability_guarantee = 0.4;
	int parameterCnt = 1;//权重向量属性个数,默认为1
	
	vector<string> datesets;
	unordered_map<int, attribute> weight_table;
	unordered_map<int, attribute> attri_table1;
	unordered_map<int, attribute> attri_table2;
	GetReady(datesets);
	// get a list of datasets.Content of GetReady()
	/*ifstream infile("datasets_list.txt");
	if (!infile)
	{
	cout << "Open datasets_list.txt failed!" << endl;
	exit(1);
	}

	vector<string> datesets;
	string file_name;
	string prefix(".\\dataset\\");
	while (getline(infile, file_name))
	{
	file_name = prefix + file_name;
	datesets.push_back(file_name);
	}

	for (i = 0; i < datesets.size(); i++)
	{
	cout << '[' << i+1 << ']' << '\t' << datesets[i] << endl;
	}
	cout << "Please select which dataset to search('0' to quit):" ;*/
	
	int dateset_nos[] = {
		3, 3, 3, 3, 3,             // |W| = 1, 2, 3, 4, 5             每个q的rank 相同
		/*	13, 13, 13, 13, 13,
		23, 23, 23, 23, 23,
		33, 33, 33, 33, 33,*/
		2, 4, 5,                  // 维度
		//12, 14, 15,          //5
		//22, 24, 25,            //5
		//32, 34, 35,                                 //1
		6, 7, 10, 11,			//数据集大小
		/*16, 17, 20, 21,
		26, 27, 30, 31,
		36, 37, 40, 41,		*/
		3, 3, 3, 					// q 在why not 下的 rank
		/*13, 13, 13,
		23, 23, 23,
		33, 33, 33,		*/
		3, 3, 3, 3, 					//k
		/*	13, 13, 13, 13,
		23, 23, 23, 23,
		33, 33, 33, 33,*/
		3, 3, 3, 3,   				//5  ts
		/*13, 13, 13, 13,
		23, 23, 23, 23,
		33, 33, 33, 33, */
		3, 3, 3, 3,   				//5  pr
		/*13, 13, 13, 13,
		23, 23, 23, 23,
		33, 33, 33, 33*/
	};

	int ks[] = {
		10, 10, 10, 10, 10,                 // |W| = 1, 2, 3, 4, 5             每个q的rank 相同
		/*10, 10, 10,	10, 10,
		10, 10, 10,	10, 10,
		10, 10, 10,	10, 10,*/

		//=======================================================

		//============for real data begin=======================
		10, 10, 10,             // 维度 
		/*10, 10, 10,
		10, 10, 10,
		10, 10, 10,*/
		10, 10, 10, 10,             //数据集大小
		/*10, 10, 10, 10,
		10, 10, 10, 10,
		10, 10, 10, 10, */
		//============for real data end=======================
		10, 10, 10,          // q 在why not 下的 rank
		/*10, 10, 10,
		10, 10, 10,
		10, 10, 10, */
		20, 30, 40, 50,  // k
		/*20, 30, 40, 50,
		20, 30, 40, 50,
		20, 30, 40, 50,*/
		10, 10, 10, 10,//ts
		/*10, 10, 10, 10,
		10, 10, 10, 10,
		10, 10, 10, 10,*/
		10, 10, 10, 10,  //pr
		/*10, 10, 10, 10,
		10, 10, 10, 10,
		10, 10, 10, 10*/
	};// 反top-k 的k 值
	int* ks_p = ks;
	int qry_ranks[] = {
		101,                             // |W| = 1, 2, 3, 4, 5             每个q的rank 相同
		101, 101,
		101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101, 101,
		/*101,
		101, 101,
		101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101, 101,
		101,
		101, 101,
		101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101, 101,
		101,
		101, 101,
		101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101, 101,*/
		////////////////////////////////

		//============for real data begin=======================
		101, 101, 101,                     // 维度
		/*101, 101, 101,
		101, 101, 101,
		101, 101, 101,*/
		101, 101, 101, 101,          // 数据集大小
		/*101, 101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101, */
		//============for real data end=======================
		11, 501, 1001,               // q 在why not 下的 rank
		/* 11, 501, 1001,
		11, 501, 1001,
		11, 501, 1001,*/
		201, 301, 401, 501,          //k
		/*201, 301, 401, 501,
		201, 301, 401, 501,
		201, 301, 401, 501, */
		101, 101, 101, 101,      // ts
		/*101, 101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101, */
		101, 101, 101, 101,   // pr
		/*101, 101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101*/
	};// 查询点在丢失向量top-k查询中的rank
	int* qry_ranks_p = qry_ranks;
	float Prs[] = {
		0.4, 0.4, 0.4, 0.4, 0.4,       // |W| = 1, 2, 3, 4, 5             每个q的rank 相同
		/*0.8, 0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8, 0.8,*/
		//=======================================================

		//============for real data begin=======================
		0.8, 0.8, 0.8,            // 维度 
		/*0.8, 0.8, 0.8,
		0.8, 0.8, 0.8,
		0.8, 0.8, 0.8,*/
		0.8, 0.8, 0.8, 0.8,         // 数据集大小
		/*0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,*/
		//============for real data end=======================
		0.8, 0.8, 0.8,            // q 在why not 下的 rank
		/*0.8, 0.8, 0.8,
		0.8, 0.8, 0.8,  */
		0.8, 0.8, 0.8, 0.8,         // k
		/*0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,  */
		0.8, 0.8, 0.8, 0.8,                 //测试 Ts
		/*0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,*/
		0.5, 0.6, 0.7, 0.9,  // 测试pr
		/*0.5, 0.6, 0.7, 0.9,
		0.5, 0.6, 0.7, 0.9,
		0.5, 0.6, 0.7, 0.9*/
	};
	float *Prs_p = Prs;
	float Ts[] = {
		0.6, 0.6, 0.6, 0.6, 0.6,         // |W| = 1, 2, 3, 4, 5             每个q的rank 相同
		/*0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2, */
		//=======================================================
		//============for real data begin=======================
		0.2, 0.2, 0.2,            // 维度
		/*0.2, 0.2, 0.2,
		0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, */
		0.2, 0.2, 0.2, 0.2,            // 数据集大小
		/*	0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,*/
		//============for real data begin=======================
		0.2, 0.2, 0.2,                // q 在why not 下的 rank
		/*0.2, 0.2, 0.2,
		0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, */
		0.2, 0.2, 0.2, 0.2,                 // k
		/*0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, */
		0.1, 0.15, 0.25, 0.3,   // 测试 ts
		/*0.1, 0.15, 0.25, 0.3,
		0.1, 0.15, 0.25, 0.3,
		0.1, 0.15, 0.25, 0.3,   */
		0.2, 0.2, 0.2, 0.2,       // 测试 pr
		/*0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2*/
	};
	float *Ts_p = Ts;
	int abW[] = {
		1, 2, 3, 4, 5,                                    // |W| 每个 q的rank 不同
		//======================
		/*1, 2, 3, 4, 5,
		1, 2, 3, 4, 5,
		1, 2, 3, 4, 5,*/
		//======================
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1,
		1, 1, 1, 1, 1
	};
	int* abW_p = abW;
	int paCnt[] = {          //权重属性个数，默认为1
		1, 1, 1, 1,1,1
	};
	int * pa_Cnt = paCnt;
	//queue<int> qryPosTmp;  // 符合条件的查询点在查询点向量集合中的起始位置
	//queue< float* > weightSetStackTmp; // 丢失的向量集合

	for (int dataset_pos = 0; dataset_pos < 1/*27*/; dataset_pos++)//循环次数
	{
		missWeightCnt = 5;
		try{
			//================
			//dateset_no = 3;
			//==============================
			dateset_no = dateset_nos[dataset_pos]; //每一次循环的文件名
			//====================================			
			//cin >> dateset_no;
			//while (dateset_no > 0 && dateset_no <= datesets.size())
			//{
			memset(fname, 0, MAX_FILE_NAME);
			memcpy(fname, datesets[dateset_no - 1].c_str(), datesets[dateset_no - 1].size() + 1);//把datesets[dateset_no-1].size()+1个datesets[dateset_no-1].c_str()中的字符打印复制到fname，即当前数据集名称
			//strcpy(fname, dateset[dateset_no-1]);

			if (datesets[dateset_no - 1].find("_2d_") != string::npos)//datesets[dateset_no-1]为当前数据集文件名，判断维度
				dimension = 2;
			else if (datesets[dateset_no - 1].find("_3d_") != string::npos)
				dimension = 3;
			else if (datesets[dateset_no - 1].find("_4d_") != string::npos)
				dimension = 4;
			else if (datesets[dateset_no - 1].find("_5d_") != string::npos)
				dimension = 5;
			else if (datesets[dateset_no - 1].find("_6d_") != string::npos)
				dimension = 6;
			else if (datesets[dateset_no - 1].find("_8d_") != string::npos)
				dimension = 8;
			else if (datesets[dateset_no - 1].find("_10d_") != string::npos)
				 
				dimension = 10;
			else if (datesets[dateset_no - 1].find("_13d_") != string::npos)
				dimension = 13;

			c = new Cache(0, 1024);//分配内存
			rtree = buildtree(fname, dimension, c);//对当前数据集建RTree，利用fname路径下的文件
			cout << rtree->dimension << endl;


			/*while (query_number != -1)
			{*/
			clock_t begin_tick, end_tick;
			double duration = 0.0;//算法持续时间

			//float *rslt = new float[2 * MAX_SKYLINE_RESULT];
			//int rsltcnt = 0;
			Heap *hp = new Heap();
			hp->init(dimension);//堆初始化
			//string which_method;//什么方法

			//////////////////////////生成向量集合，并写入文件//////////////////////////////////////////////////////////////////////
			//for( int d = 2; d <= 5; d++ )
			//{
			//	float* weights = rtree->SampleWeights(d, 1000, true);
			//	char fname[100];
			//	memset(fname, '\0', 100);
			//	sprintf(fname,".\\SampleWeights_%dd.txt", d );
			//	writeSet(weights, d, 1000, fname); 
			//}


			//---------------------begin--------------------------------------------
			// 寻找查询点
			//-----------------------------------------------------------------
			
			float* queryPoints = new float[dimension * 1000];//查询点，每一列是一个查询点，dimension维
			float* weightSet = new float[dimension * 1000];//缺失权重向量，每一列是一个权重向量，dimension维
			int* missIdSet = new int[parameterCnt * 1000];
			int qryCnt = 0; // 满足条件的查询点个数 
			//missWeightCnt = *abW_p;// 丢失的向量个数

			//char setFileName[100];

			memset(setFileName, '\0', 100);
			sprintf(setFileName, ".\\querypoints\\%d_%d_%s", *qry_ranks_p, missWeightCnt, &fname[10]);//给查询点文件名赋值，查询点排名_缺失向量数_当前使用数据集名称

			qryCnt = readSet(queryPoints, weightSet, setFileName, dimension, missWeightCnt,missIdSet);//读取querypoints的一个文件，把缺失向量放在weightSet里，把查询点放在queryPoints里
			float* missWeightSet = new float[missWeightCnt * dimension];//对当前数据集的单个查询点的一组丢失权重向量，可能是一个，可能是多个，根据数据集情况，每一行是一个权重
			float *_qry_pt = new float[dimension];//对当前数据集的一个查询点
			//---------------------end--------------------------------------------


			//int qry_pos = 163;
			//for(; qry_pos < 164 /*&& qryCnt <= 10*/; qry_pos++ )
			//{
			//	cout<<"q: "<<qry_pos<<endl;
			//	missWeightSet = rtree->checkQueryPoint ( hp, &queryPoints[ qry_pos * dimension ], weightSet, 1000, dimension, missWeightCnt, qry_ranks_p, &fname[10] ) ;
			//	if( missWeightSet != NULL )
			//	{
			//		qryCnt++;
			//		//qryPosQueue.push( qry_pos );
			//		//weightSetQueue.push( missWeightSet );
			//		//qry_ranks_p = qry_ranks_p + *abW_p;// 当在这个数据集找到了合适的查询点，那么把指针指向下一个查询点值
			//		abW_p = abW_p + 1;
			//		ks_p = ks_p + 1;
			//		Prs_p = Prs_p + 1;
			//		Ts_p = Ts_p + 1;
			//		//memcpy( _qry_pt, &queryPoints[qry_pos * dimension], sizeof(float) * dimension );
			//		cout<<"queryPoints";
			//		for( int d = 0; d < dimension; d++ )
			//			cout<<queryPoints[qry_pos * dimension + d]<<",\t";
			//		cout<<endl;
			//		
			//	}
			//}
			//// 这个集合中没有找到相应的why -not 向量，继续下一个
			//if( qry_pos == 1000 )
			//{
			//	abW_p = abW_p + 1;
			//	Prs_p = Prs_p + 1;
			//	Ts_p = Ts_p + 1;
			//	ks_p = ks_p + 1;
			//	qry_ranks_p = qry_ranks_p + 1;
			//	continue;
			//}
			//--------------------------end-----------------------------------------


			float *modifyweight = new float[dimension * missWeightCnt];//修改后的权重向量，每一列是一个权重
			int  modifyk;
			float pruning = 0.0, prunings = 0.0;
			float penalty = 0.0, penaltys = 0.0;
			ofstream outfile;
			float ticksCnt = 0, page_faults = 0;
			long samplesize = (long)(log(1 - probability_guarantee) / log(1 - quality_of_answer / 100)) + 1;//计算抽样的空间大小S
			memset(setFileName, '\0', 100);
			sprintf(setFileName, ".\\results\\%d_%d_%s",queryRank,missWeightCnt, &fname[10]);//把输出结果文件名赋值给setFilename
			outfile.open(setFileName, ios::out | ios::app);//如果没有文件，那么生成空文件；如果有文件，那么在文件尾追加。

			vector<qPoint> qvector;
			cout << "建好了vector" << endl;
			ifReverse(qvector, qryCnt, queryPoints, weightSet, missWeightCnt, dimension, c, rtree, hp, k, missIdSet);
			cout << "完成了反Top-k查询" << endl;
			cout << "当前qvector大小为" << qvector.size() << endl;
			
			modifyKW(qvector, fname, qryCnt, missWeightCnt, queryPoints, weightSet, c, rtree, dimension, hp, samplesize, quality_of_answer, probability_guarantee, k, modifyk);
			cout << "kw修改完成" << endl;
			queue< float* > qmodifiedQueue;
			modifyQ(qmodifiedQueue, qvector, dimension, c, rtree, fname, qryCnt, missWeightCnt, queryPoints, weightSet, hp, k);
			cout << "q修改完成" << endl;
			modifyKWQ(qmodifiedQueue, qvector, dimension, c, rtree, fname, qryCnt, missWeightCnt, queryPoints, weightSet, hp, k, samplesize, modifyweight, modifyk, quality_of_answer, probability_guarantee);
			cout << "kwq修改完成" << endl;


			//------------------------------begin-------------------------------------------------------------------
			// 执行近似修改k w 的算法
			//---------------------------------------------------------------------	
			//outfile<<"========== 修改 k w 的近似算法 ==========\n";
			//outfile<<"共有\t"<<qryCnt<<"\t个查询点\n";
			//outfile<<"数据集："<<fname<<endl;
			//for( int cnt = 0; cnt < qryCnt; cnt++ ) //按照查询点个数进行循环
			//{
			//	outfile<<"第\t"<<cnt<<"\t个查询点:\n";
			//	memcpy( _qry_pt,  &queryPoints[cnt * dimension], sizeof( float ) * dimension );//赋值给查询点  _qry_pt
			//	memcpy( missWeightSet,  &weightSet[cnt * dimension * missWeightCnt], sizeof( float ) * dimension * missWeightCnt); //赋值给缺失权重
			//	c->page_faults = 0;//cache c
			//	begin_tick = clock();
			//	rtree->modifyWandK_approximate_reuse(hp, modifyweight, modifyk, penalty, _qry_pt,  missWeightSet, missWeightCnt, k,  quality_of_answer, probability_guarantee, pruning);//修改kw
			//	end_tick = clock();
			//	//将结果写入文件
			//	outfile<<"---------- 参数设置----------\n";
			//	outfile<<"查询点q：( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<_qry_pt[i]<<" , ";
			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;   //输出当前查询点的各个维度的值
			//	outfile<<"反top-k查询的k值： "<<k<<endl;  //原k值
			//	outfile<<"维度d = "<<dimension<<endl;
			//	outfile<<"丢失的向量集合大小|W| = "<<missWeightCnt<<endl; //缺失权重个数
			//	outfile<<"丢失的向量集合:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ ) 
			//	{
			//		outfile<<"在why-not向量下的rank = "<<*(qry_ranks_p - i )<<" : ( ";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<missWeightSet[ i * dimension + j]<<" , ";
			//		outfile<<missWeightSet[j]<<" )"<<endl;
			//	}
			//	outfile<<"Pr = "<<probability_guarantee<<endl;
			//	outfile<<"T = "<<quality_of_answer<<endl;
			//	outfile<<"\n----------执行结果----------"<<endl;
			//	outfile<<"向量集合:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<i<<". (";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<modifyweight[i*dimension + j]<<" , ";
			//		outfile<<modifyweight[j]<<" )"<<endl;
			//	}
			//	outfile<<"k = "<<modifyk<<endl;
			//	outfile<<"IO代价 = "<< c->page_faults / ( float ) 100<<endl;
			//	outfile<<"CPU代价 =  "<<( end_tick - begin_tick ) / ( float )1000<<endl;		
			//	outfile<<"TotalTime = "<<c->page_faults /( float ) 100 +( end_tick - begin_tick ) / ( float )1000<<endl;
			//	outfile<<"penalty = "<<penalty<<endl;
			//	outfile<<"剪枝效果 = "<<pruning<<endl;
			//	outfile<<"剪枝个数 = "<<(int)(samplesize * pruning)<<endl;
			//	ticksCnt += ( end_tick - begin_tick );
			//	page_faults += ( c->page_faults );
			//	prunings += pruning;
			//	penaltys += penalty;
			//	outfile<<"--------------------\n";
			//}
			//ticksCnt  =  ticksCnt / ( float) 1000 / ( float ) qryCnt;
			//page_faults = page_faults / ( float ) 100 / ( float ) qryCnt;
			//penaltys = penaltys / ( float ) qryCnt;
			//prunings = prunings / ( float ) qryCnt;
			//outfile<<"==========KW平均结果==========\n";
			//outfile<<"IO代价 = "<< page_faults<<endl;
			//outfile<<"CPU代价 =  "<< ticksCnt<<endl;		
			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
			//outfile<<"penalty = "<<penaltys <<endl;
			//outfile<<"剪枝效果 = "<<prunings<<endl;
			//outfile<<"取样个数 = "<<samplesize<<endl;
			//outfile<<"剪枝个数 = "<<(int)( samplesize * prunings )<<endl;
			//outfile<<"========================================\n\n";
			//-----------------------------------end-----------------------------------------------------------------------

			//------------------------------begin-------------------------------------------------------------------
			// 执行精确修改Q 的算法
			//---------------------------------------------------------------------
			//float* qmodified = new float[dimension];
			//queue< float* > qmodifiedQueue;
			//float qabsolute;
			//ticksCnt = 0.0, page_faults = 0.0;
			//penalty = 0.0, penaltys = 0.0;
			//memset(setFileName, '\0', 100);
			//if(outfile.is_open())
			//	cout<<"open";
			//outfile<<"========== 修改 q 的精确算法 ==========\n---------- 参数设置----------\n";
			//outfile<<"共有\t"<<qryCnt<<"\t个查询点\n";
			//outfile<<"数据集："<<fname<<endl;
			//for( int cnt = 0; cnt < qryCnt; cnt++ )
			//{
			//	outfile<<"第\t"<<cnt<<"\t个查询点:\n";
			//	memcpy( _qry_pt,  &queryPoints[cnt * dimension], sizeof( float ) * dimension );//赋值给查询点
			//	memcpy( missWeightSet,  &weightSet[cnt * dimension * missWeightCnt], sizeof( float ) * dimension * missWeightCnt);//赋值给缺失向量
			//	c->page_faults = 0.0;
			//	begin_tick = clock();
			//	rtree->modifyQ_accurate(hp, _qry_pt,  missWeightSet, missWeightCnt, k, qmodified);//精确修改Q
			//	end_tick = clock();
			//	qmodifiedQueue.push( qmodified );			
			//	outfile<<"查询点q：( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<_qry_pt[i]<<" , ";
			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;
			//	outfile<<"反top-k查询的k值： "<<k<<endl;
			//	outfile<<"维度d = "<<dimension<<endl;
			//	outfile<<"丢失的向量集合大小|W| = "<<missWeightCnt<<endl;
			//	outfile<<"丢失的向量集合:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<"在why-not向量下的rank = "<<* ( qry_ranks_p - i )<<" : ( ";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<missWeightSet[ i * dimension + j ]<<" , ";
			//		outfile<<missWeightSet[j]<<" )"<<endl;
			//	}
			//	outfile<<"\n----------执行结果----------"<<endl;
			//	outfile<<"查询点修改: (";
			//	
			//	qabsolute = penalty = 0.0;
			//	for( int i = 0; i < dimension - 1; i++ )
			//	{
			//		outfile<<qmodified[i]<<" , ";
			//		penalty += (qmodified[i] - _qry_pt[i]) * (qmodified[i] - _qry_pt[i]);
			//		qabsolute += _qry_pt[i] * _qry_pt[i];
			//	}
			//	outfile<<qmodified[dimension - 1]<<" )"<<endl;
			//	penalty += (qmodified[dimension - 1] - _qry_pt[dimension - 1]) * (qmodified[dimension - 1] - _qry_pt[dimension - 1]);
			//	qabsolute += _qry_pt[dimension - 1] * _qry_pt[dimension - 1];
			//	penalty = penalty / qabsolute;
			//	outfile<<"IO代价 = "<< c->page_faults / ( float )100<<endl;
			//	outfile<<"CPU代价 =  "<<( end_tick - begin_tick ) / ( float ) 1000<<endl;		
			//	outfile<<"TotalTime = "<< c->page_faults / ( float ) 100 +( end_tick - begin_tick ) / ( float ) 1000<<endl;
			//	outfile<<"penalty = "<<penalty <<endl;
			//	ticksCnt += ( end_tick - begin_tick );
			//	page_faults += ( c->page_faults );
			//	penaltys += penalty;
			//	outfile<<"--------------------\n";
			//}
			//ticksCnt  =  ticksCnt / ( float) 1000 / ( float ) qryCnt;
			//page_faults = page_faults / ( float ) 100 / ( float ) qryCnt;
			//penaltys = penaltys / ( float ) qryCnt;
			//outfile<<"==========AccurateQ平均结果==========\n";
			//outfile<<"IO代价 = "<< page_faults<<endl;
			//outfile<<"CPU代价 =  "<< ticksCnt<<endl;		
			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
			//outfile<<"penalty = "<<penaltys <<endl;
			//outfile<<"========================================\n\n";
			//-----------------------------------end-----------------------------------------------------------------------

			//------------------------------begin-------------------------------------------------------------------
			// 执行修改K W Q 的算法
			//---------------------------------------------------------------------
			//float* qmodified_kwq = new float[dimension];
			//ticksCnt = 0.0, page_faults = 0.0;
			//pruning = 0.0, prunings = 0.0;
			//penalty = 0.0, penaltys = 0.0;
			//memset(setFileName, '\0', 100);
			//if(outfile.is_open())
			//	cout<<"open";
			//outfile<<"========== 修改 q k w 的算法 ==========\n---------- 参数设置----------\n";
			//outfile<<"共有\t"<<qryCnt<<"\t个查询点\n";
			//outfile<<"数据集："<<fname<<endl;
			//int hpUsed = 0;
			//for( int cnt = 0; cnt < qryCnt; cnt++ )
			//{
			//	outfile<<"第\t"<<cnt<<"\t个查询点:\n";
			//	memcpy( _qry_pt,  &queryPoints[cnt * dimension], sizeof( float ) * dimension );
			//	memcpy( missWeightSet,  &weightSet[cnt * dimension * missWeightCnt], sizeof( float ) * dimension * missWeightCnt);
			//	pruning = 0.0;
			//	qmodified = qmodifiedQueue.front();
			//	qmodifiedQueue.pop();
			//	c->page_faults = 0;
			//	hp->used = 0;
			//	begin_tick = clock();
			//	cout << "xiugaiqian" << endl;
			//	rtree->modify_k_W_Q_reuse(hp, _qry_pt,  qmodified, missWeightSet, missWeightCnt, k, qmodified_kwq, modifyweight, modifyk, quality_of_answer, probability_guarantee, penalty, pruning);
			//	end_tick = clock();
			//	cout << "paozhene" << endl;
			//	 //将结果写入文件
			//	outfile<<"查询点q：( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<_qry_pt[i]<<" , ";
			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;
			//	hpUsed += hp->used;
			//	outfile<<"修改Q方案得到的查询点: ( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<qmodified[i]<<" , ";
			//	outfile<<qmodified[dimension - 1]<<" )"<<endl;
			//	outfile<<"反top-k查询的k值： "<<k<<endl;
			//	outfile<<"维度d = "<<dimension<<endl;
			//	outfile<<"丢失的向量集合大小|W| = "<<missWeightCnt<<endl;
			//	outfile<<"丢失的向量集合:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<"在why-not向量下的rank = "<<*(qry_ranks_p - i )<<" : ( ";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<missWeightSet[i*dimension + j]<<" , ";
			//		outfile<<missWeightSet[j]<<" )"<<endl;
			//	}
			//	outfile<<"T = "<<quality_of_answer<<endl;
			//	outfile<<"Pr = "<<probability_guarantee<<endl;
			//	outfile<<"\n----------执行结果----------"<<endl;
			//	outfile<<"向量集合:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<i<<". (";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<modifyweight[i*dimension + j]<<" , ";
			//		outfile<<modifyweight[j]<<" )"<<endl;
			//	}
			//	outfile<<"k = "<<modifyk<<endl;
			//	outfile<<"修改查询点q：( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<qmodified_kwq[i]<<" , ";
			//	outfile<<qmodified_kwq[dimension - 1]<<" )"<<endl;
			//	outfile<<"IO代价 = "<< c->page_faults / ( float )100<<endl;
			//	outfile<<"CPU代价 =  "<< ( end_tick - begin_tick ) / ( float )1000<<endl;		
			//	outfile<<"TotalTime = "<< c -> page_faults / ( float )100 + ( end_tick - begin_tick ) / ( float )1000<<endl;
			//	outfile<<"penalty = "<<penalty<<endl;
			//	outfile<<"剪枝效果 = "<<pruning<<endl;
			//	outfile<<"hp resued = "<<hp->used<<endl;
			//	outfile<<"剪枝个数 = "<<(int)(samplesize * pruning)<<endl;
			//	ticksCnt += ( end_tick - begin_tick );
			//	page_faults += ( c->page_faults );
			//	prunings += pruning;
			//	penaltys += penalty;
			//	outfile<<"--------------------\n";
			//}
			//ticksCnt  =  ticksCnt / ( float) 1000 / ( float ) qryCnt;
			//page_faults = page_faults / ( float ) 100 / ( float ) qryCnt;
			//penaltys = penaltys / ( float ) qryCnt;
			//prunings = prunings / ( float ) qryCnt;
			//outfile<<"==========KWQ平均结果==========\n";
			//outfile<<"IO代价 = "<< page_faults<<endl;
			//outfile<<"CPU代价 =  "<< ticksCnt<<endl;		
			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
			//outfile<<"penalty = "<<penaltys <<endl;
			//outfile<<"剪枝效果 = "<<prunings<<endl;
			//outfile<<"hp resued = "<<hpUsed / qryCnt<<endl;
			//outfile<<"取样个数 = "<<samplesize<<endl;
			//outfile<<"剪枝个数 = "<<(int)( samplesize * prunings )<<endl;
			//outfile<<"========================================\n\n";
			//-----------------------------------end-----------------------------------------------------------------------		

			outfile.close();
			//---------------------------------------------
			/*qry_ranks_p = qry_ranks_p + *abW_p;
			Prs_p = Prs_p + 1;
			Ts_p = Ts_p + 1;
			ks_p = ks_p + 1;
			abW_p = abW_p + 1;
			pa_Cnt = pa_Cnt + 1;*/
			//---------------------------------------------
			if (rtree)
			{
				delete rtree;
				rtree = NULL;
			}
		}
		//-----------------------------------end-----------------------------------------------------------------------
		catch (exception e)
		{
			if (rtree)
			{
				delete rtree;
				rtree = NULL;
			}
		}
		delete c;

		/*do
		{
		cout << "Search another dataset?\n(dataset number, or '0' to quit, or '-1' to list all files)" << endl;
		cin >> dateset_no;
		if (dateset_no == -1)
		{
		for (i = 0; i < datesets.size(); i++)
		{
		cout << '[' << i+1 << ']' << '\t' << datesets[i] << endl;
		}
		}
		} while(dateset_no == -1);*/
	}
	delete[] fname;

	return 0;
}

//RTree* buildtree(char *_fname, int _dimension, Cache* c)
//{
//	RTree *_rtree;
//	string map_name(_fname);
//	string base_name;
//	string::size_type begin_idx = map_name.find_last_of('\\');
//	string::size_type end_idx = map_name.find_last_of('.');
//	if (end_idx == string::npos)
//		base_name = map_name.substr(begin_idx+1);
//	else
//		base_name = map_name.substr(begin_idx+1, end_idx-begin_idx-1);
//
//	cout << "[" << base_name << "]" << endl << endl;
//
//	rslt_fname = "c:\\";
//	rslt_fname = rslt_fname + base_name;
//
//	string index_name(".\\trees\\");
//	index_name += base_name + ".tree";
//
//	char *iname = new char[index_name.size()+1];
//	memcpy(iname, index_name.c_str(), index_name.size()+1);
//	ifstream index_file(index_name.c_str());
//	if (!index_file)
//	{
//		cout << "There is no rtree built already!" << endl;
//		cout << "Building rtree ..." << endl;
// 		//c = new Cache(0, 1024);
//		_rtree = new RTree(_fname, iname, 1024, c, _dimension);
//		cout << "OK, a rtree is built." << endl << endl;
//	}
//	else
//	{
//		cout << "There is a rtree already" << endl << endl;
//		//c = new Cache(0, 1024);
//		_rtree = new RTree(iname, c);
//	}	
//
//	return _rtree;
//}

RTree* buildtree(char *_fname, int _dimension, Cache* c)//有树就出来，没树就建树
{
	RTree *_rtree;
	string map_name(_fname);
	string base_name;
	string::size_type begin_idx = map_name.find_last_of('\\');
	string::size_type end_idx = map_name.find_last_of('.');
	if (end_idx == string::npos)
		base_name = map_name.substr(begin_idx + 1);
	else
		base_name = map_name.substr(begin_idx + 1, end_idx - begin_idx - 1);

	cout << "[" << base_name << "]" << endl << endl;

	rslt_fname = "c:\\";
	rslt_fname = rslt_fname + base_name;

	string index_name(".\\trees\\");
	index_name += base_name + ".tree";//建树的名称放在了index_name底下

	char *iname = new char[index_name.size() + 1];
	memcpy(iname, index_name.c_str(), index_name.size() + 1);//用index_name给iname赋值
	ifstream index_file(index_name.c_str());
	if (!index_file)
	{
		cout << "There is no rtree built already!" << endl;
		cout << "Building rtree ..." << endl;
		//c = new Cache(0, 1024);
		_rtree = new RTree(_fname, iname, 1024, c, _dimension);
		cout << "OK, a rtree is built." << endl << endl;
	}
	else
	{
		cout << "There is a rtree already" << endl << endl;
		//c = new Cache(0, 1024);
		_rtree = new RTree(iname, c);
	}
	index_file.close();
	return _rtree;
}

/*
* Description: 把矢量数据写入 filePath 的文件中
* Parameters:
*  set: 要写入的集合
*  dimension: 矢量的维度
*  size: 矢量个数
*  filePath: 文件的路径
*/
template <class T>
void writeSet(T* set, int dimension, int size, char* filePath)
{
	ofstream writeFile(filePath);
	for (int i = 0; i < size; i++)
	{
		writeFile << ' ' << setw(4) << right << i;
		for (int j = 0; j < dimension; j++)
			writeFile << ' ' << setw(12) << set[i * dimension + j];
		writeFile << endl;
	}
	writeFile.close();
}

template <class T>
void readSet(T* set, char* filePath, int dimension)
{
	ifstream readFile(filePath);
	int num;
	readFile >> num;
	while (num >= 0)
	{
		for (int i = 0; i < dimension; i++)
			readFile >> set[num* dimension + i];

		readFile >> num;
		if (readFile.eof())
			break;
	}
	readFile.close();
}

int readSet(float* queryPoints, float* weightSet, char* filePath, int dimension, int missCnt, int * missIdSet)//读取querypoints的一个文件，把缺失向量放在weightSet里，把查询点放在queryPoints里
{

	ifstream infile(filePath);
	char str[30];
	int data;
	int cnt = 0;
	int rank;
	int tablecnt;
	while (infile.peek() != EOF)
	{

		infile >> str;
		if (infile.peek() != EOF)
		{
			for (int i = 0; i < dimension; i++)
				infile >> queryPoints[cnt * dimension + i];
			infile >> tablecnt;
			for (int i = 0; i < tablecnt; i++)
				;
			
			for (int j = 0; j < missCnt; j++)
			{
				infile >> str;
				for (int i = 0; i < dimension; i++)
				{
					infile >> weightSet[cnt * dimension * missCnt + i + j * dimension];
				}
				for (int i = 0; i < 1; i++)
				{
					infile >> missIdSet[cnt*missCnt + i + j];
					cout <<	"缺失值ID："<<missIdSet[cnt*missCnt + i + j]<<endl;
				}
				infile.getline(str, 100, '\n');
			}
			cnt++;
		}
	}
	return cnt;
}
