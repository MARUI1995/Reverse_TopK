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
	int m_tableCnt; //���ӱ�ĸ�����Ϊʵ���������
	int m_attriCnt; //���Ը�����Ϊ��������
	
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
	// ִ�в鿴��ѯ���Ƿ��Ƿ�Top-k���
	//---------------------------------------------------------------------	
	float* missWeightSet = new float[missWeightCnt * dimension];//�Ե�ǰ���ݼ��ĵ�����ѯ���һ�鶪ʧȨ��������������һ���������Ƕ�����������ݼ������ÿһ����һ��Ȩ��
	float *_qry_pt = new float[dimension];//�Ե�ǰ���ݼ���һ����ѯ��
	qPoint *qpoint;//һ����ѯ��
	int * missId;//һ����ѯ���ȱʧ��������ID
	//float *parameter;
	clock_t begin_tick, end_tick;
	for (int cnt = 0; cnt < qryCnt; cnt++) //���ղ�ѯ���������ѭ��
	{
		//parameter = new float[100];
		//outfile << "��\t" << cnt << "\t����ѯ��:\n";
		missId = new int[missWeightCnt];
		memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��  _qry_pt
		memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //��ֵ��ȱʧȨ��
		//memcpy(parameter, &missIdSet[cnt*1*missWeightCnt], sizeof(float) * 1 * missWeightCnt);//��ֵ��Ȩ������
		memcpy(missId, &missIdSet[cnt * 1 * missWeightCnt], sizeof(int) * 1 * missWeightCnt);//��ֵ��Ȩ��ID
		c->page_faults = 0;//cache c
		qpoint = new qPoint(_qry_pt, dimension, missWeightCnt);
		qpoint->m_missID = missId;
		//------------------------------begin-------------------------------------------------------------------
		// ��ֵ����ѯ��ȱʧȨ�ص�����
		//---------------------------------------------------------------------	
		/*begin_tick = clock();
		float *currentPa = new float[1];
		float *minPa = new float[1];
		memcpy(minPa, parameter, sizeof(float) * dimension);//��ʼ��
		for (int missweightcnt = 0; missweightcnt < missWeightCnt; missweightcnt++)
		{
			memcpy(currentPa, &parameter[missweightcnt *1], sizeof(float) * dimension);//��ֵ����ѯ��  _qry_pt
			for (int j = 0; j < 1; j++)
			{
				if (currentPa[j] < minPa[j])
					minPa[j] = currentPa[j];
			}
		}
		end_tick = clock();*/
		//memcpy(qpoint->m_minParameter, minPa, sizeof(float) * dimension);//��ʼ��
		//qpoint->m_spjTime = (end_tick - begin_tick) / (float)1000;


		cout << "��ֵ����ѯ��ʱ��";
		for (int i = 0; i<dimension; i++)
		{
			cout << _qry_pt[i] << " ";

		}
		cout << endl;
		cout << "������Ӧqpoint��";
		for (int i = 0; i<dimension; i++)
		{
			cout << qpoint->m_qry_pt[i] << " ";

		}

		
		cout << endl;

		/*���㵱ǰ��ѯ���������Ͳ��ɱȵ�*/
		int ranking, rank;
		int bedomnum, incomnum;  //incomnum��incomparable��ĸ�����bedomnum��dominate��ѯ��q�ĵ�ĸ�����
		float *incomparable = new float[3000000];
		rtree->incomparable_points_reuse(hp, incomparable, incomnum, _qry_pt, bedomnum);   //��incomparable��ļ��ϣ���������incomparable�У������incomnum�� bedomnum��ֵ
		//outfile << "��������" << bedomnum - 1 << "��" << endl;//��Ϊ����������˲�ѯ���Լ�
		//outfile << "���ɱȵ���" << incomnum << "��" << endl;
		for (int j = 0; j < incomnum; j++)
		{


			for (int i = 0; i < dimension; i++)
			{
				//outfile << incomparable[j*dimension + i] << "��";

			}
			//outfile << endl;
		}

		float * currentWeight = new float[dimension];

		for (int misweightCnt = 0; misweightCnt < missWeightCnt; misweightCnt++)
		{
			memcpy(currentWeight, &missWeightSet[misweightCnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��  _qry_pt
			rtree->bbs_ranking_querypoint(incomparable, incomnum, _qry_pt, currentWeight, ranking);//�˴�����currentWeight�����ǵ���Ȩ��
			rank = ranking + bedomnum;
			//outfile << "��ǰrank=" << rank << "��" << "��ǰTopkֵ=" << k << endl;

			if (rank <= k)
			{
				//outfile << "==========�Ƿ�Top-k��������ø���==========" << "rank=" << rank << endl;
				memcpy(qpoint->m_reWeight + qpoint->m_reWeightcnt*dimension, currentWeight, sizeof(float) * dimension);
				qpoint->m_reWeightcnt++;

			}
			else
			{
				//outfile << "==========���Ƿ�Top-k�������Ҫ����==========" << "rank=" << rank << endl;
				memcpy(qpoint->m_notReWeight + qpoint->m_notReWeightcnt*dimension, currentWeight, sizeof(float) * dimension);
				qpoint->m_notReWeightcnt++;
			}
		}
		qvector.push_back(*qpoint);

		cout << "����һ��qpoint" << endl;
		cout << "�˲�ѯ����" << qpoint->m_missCnt << "��ȱʧ���" << "    " << "ȱʧ�����ID�ֱ�Ϊ";
		for (int i = 0; i < qpoint->m_missCnt; i++)
		{
			cout << qpoint->m_missID[i] << " " << endl;
		}
		cout << "��ǰqvector��СΪ" << qvector.size() << endl;
		cout << endl;
	}
}
void modifyKW(vector<qPoint> &qvector, char *fname, int qryCnt, int missWeightCnt, float* queryPoints, float* weightSet, Cache *c, RTree *rtree, int dimension, Heap *hp, long samplesize, float  quality_of_answer, float probability_guarantee,int k, int &modifyk)
{
	float* missWeightSet = new float[missWeightCnt * dimension];//�Ե�ǰ���ݼ��ĵ�����ѯ���һ�鶪ʧȨ��������������һ���������Ƕ�����������ݼ������ÿһ����һ��Ȩ��
	float *_qry_pt = new float[dimension];//�Ե�ǰ���ݼ���һ����ѯ��
	float *modifyweight = new float[dimension * missWeightCnt];//�޸ĺ��Ȩ��������ÿһ����һ��Ȩ��
	outfile.open(setFileName, ios::out | ios::app);
	cout << "�����޸�KW" << endl;
	outfile << "========== �޸� k w �Ľ����㷨 ==========\n";
	outfile << "����\t" << qryCnt << "\t����ѯ��\n";
	outfile << "���ݼ���" << fname << endl;
	outfile << endl;
	clock_t begin_tick, end_tick;
	float pruning = 0.0, prunings = 0.0;
	float penalty = 0.0, penaltys = 0.0;
	float ticksCnt = 0, page_faults = 0;
	int totalMissCnt = missWeightCnt;
	//float * minPa;
	//float spjTime = 0;
	//if(qvector.size()==qryCnt)
	//cout<<"qvector��С���ѯ����Ŀ��ͬ"<<endl;

	for (int cnt = 0; cnt < qryCnt; cnt++) //���ղ�ѯ���������ѭ��
	{
		outfile << "��\t" << cnt << "\t����ѯ��:\n";
		//memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��  _qry_pt
		//memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //��ֵ��ȱʧȨ�أ��˴�����ֵ��Ȩ�ر����ǷǷ�Top-k�����
		_qry_pt = qvector[cnt].m_qry_pt;
		missWeightSet = qvector[cnt].m_notReWeight;
		missWeightCnt = qvector[cnt].m_notReWeightcnt;
		//minPa = qvector[cnt].m_minParameter;
		//spjTime = qvector[cnt].m_spjTime;
		c->page_faults = 0;//cache c
		begin_tick = clock();
		rtree->modifyWandK_approximate_reuse(hp, modifyweight, modifyk, penalty, _qry_pt, missWeightSet, missWeightCnt, k, quality_of_answer, probability_guarantee, pruning);//�޸�kw
		end_tick = clock();

		//�����д���ļ�
		outfile << "---------- ��������----------\n";
		outfile << "��ѯ��q��( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << _qry_pt[i] << " , ";
		outfile << _qry_pt[dimension - 1] << " )" << endl;   //�����ǰ��ѯ��ĸ���ά�ȵ�ֵ

		outfile << "��top-k��ѯ��kֵ�� " << k << endl;  //ԭkֵ
		outfile << "ά��d = " << dimension << endl;
		outfile << "��ʧ��Ȩ���������ϴ�С=" << totalMissCnt << endl;
		outfile << "��ʧ�ķǷ�Top-kȨ�ؼ��ϴ�С|W| = " << missWeightCnt << ",��Ҫ����" << endl; //ȱʧȨ�ظ���
		outfile << "��Ҫ��������������:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			//outfile << "��why-not�����µ�rank = " << rank << " : ( ";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << missWeightSet[i * dimension + j] << " , ";
			outfile << missWeightSet[j] << " )" << endl;
		}

		outfile << "Pr = " << probability_guarantee << endl;
		outfile << "T = " << quality_of_answer << endl;

		outfile << "\n----------ִ�н��----------" << endl;
		outfile << "��������������:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			outfile << i << ". (";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << modifyweight[i*dimension + j] << " , ";
			outfile << modifyweight[j] << " )" << endl;
		}
		outfile << "k = " << modifyk << endl;
		//outfile << "IO���� = " << c->page_faults / (float)100 << endl;
		//outfile << "CPU���� =  " << (end_tick - begin_tick) / (float)1000 << endl;

		//ɾ����SPJԭʼ�޸�
		/*outfile << "������spj�������ԣ�" << endl;
		for (int i = 0; i < parameterCnt; i++)
		{
			outfile << minPa[i] << "	";
		}
		*/
		outfile << endl;
		outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
		outfile << "penalty = " << penalty << endl;
		//outfile << "��֦Ч�� = " << pruning << endl;
		//outfile << "��֦���� = " << (int)(samplesize * pruning) << endl;
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
	outfile << "==========KWƽ�����==========\n";
	//outfile << "IO���� = " << page_faults << endl;
	//outfile << "CPU���� =  " << ticksCnt << endl;
	outfile << "TotalTime = " << page_faults + ticksCnt << endl;
	outfile << "penalty = " << penaltys << endl;
	//outfile << "��֦Ч�� = " << prunings << endl;
	outfile << "ȡ������ = " << samplesize << endl;
	//outfile << "��֦���� = " << (int)(samplesize * prunings) << endl;
	outfile << "========================================\n\n";
	outfile.flush();
}

void modifyQ(queue< float* > &qmodifiedQueue, vector<qPoint> &qvector, int dimension, Cache *c, RTree *rtree, char *fname, int qryCnt, int missWeightCnt, float* queryPoints, float* weightSet, Heap *hp, int k)
{
	// outfile.open(setFileName, ios::out | ios::app);
	cout << "�����޸�q" << endl;
	float* missWeightSet = new float[missWeightCnt * dimension];//�Ե�ǰ���ݼ��ĵ�����ѯ���һ�鶪ʧȨ��������������һ���������Ƕ�����������ݼ������ÿһ����һ��Ȩ��
	float *_qry_pt = new float[dimension];//�Ե�ǰ���ݼ���һ����ѯ��
	float* qmodified = new float[dimension];//�޸ĺ��q����qmodified
	float qabsolute;
	float ticksCnt = 0.0, page_faults = 0.0;
	float penalty = 0.0, penaltys = 0.0;
	clock_t begin_tick, end_tick;
	//float * minPa;
	//float spjTime = 0;
	//memset(setFileName, '\0', 100);
	if (outfile.is_open())
		cout << "open";
	outfile << "========== �޸� q �ľ�ȷ�㷨 ==========\n---------- ��������----------\n";
	outfile << "����\t" << qryCnt << "\t����ѯ��\n";
	outfile << "���ݼ���" << fname << endl;

	for (int cnt = 0; cnt < qryCnt; cnt++)
	{
		outfile << "��\t" << cnt << "\t����ѯ��:\n";
		//memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��
		//memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt);//��ֵ��ȱʧ����
		_qry_pt = qvector[cnt].m_qry_pt;
		missWeightSet = qvector[cnt].m_notReWeight;
		missWeightCnt = qvector[cnt].m_notReWeightcnt;
		//minPa = qvector[cnt].m_minParameter;
		//spjTime = qvector[cnt].m_spjTime;
		c->page_faults = 0;

		if (missWeightCnt>0)
		{
			begin_tick = clock();
			rtree->modifyQ_accurate(hp, _qry_pt, missWeightSet, missWeightCnt, k, qmodified);//��ȷ�޸�Q
			end_tick = clock();
		}
		else
		{
			begin_tick = clock();
			memcpy(qmodified, _qry_pt, sizeof(float) * dimension);
			end_tick = clock();
		}

		qmodifiedQueue.push(qmodified);

		outfile << "��ѯ��q��( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << _qry_pt[i] << " , ";
		outfile << _qry_pt[dimension - 1] << " )" << endl;

		outfile << "��top-k��ѯ��kֵ�� " << k << endl;
		outfile << "ά��d = " << dimension << endl;
		outfile << "��ʧ�ķǷ�Top-k�������ϴ�С|W| = " << missWeightCnt << endl;
		outfile << "��ʧ����������:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			//outfile << "��why-not�����µ�rank = " << rank << " : ( ";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << missWeightSet[i * dimension + j] << " , ";
			outfile << missWeightSet[j] << " )" << endl;
		}

		outfile << "\n----------ִ�н��----------" << endl;
		outfile << "��ѯ���޸ĺ�: (";

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
		
		//ɾ����ԭʼSPJ����
		/*outfile << "������spj�������ԣ�" << endl;
		for (int i = 0; i < parameterCnt; i++)
		{
			outfile << minPa[i] << "	";
		}
		outfile << endl;
		*/

		//outfile << "IO���� = " << c->page_faults / (float)100 << endl;
		//outfile << "CPU���� =  " << (end_tick - begin_tick) / (float)1000 << endl;
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
	outfile << "==========AccurateQƽ�����==========\n";
	//outfile<<"IO���� = "<< page_faults<<endl;
	//outfile<<"CPU���� =  "<< ticksCnt<<endl;		
	outfile << "TotalTime = " << page_faults + ticksCnt << endl;
	outfile << "penalty = " << penaltys << endl;
	outfile << "========================================\n\n";
}

void modifyKWQ(queue< float* > &qmodifiedQueue, vector<qPoint> &qvector,  int dimension, Cache *c, RTree *rtree, char *fname, int qryCnt, int missWeightCnt, float* queryPoints, float* weightSet, Heap *hp, int k, long samplesize, float *modifyweight, int modifyk, float  quality_of_answer, float probability_guarantee)
{
	float* missWeightSet = new float[missWeightCnt * dimension];//�Ե�ǰ���ݼ��ĵ�����ѯ���һ�鶪ʧȨ��������������һ���������Ƕ�����������ݼ������ÿһ����һ��Ȩ��
	float *_qry_pt = new float[dimension];//�Ե�ǰ���ݼ���һ����ѯ��
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
	cout << "�����޸�KWQ" << endl;
	if (outfile.is_open())
		cout << "open";
	outfile << "========== �޸� q k w ���㷨 ==========\n---------- ��������----------\n";
	outfile << "����\t" << qryCnt << "\t����ѯ��\n";
	outfile << "���ݼ���" << fname << endl;
	int hpUsed = 0;

	for (int cnt = 0; cnt < qryCnt; cnt++)
	{
		outfile << "��\t" << cnt << "\t����ѯ��:\n";
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

		//�����д���ļ�
		outfile << "��ѯ��q��( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << _qry_pt[i] << " , ";
		outfile << _qry_pt[dimension - 1] << " )" << endl;
		hpUsed += hp->used;
		outfile << "�޸�Q�����õ��Ĳ�ѯ��: ( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << qmodified[i] << " , ";
		outfile << qmodified[dimension - 1] << " )" << endl;

		outfile << "��top-k��ѯ��kֵ�� " << k << endl;
		outfile << "ά��d = " << dimension << endl;
		outfile << "��ʧ��Ȩ���������ϴ�С=" << totalMissCnt << endl;
		outfile << "��ʧ�ķǷ�Top-kȨ�ؼ��ϴ�С|W| = " << missWeightCnt << ",��Ҫ����" << endl; //ȱʧȨ�ظ���
		outfile << "��Ҫ��������������:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			//outfile<<"��why-not�����µ�rank = "<<*(qry_ranks_p - i )<<" : ( ";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << missWeightSet[i*dimension + j] << " , ";
			outfile << missWeightSet[j] << " )" << endl;
		}

		outfile << "T = " << quality_of_answer << endl;
		outfile << "Pr = " << probability_guarantee << endl;

		outfile << "\n----------ִ�н��----------" << endl;
		outfile << "��������������:" << endl;
		for (int i = 0; i < missWeightCnt; i++)
		{
			outfile << i << ". (";
			int j;
			for (j = 0; j < dimension - 1; j++)
				outfile << modifyweight[i*dimension + j] << " , ";
			outfile << modifyweight[j] << " )" << endl;
		}
		outfile << "k = " << modifyk << endl;
		outfile << "�޸Ĳ�ѯ��q��( ";
		for (int i = 0; i < dimension - 1; i++)
			outfile << qmodified_kwq[i] << " , ";
		outfile << qmodified_kwq[dimension - 1] << " )" << endl;
		
		//ɾ��ԭʼSPJ�޸�
		/*outfile << "������spj�������ԣ�" << endl;
		for (int i = 0; i < parameterCnt; i++)
		{
			outfile << minPa[i] << "	";
		}
		outfile << endl;
		*/

		//outfile<<"IO���� = "<< c->page_faults / ( float )100<<endl;
		//outfile<<"CPU���� =  "<< ( end_tick - begin_tick ) / ( float )1000<<endl;		
		outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
		outfile << "penalty = " << penalty << endl;
		//outfile<<"��֦Ч�� = "<<pruning<<endl;
		//outfile<<"hp resued = "<<hp->used<<endl;
		//outfile<<"��֦���� = "<<(int)(samplesize * pruning)<<endl;
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
	outfile << "==========KWQƽ�����==========\n";
	//outfile<<"IO���� = "<< page_faults<<endl;
	//outfile<<"CPU���� =  "<< ticksCnt<<endl;		
	outfile << "TotalTime = " << page_faults + ticksCnt << endl;
	outfile << "penalty = " << penaltys << endl;
	//outfile<<"��֦Ч�� = "<<prunings<<endl;
	//outfile<<"hp resued = "<<hpUsed / qryCnt<<endl;
	outfile << "ȡ������ = " << samplesize << endl;
	//outfile<<"��֦���� = "<<(int)( samplesize * prunings )<<endl;
	outfile << "========================================\n\n";
}

// ������ɲ�ѯ�㣬��ΧΪmin,max,��д���ļ�, min = 0 
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
	//��Ȩ���ļ�������datasets�С�
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
//		c = new Cache(0, 1024);//�����ڴ�
//		rtree = buildtree(fname, dimension, c);//�Ե�ǰ���ݼ���RTree������fname·���µ��ļ�
//		cout << rtree->dimension << endl;
//		clock_t begin_tick, end_tick;
//		double duration = 0.0;//�㷨����ʱ��
//		Heap *hp = new Heap();
//		hp->init(dimension);//�ѳ�ʼ��
//
//		//---------------------begin--------------------------------------------
//		// Ѱ�Ҳ�ѯ��
//		//-----------------------------------------------------------------
//		float* queryPoints = new float[dimension * 1000];//��ѯ�㣬ÿһ����һ����ѯ�㣬dimensionά
//		float* weightSet = new float[dimension * 1000];//ȱʧȨ��������ÿһ����һ��Ȩ��������dimensionά
//		int qryCnt = 1; // ���������Ĳ�ѯ����� 
//		int missWeightCnt = 1;// ��ʧ����������
//		char setFileName[100];
//
//		memset(setFileName, '\0', 100);
//		sprintf(setFileName, ".\\querypoints\\4_1_Test.txt");//����ѯ���ļ�����ֵ����ѯ������_ȱʧ������_��ǰʹ�����ݼ�����
//
//		qryCnt = readSet(queryPoints, weightSet, setFileName, dimension, missWeightCnt);//��ȡquerypoints��һ���ļ�����ȱʧ��������weightSet��Ѳ�ѯ�����queryPoints��
//		cout << "����" << qryCnt << "����ѯ��" << endl;
//
//		float* missWeightSet = new float[missWeightCnt * dimension];//�Ե�ǰ���ݼ��ĵ�����ѯ���һ�鶪ʧȨ��������������һ���������Ƕ�����������ݼ������ÿһ����һ��Ȩ��
//		float *_qry_pt = new float[dimension];//�Ե�ǰ���ݼ���һ����ѯ��
//
//		float *modifyweight = new float[dimension * missWeightCnt];//�޸ĺ��Ȩ��������ÿһ����һ��Ȩ��
//		float  quality_of_answer = 0.5, probability_guarantee = 0.5;//��������Ŀռ��СS
//		int  k = 3, modifyk=k;
//		float pruning = 0.0, prunings = 0.0;
//		float penalty = 0.0, penaltys = 0.0;
//		//ofstream outfile;
//		float ticksCnt = 0, page_faults = 0;
//		long samplesize = (long)(log(1 - probability_guarantee) / log(1 - quality_of_answer / 100)) + 1;//��������Ŀռ��СS
//		cout << "�����ռ��СΪ" << samplesize << endl;
//		memset(setFileName, '\0', 100);
//		sprintf(setFileName, ".\\results\\%s", &fname[10]);//���������ļ�����ֵ��setFilename
//		outfile.open(setFileName, ios::out | ios::app);//���û���ļ�����ô���ɿ��ļ���������ļ�����ô���ļ�β׷�ӡ�
//
//		//------------------------------begin-------------------------------------------------------------------
//		// ִ�в鿴��ѯ���Ƿ��Ƿ�Top-k���
//		//---------------------------------------------------------------------	
//		//for (int cnt = 0; cnt < qryCnt; cnt++) //���ղ�ѯ���������ѭ��
//		//{
//		//	outfile << "��\t" << cnt << "\t����ѯ��:\n";
//		//	memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��  _qry_pt
//		//	memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //��ֵ��ȱʧȨ��
//		//	c->page_faults = 0;//cache c
//		//	/*���㵱ǰ��ѯ���������Ͳ��ɱȵ�*/
//		//	int ranking, rank;
//		//	int bedomnum, incomnum;  //incomnum��incomparable��ĸ�����bedomnum��dominate��ѯ��q�ĵ�ĸ�����
//		//	float *incomparable = new float[3000000];
//		//	rtree->incomparable_points_reuse(hp, incomparable, incomnum, _qry_pt, bedomnum);   //��incomparable��ļ��ϣ���������incomparable�У������incomnum�� bedomnum��ֵ
//		//	outfile << "��������" << bedomnum - 1 << "��" << endl;//��Ϊ����������˲�ѯ���Լ�
//		//	outfile << "���ɱȵ���" << incomnum << "��" << endl;
//		//	for (int j = 0; j < incomnum; j++)
//		//	{
//		//		for (int i = 0; i < dimension; i++)
//		//		{
//		//			outfile << incomparable[j*dimension + i] << "��";
//		//		}
//		//		outfile << endl;
//		//	}
//		//	float * currentWeight = new float[dimension];
//		//	for (int misweightCnt = 0; misweightCnt < missWeightCnt; misweightCnt++)
//		//	{
//		//		memcpy(currentWeight, &missWeightSet[misweightCnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��  _qry_pt
//		//		rtree->bbs_ranking_querypoint(incomparable, incomnum, _qry_pt, currentWeight, ranking);//�˴�����currentWeight�����ǵ���Ȩ��
//		//		rank = ranking + bedomnum;
//		//		outfile << "��ǰrank=" << rank << "��" << "��ǰTopkֵ=" << k << endl;
//		//		if (rank <= k)
//		//		{
//		//			outfile << "==========�Ƿ�Top-k��������ø���==========" << "rank=" << rank << endl;
//		//		}
//		//		else outfile << "==========���Ƿ�Top-k�������Ҫ����==========" << "rank=" << rank << endl;
//		//	}
//		//}
//	
//		vector<qPoint> qvector;
//		cout<<"������vector"<<endl;
//		ifReverse(qvector, qryCnt, _qry_pt, queryPoints, missWeightSet, weightSet, missWeightCnt, dimension, c, rtree, hp, k);
//		cout<<"����˷�Top-k��ѯ"<<endl;
//		cout<<"��ǰqvector��СΪ"<<qvector.size()<<endl;
//		/*for(int i=0;i<qryCnt;i++)
//		{
//			int missCnt=qvector[i].m_missCnt;
//			cout<<"��ǰΪ��"<<i<<"����ѯ��:	";
//			for(int j=0;j<dimension;j++)
//				cout<<qvector[i].m_qry_pt[j]<<"	";
//			cout<<endl<<"����"<<missCnt<<"��ȱʧȨ��	"<<endl;
//			cout<<qvector[i].m_notReWeightcnt<<"���Ƿ�Top-k�����Ȩ��:"<<endl;
//			for(int j=0;j<qvector[i].m_notReWeightcnt;j++)
//			{
//				for(int z=0;z<dimension;z++)
//				{
//					cout<<qvector[i].m_notReWeight[j*dimension+z]<<",";
//				}
//				cout<<endl;
//			}
//			cout<<qvector[i].m_reWeightcnt<<"���Ƿ�Top-k�����Ȩ��:"<<endl;
//			for(int j=0;j<qvector[i].m_reWeightcnt;j++)
//			{
//				for(int z=0;z<dimension;z++)
//				{
//					cout<<qvector[i].m_reWeight[j*dimension+z]<<",";
//				}
//				cout<<endl;
//			}
//		}
//		cout<<"�鿴�Ƿ��Ƿ�Top-k���"<<endl;*/
//		
//		modifyKW(qvector,_qry_pt,fname, qryCnt, missWeightCnt, queryPoints, missWeightSet, weightSet, c, rtree, dimension, hp, modifyweight, samplesize, quality_of_answer, probability_guarantee, k, modifyk);
//		cout<<"kw�޸����"<<endl;
//		queue< float* > qmodifiedQueue;
//		modifyQ(qmodifiedQueue,qvector,_qry_pt,dimension,c,rtree,fname,qryCnt,missWeightCnt,queryPoints,missWeightSet,weightSet,hp,k);
//		cout<<"q�޸����"<<endl;
//		modifyKWQ(qmodifiedQueue,qvector,_qry_pt,dimension,c,rtree,fname,qryCnt,missWeightCnt,queryPoints,missWeightSet,weightSet,hp,k,samplesize,modifyweight,modifyk,quality_of_answer,probability_guarantee);
//		cout<<"kwq�޸����"<<endl;
//
//		//Sleep(100000);
//		//------------------------------begin-------------------------------------------------------------------
//		// ִ�н����޸�k w ���㷨
//		//---------------------------------------------------------------------	
//		//outfile << "========== �޸� k w �Ľ����㷨 ==========\n";
//		//outfile << "����\t" << qryCnt << "\t����ѯ��\n";
//		//outfile << "���ݼ���" << fname << endl;
//		//outfile << endl;
//		//for (int cnt = 0; cnt < qryCnt; cnt++) //���ղ�ѯ���������ѭ��
//		//{
//		//	outfile << "��\t" << cnt << "\t����ѯ��:\n";
//		//	memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��  _qry_pt
//		//	memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt); //��ֵ��ȱʧȨ��
//		//	c->page_faults = 0;//cache c			
//		//			begin_tick = clock();
//		//			rtree->modifyWandK_approximate_reuse(hp, modifyweight, modifyk, penalty, _qry_pt, missWeightSet, missWeightCnt, k, quality_of_answer, probability_guarantee, pruning);//�޸�kw
//		//			end_tick = clock();
//		//			//�����д���ļ�
//		//			outfile << "---------- ��������----------\n";
//		//			outfile << "��ѯ��q��( ";
//		//			for (int i = 0; i < dimension - 1; i++)
//		//				outfile << _qry_pt[i] << " , ";
//		//			outfile << _qry_pt[dimension - 1] << " )" << endl;   //�����ǰ��ѯ��ĸ���ά�ȵ�ֵ
//		//			outfile << "��top-k��ѯ��kֵ�� " << k << endl;  //ԭkֵ
//		//			outfile << "ά��d = " << dimension << endl;
//		//			outfile << "��ʧ���������ϴ�С|W| = " << missWeightCnt << endl; //ȱʧȨ�ظ���
//		//			outfile << "��ʧ����������:" << endl;
//		//			for (int i = 0; i < missWeightCnt; i++)
//		//			{
//		//				//outfile << "��why-not�����µ�rank = " << rank << " : ( ";
//		//				int j;
//		//				for (j = 0; j < dimension - 1; j++)
//		//					outfile << missWeightSet[i * dimension + j] << " , ";
//		//				outfile << missWeightSet[j] << " )" << endl;
//		//			}
//		//			outfile << "Pr = " << probability_guarantee << endl;
//		//			outfile << "T = " << quality_of_answer << endl;
//		//			outfile << "\n----------ִ�н��----------" << endl;
//		//			outfile << "��������:" << endl;
//		//			for (int i = 0; i < missWeightCnt; i++)
//		//			{
//		//				outfile << i << ". (";
//		//				int j;
//		//				for (j = 0; j < dimension - 1; j++)
//		//					outfile << modifyweight[i*dimension + j] << " , ";
//		//				outfile << modifyweight[j] << " )" << endl;
//		//			}
//		//			outfile << "k = " << modifyk << endl;
//		//			outfile << "IO���� = " << c->page_faults / (float)100 << endl;
//		//			outfile << "CPU���� =  " << (end_tick - begin_tick) / (float)1000 << endl;
//		//			outfile << "TotalTime = " << c->page_faults / (float)100 + (end_tick - begin_tick) / (float)1000 << endl;
//		//			outfile << "penalty = " << penalty << endl;
//		//			outfile << "��֦Ч�� = " << pruning << endl;
//		//			outfile << "��֦���� = " << (int)(samplesize * pruning) << endl;
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
//		//outfile << "==========KWƽ�����==========\n";
//		//outfile << "IO���� = " << page_faults << endl;
//		//outfile << "CPU���� =  " << ticksCnt << endl;
//		//outfile << "TotalTime = " << page_faults + ticksCnt << endl;
//		//outfile << "penalty = " << penaltys << endl;
//		//outfile << "��֦Ч�� = " << prunings << endl;
//		//outfile << "ȡ������ = " << samplesize << endl;
//		//outfile << "��֦���� = " << (int)(samplesize * prunings) << endl;
//		//outfile << "========================================\n\n";
//		//Sleep(100000);
//		//-----------------------------------end-----------------------------------------------------------------------
//
//		//------------------------------begin-------------------------------------------------------------------
//		// ִ�о�ȷ�޸�Q ���㷨
//		//---------------------------------------------------------------------
//		//float* qmodified = new float[dimension];//�޸ĺ��q����qmodified
//		////queue< float* > qmodifiedQueue;
//		//float qabsolute;
//		//ticksCnt = 0.0, page_faults = 0.0;
//		//penalty = 0.0, penaltys = 0.0;
//		//memset(setFileName, '\0', 100);
//		//if(outfile.is_open())
//		//	cout<<"open";
//		//outfile<<"========== �޸� q �ľ�ȷ�㷨 ==========\n---------- ��������----------\n";
//		//outfile<<"����\t"<<qryCnt<<"\t����ѯ��\n";
//		//outfile<<"���ݼ���"<<fname<<endl;
//		//for (int cnt = 0; cnt < qryCnt; cnt++)
//		//{
//		//	outfile << "��\t" << cnt << "\t����ѯ��:\n";
//		//	memcpy(_qry_pt, &queryPoints[cnt * dimension], sizeof(float) * dimension);//��ֵ����ѯ��
//		//	memcpy(missWeightSet, &weightSet[cnt * dimension * missWeightCnt], sizeof(float) * dimension * missWeightCnt);//��ֵ��ȱʧ����
//		//	c->page_faults = 0;			
//		//	begin_tick = clock();
//		//	rtree->modifyQ_accurate(hp, _qry_pt, missWeightSet, missWeightCnt, k, qmodified);//��ȷ�޸�Q
//		//	end_tick = clock();
//		//	qmodifiedQueue.push(qmodified);
//		//	outfile << "��ѯ��q��( ";
//		//	for (int i = 0; i < dimension - 1; i++)
//		//		outfile << _qry_pt[i] << " , ";
//		//	outfile << _qry_pt[dimension - 1] << " )" << endl;
//		//	outfile << "��top-k��ѯ��kֵ�� " << k << endl;
//		//	outfile << "ά��d = " << dimension << endl;
//		//	outfile << "��ʧ���������ϴ�С|W| = " << missWeightCnt << endl;
//		//	outfile << "��ʧ����������:" << endl;
//		//	for (int i = 0; i < missWeightCnt; i++)
//		//	{
//		//		//outfile << "��why-not�����µ�rank = " << rank << " : ( ";
//		//		int j;
//		//		for (j = 0; j < dimension - 1; j++)
//		//			outfile << missWeightSet[i * dimension + j] << " , ";
//		//		outfile << missWeightSet[j] << " )" << endl;
//		//	}
//		//	outfile << "\n----------ִ�н��----------" << endl;
//		//	outfile << "��ѯ���޸�: (";
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
//		//	outfile << "IO���� = " << c->page_faults / (float)100 << endl;
//		//	outfile << "CPU���� =  " << (end_tick - begin_tick) / (float)1000 << endl;
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
//		//outfile<<"==========AccurateQƽ�����==========\n";
//		//outfile<<"IO���� = "<< page_faults<<endl;
//		//outfile<<"CPU���� =  "<< ticksCnt<<endl;		
//		//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
//		//outfile<<"penalty = "<<penaltys <<endl;
//		//outfile<<"========================================\n\n";
//		//-----------------------------------end-----------------------------------------------------------------------
//
//
////------------------------------begin-------------------------------------------------------------------
//// ִ���޸�K W Q ���㷨
////---------------------------------------------------------------------
//			//float* qmodified_kwq = new float[dimension];
//			//ticksCnt = 0.0, page_faults = 0.0;
//			//pruning = 0.0, prunings = 0.0;
//			//penalty = 0.0, penaltys = 0.0;
//			//memset(setFileName, '\0', 100);
//			//if(outfile.is_open())
//			//	cout<<"open";
//			//outfile<<"========== �޸� q k w ���㷨 ==========\n---------- ��������----------\n";
//			//outfile<<"����\t"<<qryCnt<<"\t����ѯ��\n";
//			//outfile<<"���ݼ���"<<fname<<endl;
//			//int hpUsed = 0;
//			//for( int cnt = 0; cnt < qryCnt; cnt++ )
//			//{
//			//	outfile<<"��\t"<<cnt<<"\t����ѯ��:\n";
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
//			//	 //�����д���ļ�
//			//	outfile<<"��ѯ��q��( ";
//			//	for( int i = 0; i < dimension - 1; i++ )
//			//		outfile<<_qry_pt[i]<<" , ";
//			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;
//			//	hpUsed += hp->used;
//			//	outfile<<"�޸�Q�����õ��Ĳ�ѯ��: ( ";
//			//	for( int i = 0; i < dimension - 1; i++ )
//			//		outfile<<qmodified[i]<<" , ";
//			//	outfile<<qmodified[dimension - 1]<<" )"<<endl;
//			//	outfile<<"��top-k��ѯ��kֵ�� "<<k<<endl;
//			//	outfile<<"ά��d = "<<dimension<<endl;
//			//	outfile<<"��ʧ���������ϴ�С|W| = "<<missWeightCnt<<endl;
//			//	outfile<<"��ʧ����������:"<<endl;
//			//	for( int i = 0; i < missWeightCnt; i++ )
//			//	{
//			//		//outfile<<"��why-not�����µ�rank = "<<*(qry_ranks_p - i )<<" : ( ";
//			//		int j;
//			//		for( j = 0; j < dimension - 1; j++ )
//			//			outfile<<missWeightSet[i*dimension + j]<<" , ";
//			//		outfile<<missWeightSet[j]<<" )"<<endl;
//			//	}
//			//	outfile<<"T = "<<quality_of_answer<<endl;
//			//	outfile<<"Pr = "<<probability_guarantee<<endl;
//			//	outfile<<"\n----------ִ�н��----------"<<endl;
//			//	outfile<<"��������:"<<endl;
//			//	for( int i = 0; i < missWeightCnt; i++ )
//			//	{
//			//		outfile<<i<<". (";
//			//		int j;
//			//		for( j = 0; j < dimension - 1; j++ )
//			//			outfile<<modifyweight[i*dimension + j]<<" , ";
//			//		outfile<<modifyweight[j]<<" )"<<endl;
//			//	}
//			//	outfile<<"k = "<<modifyk<<endl;
//			//	outfile<<"�޸Ĳ�ѯ��q��( ";
//			//	for( int i = 0; i < dimension - 1; i++ )
//			//		outfile<<qmodified_kwq[i]<<" , ";
//			//	outfile<<qmodified_kwq[dimension - 1]<<" )"<<endl;
//			//	outfile<<"IO���� = "<< c->page_faults / ( float )100<<endl;
//			//	outfile<<"CPU���� =  "<< ( end_tick - begin_tick ) / ( float )1000<<endl;		
//			//	outfile<<"TotalTime = "<< c -> page_faults / ( float )100 + ( end_tick - begin_tick ) / ( float )1000<<endl;
//			//	outfile<<"penalty = "<<penalty<<endl;
//			//	outfile<<"��֦Ч�� = "<<pruning<<endl;
//			//	outfile<<"hp resued = "<<hp->used<<endl;
//			//	outfile<<"��֦���� = "<<(int)(samplesize * pruning)<<endl;
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
//			//outfile<<"==========KWQƽ�����==========\n";
//			//outfile<<"IO���� = "<< page_faults<<endl;
//			//outfile<<"CPU���� =  "<< ticksCnt<<endl;		
//			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
//			//outfile<<"penalty = "<<penaltys <<endl;
//			//outfile<<"��֦Ч�� = "<<prunings<<endl;
//			//outfile<<"hp resued = "<<hpUsed / qryCnt<<endl;
//			//outfile<<"ȡ������ = "<<samplesize<<endl;
//			//outfile<<"��֦���� = "<<(int)( samplesize * prunings )<<endl;
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
	int parameterCnt = 1;//Ȩ���������Ը���,Ĭ��Ϊ1
	
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
		3, 3, 3, 3, 3,             // |W| = 1, 2, 3, 4, 5             ÿ��q��rank ��ͬ
		/*	13, 13, 13, 13, 13,
		23, 23, 23, 23, 23,
		33, 33, 33, 33, 33,*/
		2, 4, 5,                  // ά��
		//12, 14, 15,          //5
		//22, 24, 25,            //5
		//32, 34, 35,                                 //1
		6, 7, 10, 11,			//���ݼ���С
		/*16, 17, 20, 21,
		26, 27, 30, 31,
		36, 37, 40, 41,		*/
		3, 3, 3, 					// q ��why not �µ� rank
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
		10, 10, 10, 10, 10,                 // |W| = 1, 2, 3, 4, 5             ÿ��q��rank ��ͬ
		/*10, 10, 10,	10, 10,
		10, 10, 10,	10, 10,
		10, 10, 10,	10, 10,*/

		//=======================================================

		//============for real data begin=======================
		10, 10, 10,             // ά�� 
		/*10, 10, 10,
		10, 10, 10,
		10, 10, 10,*/
		10, 10, 10, 10,             //���ݼ���С
		/*10, 10, 10, 10,
		10, 10, 10, 10,
		10, 10, 10, 10, */
		//============for real data end=======================
		10, 10, 10,          // q ��why not �µ� rank
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
	};// ��top-k ��k ֵ
	int* ks_p = ks;
	int qry_ranks[] = {
		101,                             // |W| = 1, 2, 3, 4, 5             ÿ��q��rank ��ͬ
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
		101, 101, 101,                     // ά��
		/*101, 101, 101,
		101, 101, 101,
		101, 101, 101,*/
		101, 101, 101, 101,          // ���ݼ���С
		/*101, 101, 101, 101,
		101, 101, 101, 101,
		101, 101, 101, 101, */
		//============for real data end=======================
		11, 501, 1001,               // q ��why not �µ� rank
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
	};// ��ѯ���ڶ�ʧ����top-k��ѯ�е�rank
	int* qry_ranks_p = qry_ranks;
	float Prs[] = {
		0.4, 0.4, 0.4, 0.4, 0.4,       // |W| = 1, 2, 3, 4, 5             ÿ��q��rank ��ͬ
		/*0.8, 0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8, 0.8,*/
		//=======================================================

		//============for real data begin=======================
		0.8, 0.8, 0.8,            // ά�� 
		/*0.8, 0.8, 0.8,
		0.8, 0.8, 0.8,
		0.8, 0.8, 0.8,*/
		0.8, 0.8, 0.8, 0.8,         // ���ݼ���С
		/*0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,*/
		//============for real data end=======================
		0.8, 0.8, 0.8,            // q ��why not �µ� rank
		/*0.8, 0.8, 0.8,
		0.8, 0.8, 0.8,  */
		0.8, 0.8, 0.8, 0.8,         // k
		/*0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,  */
		0.8, 0.8, 0.8, 0.8,                 //���� Ts
		/*0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,
		0.8, 0.8, 0.8, 0.8,*/
		0.5, 0.6, 0.7, 0.9,  // ����pr
		/*0.5, 0.6, 0.7, 0.9,
		0.5, 0.6, 0.7, 0.9,
		0.5, 0.6, 0.7, 0.9*/
	};
	float *Prs_p = Prs;
	float Ts[] = {
		0.6, 0.6, 0.6, 0.6, 0.6,         // |W| = 1, 2, 3, 4, 5             ÿ��q��rank ��ͬ
		/*0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2, */
		//=======================================================
		//============for real data begin=======================
		0.2, 0.2, 0.2,            // ά��
		/*0.2, 0.2, 0.2,
		0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, */
		0.2, 0.2, 0.2, 0.2,            // ���ݼ���С
		/*	0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,*/
		//============for real data begin=======================
		0.2, 0.2, 0.2,                // q ��why not �µ� rank
		/*0.2, 0.2, 0.2,
		0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, */
		0.2, 0.2, 0.2, 0.2,                 // k
		/*0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, */
		0.1, 0.15, 0.25, 0.3,   // ���� ts
		/*0.1, 0.15, 0.25, 0.3,
		0.1, 0.15, 0.25, 0.3,
		0.1, 0.15, 0.25, 0.3,   */
		0.2, 0.2, 0.2, 0.2,       // ���� pr
		/*0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2*/
	};
	float *Ts_p = Ts;
	int abW[] = {
		1, 2, 3, 4, 5,                                    // |W| ÿ�� q��rank ��ͬ
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
	int paCnt[] = {          //Ȩ�����Ը�����Ĭ��Ϊ1
		1, 1, 1, 1,1,1
	};
	int * pa_Cnt = paCnt;
	//queue<int> qryPosTmp;  // ���������Ĳ�ѯ���ڲ�ѯ�����������е���ʼλ��
	//queue< float* > weightSetStackTmp; // ��ʧ����������

	for (int dataset_pos = 0; dataset_pos < 1/*27*/; dataset_pos++)//ѭ������
	{
		missWeightCnt = 5;
		try{
			//================
			//dateset_no = 3;
			//==============================
			dateset_no = dateset_nos[dataset_pos]; //ÿһ��ѭ�����ļ���
			//====================================			
			//cin >> dateset_no;
			//while (dateset_no > 0 && dateset_no <= datesets.size())
			//{
			memset(fname, 0, MAX_FILE_NAME);
			memcpy(fname, datesets[dateset_no - 1].c_str(), datesets[dateset_no - 1].size() + 1);//��datesets[dateset_no-1].size()+1��datesets[dateset_no-1].c_str()�е��ַ���ӡ���Ƶ�fname������ǰ���ݼ�����
			//strcpy(fname, dateset[dateset_no-1]);

			if (datesets[dateset_no - 1].find("_2d_") != string::npos)//datesets[dateset_no-1]Ϊ��ǰ���ݼ��ļ������ж�ά��
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

			c = new Cache(0, 1024);//�����ڴ�
			rtree = buildtree(fname, dimension, c);//�Ե�ǰ���ݼ���RTree������fname·���µ��ļ�
			cout << rtree->dimension << endl;


			/*while (query_number != -1)
			{*/
			clock_t begin_tick, end_tick;
			double duration = 0.0;//�㷨����ʱ��

			//float *rslt = new float[2 * MAX_SKYLINE_RESULT];
			//int rsltcnt = 0;
			Heap *hp = new Heap();
			hp->init(dimension);//�ѳ�ʼ��
			//string which_method;//ʲô����

			//////////////////////////�����������ϣ���д���ļ�//////////////////////////////////////////////////////////////////////
			//for( int d = 2; d <= 5; d++ )
			//{
			//	float* weights = rtree->SampleWeights(d, 1000, true);
			//	char fname[100];
			//	memset(fname, '\0', 100);
			//	sprintf(fname,".\\SampleWeights_%dd.txt", d );
			//	writeSet(weights, d, 1000, fname); 
			//}


			//---------------------begin--------------------------------------------
			// Ѱ�Ҳ�ѯ��
			//-----------------------------------------------------------------
			
			float* queryPoints = new float[dimension * 1000];//��ѯ�㣬ÿһ����һ����ѯ�㣬dimensionά
			float* weightSet = new float[dimension * 1000];//ȱʧȨ��������ÿһ����һ��Ȩ��������dimensionά
			int* missIdSet = new int[parameterCnt * 1000];
			int qryCnt = 0; // ���������Ĳ�ѯ����� 
			//missWeightCnt = *abW_p;// ��ʧ����������

			//char setFileName[100];

			memset(setFileName, '\0', 100);
			sprintf(setFileName, ".\\querypoints\\%d_%d_%s", *qry_ranks_p, missWeightCnt, &fname[10]);//����ѯ���ļ�����ֵ����ѯ������_ȱʧ������_��ǰʹ�����ݼ�����

			qryCnt = readSet(queryPoints, weightSet, setFileName, dimension, missWeightCnt,missIdSet);//��ȡquerypoints��һ���ļ�����ȱʧ��������weightSet��Ѳ�ѯ�����queryPoints��
			float* missWeightSet = new float[missWeightCnt * dimension];//�Ե�ǰ���ݼ��ĵ�����ѯ���һ�鶪ʧȨ��������������һ���������Ƕ�����������ݼ������ÿһ����һ��Ȩ��
			float *_qry_pt = new float[dimension];//�Ե�ǰ���ݼ���һ����ѯ��
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
			//		//qry_ranks_p = qry_ranks_p + *abW_p;// ����������ݼ��ҵ��˺��ʵĲ�ѯ�㣬��ô��ָ��ָ����һ����ѯ��ֵ
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
			//// ���������û���ҵ���Ӧ��why -not ������������һ��
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


			float *modifyweight = new float[dimension * missWeightCnt];//�޸ĺ��Ȩ��������ÿһ����һ��Ȩ��
			int  modifyk;
			float pruning = 0.0, prunings = 0.0;
			float penalty = 0.0, penaltys = 0.0;
			ofstream outfile;
			float ticksCnt = 0, page_faults = 0;
			long samplesize = (long)(log(1 - probability_guarantee) / log(1 - quality_of_answer / 100)) + 1;//��������Ŀռ��СS
			memset(setFileName, '\0', 100);
			sprintf(setFileName, ".\\results\\%d_%d_%s",queryRank,missWeightCnt, &fname[10]);//���������ļ�����ֵ��setFilename
			outfile.open(setFileName, ios::out | ios::app);//���û���ļ�����ô���ɿ��ļ���������ļ�����ô���ļ�β׷�ӡ�

			vector<qPoint> qvector;
			cout << "������vector" << endl;
			ifReverse(qvector, qryCnt, queryPoints, weightSet, missWeightCnt, dimension, c, rtree, hp, k, missIdSet);
			cout << "����˷�Top-k��ѯ" << endl;
			cout << "��ǰqvector��СΪ" << qvector.size() << endl;
			
			modifyKW(qvector, fname, qryCnt, missWeightCnt, queryPoints, weightSet, c, rtree, dimension, hp, samplesize, quality_of_answer, probability_guarantee, k, modifyk);
			cout << "kw�޸����" << endl;
			queue< float* > qmodifiedQueue;
			modifyQ(qmodifiedQueue, qvector, dimension, c, rtree, fname, qryCnt, missWeightCnt, queryPoints, weightSet, hp, k);
			cout << "q�޸����" << endl;
			modifyKWQ(qmodifiedQueue, qvector, dimension, c, rtree, fname, qryCnt, missWeightCnt, queryPoints, weightSet, hp, k, samplesize, modifyweight, modifyk, quality_of_answer, probability_guarantee);
			cout << "kwq�޸����" << endl;


			//------------------------------begin-------------------------------------------------------------------
			// ִ�н����޸�k w ���㷨
			//---------------------------------------------------------------------	
			//outfile<<"========== �޸� k w �Ľ����㷨 ==========\n";
			//outfile<<"����\t"<<qryCnt<<"\t����ѯ��\n";
			//outfile<<"���ݼ���"<<fname<<endl;
			//for( int cnt = 0; cnt < qryCnt; cnt++ ) //���ղ�ѯ���������ѭ��
			//{
			//	outfile<<"��\t"<<cnt<<"\t����ѯ��:\n";
			//	memcpy( _qry_pt,  &queryPoints[cnt * dimension], sizeof( float ) * dimension );//��ֵ����ѯ��  _qry_pt
			//	memcpy( missWeightSet,  &weightSet[cnt * dimension * missWeightCnt], sizeof( float ) * dimension * missWeightCnt); //��ֵ��ȱʧȨ��
			//	c->page_faults = 0;//cache c
			//	begin_tick = clock();
			//	rtree->modifyWandK_approximate_reuse(hp, modifyweight, modifyk, penalty, _qry_pt,  missWeightSet, missWeightCnt, k,  quality_of_answer, probability_guarantee, pruning);//�޸�kw
			//	end_tick = clock();
			//	//�����д���ļ�
			//	outfile<<"---------- ��������----------\n";
			//	outfile<<"��ѯ��q��( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<_qry_pt[i]<<" , ";
			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;   //�����ǰ��ѯ��ĸ���ά�ȵ�ֵ
			//	outfile<<"��top-k��ѯ��kֵ�� "<<k<<endl;  //ԭkֵ
			//	outfile<<"ά��d = "<<dimension<<endl;
			//	outfile<<"��ʧ���������ϴ�С|W| = "<<missWeightCnt<<endl; //ȱʧȨ�ظ���
			//	outfile<<"��ʧ����������:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ ) 
			//	{
			//		outfile<<"��why-not�����µ�rank = "<<*(qry_ranks_p - i )<<" : ( ";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<missWeightSet[ i * dimension + j]<<" , ";
			//		outfile<<missWeightSet[j]<<" )"<<endl;
			//	}
			//	outfile<<"Pr = "<<probability_guarantee<<endl;
			//	outfile<<"T = "<<quality_of_answer<<endl;
			//	outfile<<"\n----------ִ�н��----------"<<endl;
			//	outfile<<"��������:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<i<<". (";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<modifyweight[i*dimension + j]<<" , ";
			//		outfile<<modifyweight[j]<<" )"<<endl;
			//	}
			//	outfile<<"k = "<<modifyk<<endl;
			//	outfile<<"IO���� = "<< c->page_faults / ( float ) 100<<endl;
			//	outfile<<"CPU���� =  "<<( end_tick - begin_tick ) / ( float )1000<<endl;		
			//	outfile<<"TotalTime = "<<c->page_faults /( float ) 100 +( end_tick - begin_tick ) / ( float )1000<<endl;
			//	outfile<<"penalty = "<<penalty<<endl;
			//	outfile<<"��֦Ч�� = "<<pruning<<endl;
			//	outfile<<"��֦���� = "<<(int)(samplesize * pruning)<<endl;
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
			//outfile<<"==========KWƽ�����==========\n";
			//outfile<<"IO���� = "<< page_faults<<endl;
			//outfile<<"CPU���� =  "<< ticksCnt<<endl;		
			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
			//outfile<<"penalty = "<<penaltys <<endl;
			//outfile<<"��֦Ч�� = "<<prunings<<endl;
			//outfile<<"ȡ������ = "<<samplesize<<endl;
			//outfile<<"��֦���� = "<<(int)( samplesize * prunings )<<endl;
			//outfile<<"========================================\n\n";
			//-----------------------------------end-----------------------------------------------------------------------

			//------------------------------begin-------------------------------------------------------------------
			// ִ�о�ȷ�޸�Q ���㷨
			//---------------------------------------------------------------------
			//float* qmodified = new float[dimension];
			//queue< float* > qmodifiedQueue;
			//float qabsolute;
			//ticksCnt = 0.0, page_faults = 0.0;
			//penalty = 0.0, penaltys = 0.0;
			//memset(setFileName, '\0', 100);
			//if(outfile.is_open())
			//	cout<<"open";
			//outfile<<"========== �޸� q �ľ�ȷ�㷨 ==========\n---------- ��������----------\n";
			//outfile<<"����\t"<<qryCnt<<"\t����ѯ��\n";
			//outfile<<"���ݼ���"<<fname<<endl;
			//for( int cnt = 0; cnt < qryCnt; cnt++ )
			//{
			//	outfile<<"��\t"<<cnt<<"\t����ѯ��:\n";
			//	memcpy( _qry_pt,  &queryPoints[cnt * dimension], sizeof( float ) * dimension );//��ֵ����ѯ��
			//	memcpy( missWeightSet,  &weightSet[cnt * dimension * missWeightCnt], sizeof( float ) * dimension * missWeightCnt);//��ֵ��ȱʧ����
			//	c->page_faults = 0.0;
			//	begin_tick = clock();
			//	rtree->modifyQ_accurate(hp, _qry_pt,  missWeightSet, missWeightCnt, k, qmodified);//��ȷ�޸�Q
			//	end_tick = clock();
			//	qmodifiedQueue.push( qmodified );			
			//	outfile<<"��ѯ��q��( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<_qry_pt[i]<<" , ";
			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;
			//	outfile<<"��top-k��ѯ��kֵ�� "<<k<<endl;
			//	outfile<<"ά��d = "<<dimension<<endl;
			//	outfile<<"��ʧ���������ϴ�С|W| = "<<missWeightCnt<<endl;
			//	outfile<<"��ʧ����������:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<"��why-not�����µ�rank = "<<* ( qry_ranks_p - i )<<" : ( ";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<missWeightSet[ i * dimension + j ]<<" , ";
			//		outfile<<missWeightSet[j]<<" )"<<endl;
			//	}
			//	outfile<<"\n----------ִ�н��----------"<<endl;
			//	outfile<<"��ѯ���޸�: (";
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
			//	outfile<<"IO���� = "<< c->page_faults / ( float )100<<endl;
			//	outfile<<"CPU���� =  "<<( end_tick - begin_tick ) / ( float ) 1000<<endl;		
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
			//outfile<<"==========AccurateQƽ�����==========\n";
			//outfile<<"IO���� = "<< page_faults<<endl;
			//outfile<<"CPU���� =  "<< ticksCnt<<endl;		
			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
			//outfile<<"penalty = "<<penaltys <<endl;
			//outfile<<"========================================\n\n";
			//-----------------------------------end-----------------------------------------------------------------------

			//------------------------------begin-------------------------------------------------------------------
			// ִ���޸�K W Q ���㷨
			//---------------------------------------------------------------------
			//float* qmodified_kwq = new float[dimension];
			//ticksCnt = 0.0, page_faults = 0.0;
			//pruning = 0.0, prunings = 0.0;
			//penalty = 0.0, penaltys = 0.0;
			//memset(setFileName, '\0', 100);
			//if(outfile.is_open())
			//	cout<<"open";
			//outfile<<"========== �޸� q k w ���㷨 ==========\n---------- ��������----------\n";
			//outfile<<"����\t"<<qryCnt<<"\t����ѯ��\n";
			//outfile<<"���ݼ���"<<fname<<endl;
			//int hpUsed = 0;
			//for( int cnt = 0; cnt < qryCnt; cnt++ )
			//{
			//	outfile<<"��\t"<<cnt<<"\t����ѯ��:\n";
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
			//	 //�����д���ļ�
			//	outfile<<"��ѯ��q��( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<_qry_pt[i]<<" , ";
			//	outfile<<_qry_pt[dimension - 1]<<" )"<<endl;
			//	hpUsed += hp->used;
			//	outfile<<"�޸�Q�����õ��Ĳ�ѯ��: ( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<qmodified[i]<<" , ";
			//	outfile<<qmodified[dimension - 1]<<" )"<<endl;
			//	outfile<<"��top-k��ѯ��kֵ�� "<<k<<endl;
			//	outfile<<"ά��d = "<<dimension<<endl;
			//	outfile<<"��ʧ���������ϴ�С|W| = "<<missWeightCnt<<endl;
			//	outfile<<"��ʧ����������:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<"��why-not�����µ�rank = "<<*(qry_ranks_p - i )<<" : ( ";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<missWeightSet[i*dimension + j]<<" , ";
			//		outfile<<missWeightSet[j]<<" )"<<endl;
			//	}
			//	outfile<<"T = "<<quality_of_answer<<endl;
			//	outfile<<"Pr = "<<probability_guarantee<<endl;
			//	outfile<<"\n----------ִ�н��----------"<<endl;
			//	outfile<<"��������:"<<endl;
			//	for( int i = 0; i < missWeightCnt; i++ )
			//	{
			//		outfile<<i<<". (";
			//		int j;
			//		for( j = 0; j < dimension - 1; j++ )
			//			outfile<<modifyweight[i*dimension + j]<<" , ";
			//		outfile<<modifyweight[j]<<" )"<<endl;
			//	}
			//	outfile<<"k = "<<modifyk<<endl;
			//	outfile<<"�޸Ĳ�ѯ��q��( ";
			//	for( int i = 0; i < dimension - 1; i++ )
			//		outfile<<qmodified_kwq[i]<<" , ";
			//	outfile<<qmodified_kwq[dimension - 1]<<" )"<<endl;
			//	outfile<<"IO���� = "<< c->page_faults / ( float )100<<endl;
			//	outfile<<"CPU���� =  "<< ( end_tick - begin_tick ) / ( float )1000<<endl;		
			//	outfile<<"TotalTime = "<< c -> page_faults / ( float )100 + ( end_tick - begin_tick ) / ( float )1000<<endl;
			//	outfile<<"penalty = "<<penalty<<endl;
			//	outfile<<"��֦Ч�� = "<<pruning<<endl;
			//	outfile<<"hp resued = "<<hp->used<<endl;
			//	outfile<<"��֦���� = "<<(int)(samplesize * pruning)<<endl;
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
			//outfile<<"==========KWQƽ�����==========\n";
			//outfile<<"IO���� = "<< page_faults<<endl;
			//outfile<<"CPU���� =  "<< ticksCnt<<endl;		
			//outfile<<"TotalTime = "<< page_faults + ticksCnt <<endl;
			//outfile<<"penalty = "<<penaltys <<endl;
			//outfile<<"��֦Ч�� = "<<prunings<<endl;
			//outfile<<"hp resued = "<<hpUsed / qryCnt<<endl;
			//outfile<<"ȡ������ = "<<samplesize<<endl;
			//outfile<<"��֦���� = "<<(int)( samplesize * prunings )<<endl;
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

RTree* buildtree(char *_fname, int _dimension, Cache* c)//�����ͳ�����û���ͽ���
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
	index_name += base_name + ".tree";//���������Ʒ�����index_name����

	char *iname = new char[index_name.size() + 1];
	memcpy(iname, index_name.c_str(), index_name.size() + 1);//��index_name��iname��ֵ
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
* Description: ��ʸ������д�� filePath ���ļ���
* Parameters:
*  set: Ҫд��ļ���
*  dimension: ʸ����ά��
*  size: ʸ������
*  filePath: �ļ���·��
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

int readSet(float* queryPoints, float* weightSet, char* filePath, int dimension, int missCnt, int * missIdSet)//��ȡquerypoints��һ���ļ�����ȱʧ��������weightSet��Ѳ�ѯ�����queryPoints��
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
					cout <<	"ȱʧֵID��"<<missIdSet[cnt*missCnt + i + j]<<endl;
				}
				infile.getline(str, 100, '\n');
			}
			cnt++;
		}
	}
	return cnt;
}
