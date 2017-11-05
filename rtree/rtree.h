/* rtree.h
   this file defines the class RTree*/

#ifndef __RTREE
#define __RTREE
//------------------------------------------------------------
#include "../func/gendef.h"
#include "../heap/heap.h"
#include <vector>
#include <list>
using namespace std;

#define MAX_SKYLINE_RESULT	102400
#define DEFAULT_VALUE		0x00000100

//------------------------------------------------------------
class LinList;
class SortedLinList;
class Cache;
class RTNode;
class Entry;
//------------------------------------------------------------

struct value_skyline_point
{
	int i;	 // index of skyline point in skyline list
	int num; // the number of points dominated by skyline i
};

class RTree : public Cacheable
{
public:
	 
//--===on disk===--
   // int total_insert;
	int dimension;                       
	int num_of_data;	                 
    int num_of_dnodes;	                 
    int num_of_inodes;	                 
	int root;                            
	bool root_is_data;                   
//--===others===--
	RTNode *root_ptr;
    bool *re_level;  
    LinList *re_data_cands; 
	LinList *deletelist;
//--===added for cost summary===--
	int na[10];
//--===functions===--
	RTree(char *fname, int _b_length, Cache* c, int _dimension);
    RTree(char *fname, Cache* c);
    RTree(char *inpname, char *fname, int _blength, Cache* c, int _dimension);
	//RTree(vector<CInfoPoint> points, char *trname, int _blength, Cache* c, int _dimension);
    ~RTree();

	bool delete_entry(Entry *d);
	bool FindLeaf(Entry *e);
    int get_num() { return num_of_data; }
	void insert(Entry *d);
	void load_root();  
	void NNQuery(float *QueryPoint, SortedLinList *res);
	void rangeQuery(float * mbr, SortedLinList * res, long type);
	void rangeQuery(float *mbr, SortedLinList *res);
	void read_header(char *buffer);      
	void write_header(char *buffer);  
//--===added for ranked queries===--
	float get_score_linear(Entry *_e, float *_weight, int _dimension);
	float get_score_linear(float *_mbr, float *_weight, int _dimension);
	float get_score_range (Entry *_e, float *_weight, int _dimension, float *_range);
	float get_score_linear(Entry *_e, float *_weight, int _dimension, 
							  float SCORE_FUNC(const float *, const float *, int));
	void rank_qry_constr_linear(float *_weight, float *_qmbr, int _k, Heap *_hp, int *_rslt);
	void rank_qry_inquiry(float *_weight, float _qscore, int *_rslt);
	void rank_qry_linear(float *_weight, int _k, Heap *_hp, int *_rslt);
	void rank_qry_monotone(float *_weight, int _k, Heap *_hp, int *_rslt, 
							float SCORE_FUNC(const float *, const float *, int));
//--===added for skyline queries===--
	bool sky_dom(int _dim, float *_rslt, int _rsltcnt, float *_pt);
	void bbs_subspace(Heap *_hp, float *_rslt, int &_rsltcnt, bool *_active_dim);

	//added by Junfeng Hu on Oct 3, 2007
	void bbs(Heap *_hp, float *_rslt, int &_rsltcnt);
	void bbs_constrained(Heap *_hp, float *_rslt, int &_rsltcnt, float *_bounces);

//--==added for reverse skyline queries===--
	bool bool_range_query(Heap *_hp, float *_mbr, float *_gpt);
	bool global_dom(int _dim, float *_skyline, int _skylinecnt, float *_mbr, float *_qry_pt);
	void bbrs(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt);

	void rssa(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt);
	bool DDR(float* _p, float *_q);
	bool DADR(float* _p, float *_q);

	void bbs_topk(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, float* weighting, int k);
	void bbs_topkth(Heap *_hp, float *_rslt, float* _qry_pt, float* weighting, int k);
	void incomparable_points(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int &bedominatnum);
	void bbs_ranking_querypoint(float *_rslt, int &_rsltcnt, float* _qry_pt, float* weighting, int &ranking_querypoint);
	void modifyQ_accurate(Heap *_hp,  float* _qry_pt, float* weightset,int weightnum, int k, float *&qpmodified );
	void modifyQ_approximate(Heap *_hp, float* _qry_pt, float* weightset,int &weightnum, int k, float *qpmodified);
	void point_nearest_q( float *_rslt, float* _qry_pt, float* weight, float *kthpoint);
	//void modifyWandK_accurate(Heap *_hp, float *modifyweight, int &modifyk, float* _qry_pt, float* weightset,int &weightnum, int k, float *qpmodified, float &penalty);
	void modifyWandK_approximate(Heap *_hp, float *modifyweight, int &modifyk, float &penalty, float* _qry_pt, float* weightset,int &weightnum, int k, float quality_of_answer, float probability_guarantee, float& pruning);
	void SampleWeights(float quality_of_answer, float probability_guarantee, float *sample);
	void modify_k_W_Q(Heap *_hp,  float* _qry_pt, float* _qry_pt1, float* weightset, int &weightnum, int k, float *qpmodified, float* modifyweight,int &modifyk, float quality_of_answer, float probability_guarantee, float &penalty, float& pruning);
	void modify_k_W_Q_reuse(Heap *_hp,  float* _qry_pt, float* _qry_pt1, float* weightset, int &weightnum, int k, float *qpmodified, float* modifyweight,int &modifyk, float quality_of_answer, float probability_guarantee, float &penalty, float& pruning);
	void SampleWeightsbyq(float quality_of_answer, float probability_guarantee, float *sample, float *quwey_lower, float *quwey_upper );
	void SampleWeights( float* sample, int dimension, int sampleSize, bool isWeight );
	void SampleWeights( int sampleSize, float *sample, float* incomPoints, int incomnum, float* _qry_pt);
    float* checkQueryPoint(Heap *_hp, float* _qry_pt, float* weightSet, int weightnum, int dimension, int notCount, int* ranks, char* dataSetName);
	// This function is to compute reverse skyline as its definition describes based on dynamic_bbs
	void reverse_skyline(Heap *_hp, float *_rslt, int &_rsltcnt, float *_qry_pt);
	void modifyWandK_approximate_reuse(Heap *_hp, float *modifyweight, int &modifyk, float &penalty, float* _qry_pt, float* weightset,int &weightnum, int k, float quality_of_answer, float probability_guarantee, float& pruning);
	void dynamic_bbs(Heap *_hp, float *_rslt, int &_rsltcnt, float *_qry_pt);
	void incomparable_points_reuse(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int &bedominatnum);
	void traverse_points(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int &bedominatnum);
//--==added for variations of reverse skyline queries==--
	// constrained reverse skyline queries
	void bbrs_constrained(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, float *_bounces);
	void rssa_constrained(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, float *_bounces);

	// enumerating reverse skyline queries
	void bbrs_enum(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k);
	void rssa_enum(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k);
	int enum_range_query(Heap *_hp, float *_mbr, float *_gpt);

	// ranked reverse skyline queries
	void bbrs_ranked(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k, float (*_mindist)(int _dim, float*_bounces, float* _qry_pt));
	void rssa_ranked(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k, float (*_mindist)(int _dim, float*_bounces, float* _qry_pt));

	// reverse skybands
	void bbrs_skyband(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k);
	void rssa_skyband(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k);
	bool global_dom_skyband(int _dim, float *_skyline, int _skylinecnt, float *_mbr, float *_qry_pt, int _k);
};	

// mindist function for bbs
float sky_mindist(int _dim, float *_bounces);

// mindist function for dynamic bbs
float dynamic_mindist(int _dim, float *_bounces, float *_qry_pt);

// domination function for dynamic bbs 
bool dynamic_dom(int _dim, float *_rslt, int _rsltcnt, float *_bounces, float *_qry_pt);

#endif // __RTREE
