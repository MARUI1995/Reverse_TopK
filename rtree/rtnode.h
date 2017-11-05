#ifndef __RTNODE
#define __RTNODE
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class SortedLinList;
class Entry;
class RTree;
class Heap;
//------------------------------------------------------------
class RTNode
{
public:
//--===on disk===--
	char level; 
	int block;
	int num_entries;
	Entry *entries;
//--===others===--
	bool dirty;
	int capacity;
    int dimension;
	RTree *my_tree;  
//--===functions===--
	RTNode(RTree *rt);
    RTNode(RTree *rt, int _block);
    ~RTNode();

    int choose_subtree(float *brm);
	R_DELETE delete_entry(Entry *e); 
	void enter(Entry *de);
	bool FindLeaf(Entry *e);
	float *get_mbr();
	int get_num_of_data();
	R_OVERFLOW insert(Entry *d, RTNode **sn);
	bool is_data_node() { return (level==0) ? TRUE : FALSE ;};
	void NNSearch(float *QueryPoint, SortedLinList *res,
				      float *nearest_distanz);
	void print();
	void rangeQuery(float * mbr, SortedLinList * res, long type);
	void rangeQuery(float *mbr, SortedLinList *res);
    void read_from_buffer(char *buffer);
	int split(float **mbr, int **distribution);
	void split(RTNode *sn);
	void write_to_buffer(char *buffer); 
	//--==added for ranked query==--
	void rank_qry_inquiry(float *_weight, float _qscore, int *_rslt);

	/*
	//--==added for skyline query==--
	bool skyline_mbr_inside_tree(float *mbr, int *dc); // by Greg
	int RTNode::countQuery(float *mbr); // by Greg
	void RTNode::boundedNNSearch_dynamic(float *QueryPoint, SortedLinList *res, float *bound, float *nearest_distanz, bool *func);
	//functions added by yufei tao
	void traverse(float *_skyarray, float *_skydomcnt, int _skyarrayused);
	bool main_memory_dom_check(int _k, float *_pt, int &_domcnt); 
	*/
};

#endif