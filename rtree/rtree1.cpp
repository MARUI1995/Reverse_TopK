///*rtree.cpp
//this file implements the RTree class*/
//
//// switch to turn on/off debugging information 
////#define FILE_DEBUG	// output to files.
//#define CMD_DEBUG // output to the screen.
//
//#define _CRTDBG_MAP_ALLOC 
//#include<stdlib.h> 
//#include<crtdbg.h> 
//
//#include <vector>
//#include <math.h>
//#include <float.h>
//#include <string>
//#include "rtree.h"
//#include "entry.h"
//#include "rtnode.h"
//#include "../blockfile/cache.h"
//#include "../blockfile/blk_file.h"
//#include "../linlist/linlist.h"
//#include "../heap/heap.h"
//#include "../func/matlabEngine.h"
//#include <algorithm>
//#include <fstream>
//#include <iostream>
//#include <iomanip>
//#include <queue>
//#include <time.h>
//using namespace std;
//
////#include <afxwin.h>
//
//extern string rslt_fname;
//
////------------------------------------------------------------
//
//RTree::RTree(char *fname, int _b_length, Cache *c, int _dimension)
//	//use this constructor to build a new tree
//{
//	file = new BlockFile(fname, _b_length);
//	cache = c;
//
//	re_data_cands = new LinList();
//	deletelist = new LinList();
//
//	dimension = _dimension;
//	root = 0;
//	root_ptr = NULL;
//	root_is_data = TRUE;
//	num_of_data = num_of_inodes = num_of_dnodes = 0;
//
//	root_ptr = new RTNode(this);
//	//note that when a tree is constructed, the root is automatically created
//	//though at this time there is no entry in the root yet.
//	num_of_dnodes++;
//	root_ptr -> level = 0;
//	root = root_ptr -> block;
//}
////------------------------------------------------------------
//RTree::RTree(char *fname, Cache *c)
//	//use this constructor to restore a tree from a file
//{
//	int j;
//
//	file = new BlockFile(fname, 0);
//	cache =c;
//
//	re_data_cands = new LinList();
//	deletelist = new LinList();
//
//	char *header = new char [file->get_blocklength()];
//	file -> read_header(header);
//	read_header(header);
//	delete [] header;
//
//	root_ptr = NULL;
//}
////------------------------------------------------------------
//RTree::RTree(char *inpname, char *fname, int _b_length, Cache *c, int _dimension)
//	// construct new R-tree from a specified input textfile with rectangles
//{
//	Entry *d;
//	FILE *fp;
//	file = new BlockFile(fname, _b_length);
//	cache = c;
//
//	re_data_cands = new LinList();
//	deletelist = new LinList();
//
//	/*
//	char *header = new char [file->get_blocklength()];
//	read_header(header);
//	delete [] header;
//	*/
//
//	dimension = _dimension;
//	root = 0;
//	root_ptr = NULL;
//	root_is_data = TRUE;
//	num_of_data = num_of_inodes = num_of_dnodes = 0;
//
//	root_ptr = new RTNode(this);
//	num_of_dnodes++;
//	root_ptr -> level = 0;
//	root = root_ptr->block;
//
//	int record_count = 0;
//
//	if((fp = fopen(inpname,"r")) == NULL)
//	{
//		delete this;
//		error("Cannot open R-Tree text file", TRUE);
//	}
//	else
//	{
//		int id=0;
//		long type;
//		float x0,y0,x1,y1;
//		while (!feof(fp))
//		{
//			record_count ++;
//			//id ++; //disable this variable if the id comes with the dataset
//
//			if (record_count%100 == 0)
//			{
//				for (int i = 0; i < 79; i ++)  //clear a line
//					printf("\b");
//				printf("inserting object %d", record_count);
//			}
//
//			fscanf(fp, "%d ", &id);
//			d = new Entry(dimension, NULL);
//			d -> son = id;
//			for (int i=0; i<dimension; i++)
//			{
//				//fscanf(fp, "%f %f ", &(d->bounces[2*i]), &(d->bounces[2*i+1]));
//				fscanf(fp, "%f ", &(d->bounces[2*i]));
//				d->bounces[2*i+1] = d->bounces[2*i];
//			}
//			d -> type = -1;
//
//			insert(d);
//			//d will be deleted in insert()
//		}
//	}
//
//	fclose(fp);
//
//	printf("\n");
//	// AfxMessageBox(_T("OK"));
//	printf("This R-Tree contains %d internal, %d data nodes and %d data\n",
//		num_of_inodes, num_of_dnodes, num_of_data);
//
//	delete root_ptr;
//	root_ptr = NULL;
//}
///*
//RTree::RTree(vector<CInfoPoint> points, char *trname, int _blength, Cache* c, int _dimension)
////construct a new tree from the vector "points"
//{
//FILE* f;
//Entry *d;
//file = new BlockFile(trname, _blength);
//cache =c;
//
//re_data_cands = new LinList();
//deletelist = new LinList();
//
//
//char *header = new char [file->get_blocklength()];
//read_header(header);
//delete [] header;
//
//
//dimension = _dimension;
//root = 0;
//root_ptr = NULL;
//root_is_data = TRUE;
//num_of_data = num_of_inodes = num_of_dnodes = 0;
//
//root_ptr = new RTNode(this);
//num_of_dnodes++;
//root_ptr -> level = 0;
//root = root_ptr->block;
//
//int record_count = 0;
//
//if(points.empty())
//{
//delete this;
//error("There is no data in vector!", TRUE);
//}
//else
//{	  
//int id=0;
////float x0,y0,x1,y1;
//vector<CInfoPoint>::iterator pos;
//for (pos = points.begin();pos != points.end(); ++pos)
//{
//record_count ++;
//
//if (record_count%100 == 0)
//{
//CWnd* processBar = ::AfxGetMainWnd()->GetWindow(GW_CHILD)->GetWindow(GW_CHILD);
//processBar->SendMessage(WM_SET_PROGRESS, record_count, points.size());
//}
//
//d = new Entry(dimension, NULL);
//d -> son = pos->id;
//d->type=pos->type;
//
//for(int i = 0; i < dimension; i++)
//{
//d->bounces[i] = pos->x;
//d->bounces[i+dimension] = pos->y;
//}
//
//insert(d);
////d will be deleted in insert()
//}
//}
//
////	fclose(fp);
////test
//for(int i = 0; i < 10; i++)
//na[i] = 0;
//
//delete root_ptr;
//root_ptr = NULL;
//}
//*/
//
////------------------------------------------------------------
//RTree::~RTree()
//{
//	int j;
//
//	char *header = new char[file -> get_blocklength()];
//	write_header(header);
//	file->set_header(header);
//	delete [] header;
//
//	if (root_ptr != NULL)
//	{
//		delete root_ptr;
//		root_ptr = NULL;
//	}
//
//	if (cache)
//		cache -> flush();
//
//	delete file;
//
//	delete re_data_cands;
//	delete deletelist;
//
//	//    printf("This R-Tree contains %d internal, %d data nodes and %d data\n",
//	//	   num_of_inodes, num_of_dnodes, num_of_data);
//}
////------------------------------------------------------------
//bool RTree::delete_entry(Entry *d)
//{
//	load_root();
//
//	R_DELETE del_ret;
//	del_ret=root_ptr->delete_entry(d);
//
//	if (del_ret == NOTFOUND) return false;
//	if (del_ret == ERASED) 
//		error("RTree::delete_entry--The root has been deleted\n",true);
//
//	if (root_ptr -> level > 0 && root_ptr -> num_entries == 1)
//		//there is only one entry in the root but the root
//		//is not leaf.  in this case, the child of the root is exhalted to root
//	{
//		root = root_ptr -> entries[0].son;
//		delete root_ptr;
//		root_ptr = NULL;
//		load_root();
//		num_of_inodes--;
//	}
//
//	num_of_data--;
//
//	//Now will reinsert the entries
//	while (deletelist -> get_num() > 0)
//	{
//		Linkable *e;
//		e = deletelist -> get_first();
//		Entry *new_e = new Entry(dimension, NULL);
//		new_e -> set_from_Linkable(e);
//		deletelist -> erase();
//		insert(new_e);
//		num_of_data--;
//	}
//
//	delete root_ptr;
//	root_ptr = NULL;
//
//	return true;
//}
////------------------------------------------------------------
//bool RTree::FindLeaf(Entry *e)
//{
//	load_root();
//	return root_ptr -> FindLeaf(e);
//}
////------------------------------------------------------------
//void RTree::insert(Entry* d)
//{
//	int i, j;
//	RTNode *sn;
//	RTNode *nroot_ptr;
//	int nroot;
//	Entry *de;
//	R_OVERFLOW split_root;
//	Entry *dc;
//	float *nmbr;
//
//	// load root into memory
//	load_root();
//
//	// no overflow occured until now
//	re_level = new bool[root_ptr -> level + 1];
//	for (i = 0; i <= root_ptr -> level; i++)
//		re_level[i] = FALSE;
//
//	// insert d into re_data_cands as the first entry to insert
//	// make a copy of d because it should be erased later
//	Linkable *new_link;
//	new_link = d -> gen_Linkable();
//	re_data_cands -> insert(new_link);
//
//	delete d;  //we follow the convention that the entry will be deleted when insertion finishes
//
//	j = -1;
//	while (re_data_cands -> get_num() > 0)
//	{
//		// first try to insert data, then directory entries
//		Linkable *d_cand;
//		d_cand = re_data_cands -> get_first();
//		if (d_cand != NULL)
//		{
//			// since "erase" deletes the data itself from the
//			// list, we should make a copy of the data before
//			// erasing it
//			dc = new Entry(dimension, NULL);
//			dc -> set_from_Linkable(d_cand);
//			re_data_cands -> erase();
//
//			// start recursive insert with root
//			split_root = root_ptr -> insert(dc, &sn);//error here! type in root_ptr doesn't changed
//		}
//		else
//			error("RTree::insert: inconsistent list re_data_cands", TRUE);
//
//		if (split_root == SPLIT)
//			// insert has lead to split --> new root-page with two sons (i.e. root and sn)
//		{
//			nroot_ptr = new RTNode(this);
//			nroot_ptr -> level = root_ptr -> level + 1;
//			num_of_inodes++;
//			nroot = nroot_ptr -> block;
//
//			de = new Entry(dimension, this);
//			nmbr = root_ptr -> get_mbr();
//			memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
//			delete [] nmbr;
//			de->son = root_ptr->block;
//			de->son_ptr = root_ptr;
//			nroot_ptr -> enter(de);
//
//			de = new Entry(dimension, this);
//			nmbr = sn -> get_mbr();
//			memcpy(de -> bounces, nmbr, 2*dimension*sizeof(float));
//			delete [] nmbr;
//			de -> son = sn -> block;
//			de -> son_ptr = sn;
//			nroot_ptr->enter(de);
//
//			root = nroot;
//			root_ptr = nroot_ptr;
//
//			root_is_data = FALSE;
//		}
//		j++;
//	}
//
//	num_of_data++;
//
//	delete [] re_level;
//
//	delete root_ptr;
//	root_ptr = NULL;
//}
////------------------------------------------------------------
//void RTree::load_root()
//{
//	if (root_ptr == NULL)
//		root_ptr = new RTNode(this, root);
//}
////------------------------------------------------------------
//void RTree::NNQuery(float *QueryPoint,
//	SortedLinList *res)
//{
//	float nearest_distanz;
//
//	// load root node into main memory
//	load_root();
//
//	nearest_distanz = MAXREAL;
//
//	root_ptr->NNSearch(QueryPoint,res,&nearest_distanz);
//
//	delete root_ptr;
//	root_ptr = NULL;
//}
////------------------------------------------------------------
//void RTree::rangeQuery(float *mbr, SortedLinList *res)
//{
//	load_root();
//
//	root_ptr -> rangeQuery(mbr,res);
//
//	delete root_ptr;
//	root_ptr = NULL;
//}
//
//void RTree::rangeQuery(float * mbr, SortedLinList * res, long type)
//{
//	load_root();
//
//	root_ptr -> rangeQuery(mbr,res,type);
//
//	delete root_ptr;
//	root_ptr = NULL;
//}
////------------------------------------------------------------
//void RTree::read_header(char *buffer)
//{
//	int i;
//
//	memcpy(&dimension, buffer, sizeof(dimension));
//	i = sizeof(dimension);
//
//	memcpy(&num_of_data, &buffer[i], sizeof(num_of_data));
//	i += sizeof(num_of_data);
//
//	memcpy(&num_of_dnodes, &buffer[i], sizeof(num_of_dnodes));
//	i += sizeof(num_of_dnodes);
//
//	memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
//	i += sizeof(num_of_inodes);
//
//	memcpy(&root_is_data, &buffer[i], sizeof(root_is_data));
//	i += sizeof(root_is_data);
//
//	memcpy(&root, &buffer[i], sizeof(root));
//	i += sizeof(root);
//}
////------------------------------------------------------------
//void RTree::write_header(char *buffer)
//{
//	int i;
//
//	memcpy(buffer, &dimension, sizeof(dimension));
//	i = sizeof(dimension);
//
//	memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
//	i += sizeof(num_of_data);
//
//	memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
//	i += sizeof(num_of_dnodes);
//
//	memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
//	i += sizeof(num_of_inodes);
//
//	memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
//	i += sizeof(root_is_data);
//
//	memcpy(&buffer[i], &root, sizeof(root));
//	i += sizeof(root);
//}
//
///*****************************************************************
//this function gets the score of an entry with respect to a linear
//preference function
//para:
//e: the leaf entry
//weight: the prefrence vector
//dimension: dimensionality
//*****************************************************************/
//
//float RTree::get_score_linear(Entry *_e, float *_weight, int _dimension)
//{
//	float max_score=-MAXREAL; //this is the score to be returned
//	//	float *tuple=new float[dimension]; //this is the instantiation 
//	//of a tuple (from one of the cornor points of the MBR)
//
//	//decide the final counter when the for loop ends-------------
//	float end_cnt=1;
//	for (int k=0; k<_dimension; k++)
//		end_cnt*=2;
//	//now we implement a for loop with dimension layers-----------
//	for (int i=0; i<end_cnt; i++)
//	{
//		//now instantiate the tuple and get its score-------------
//		int modseed=1;
//		float score=0;
//		for (int j=0; j<_dimension; j++)
//		{
//			int sub=i/modseed%2;
//			//			tuple[j]=_e[2*j+sub];
//			score+=_e->bounces[2*j+sub]*_weight[j];
//
//			modseed*=2;
//		}
//		//update the max_score if necessary-----------------------
//		if (score>max_score)
//			max_score=score;
//		//--------------------------------------------------------
//	}
//	//------------------------------------------------------------
//
//	//	delete [] tuple;
//	return max_score;
//}
//
//float RTree::get_score_linear(float *_mbr, float *_weight, int _dimension)
//{
//	float max_score=-MAXREAL; //this is the score to be returned
//	//	float *tuple=new float[dimension]; //this is the instantiation 
//	//of a tuple (from one of the cornor points of the MBR)
//
//	//decide the final counter when the for loop ends-------------
//	float end_cnt=1;
//	for (int k=0; k<_dimension; k++)
//		end_cnt*=2;
//	//now we implement a for loop with dimension layers-----------
//	for (int i=0; i<end_cnt; i++)
//	{
//		//now instantiate the tuple and get its score-------------
//		int modseed=1;
//		float score=0;
//		for (int j=0; j<_dimension; j++)
//		{
//			int sub=i/modseed%2;
//			//			tuple[j]=_e[2*j+sub];
//			score+=_mbr[2*j+sub]*_weight[j];
//
//			modseed*=2;
//		}
//		//update the max_score if necessary-----------------------
//		if (score>max_score)
//			max_score=score;
//		//--------------------------------------------------------
//	}
//	//------------------------------------------------------------
//
//	//	delete [] tuple;
//	return max_score;
//}
//
//float RTree::get_score_linear(Entry *_e, float *_weight, int _dimension, 
//	float SCORE_FUNC(const float *, const float *, int))
//{
//	float max_score=-MAXREAL; //this is the score to be returned
//	float *tuple=new float[dimension]; //this is the instantiation 
//	//of a tuple (from one of the cornor points of the MBR)
//
//	//decide the final counter when the for loop ends-------------
//	float end_cnt=1;
//	for (int k=0; k<_dimension; k++)
//		end_cnt*=2;
//	//now we implement a for loop with dimension layers-----------
//	for (int i=0; i<end_cnt; i++)
//	{
//		//now instantiate the tuple and get its score-------------
//		int modseed=1;
//		float score=0;
//		for (int j=0; j<_dimension; j++)
//		{
//			int sub=i/modseed%2;
//			tuple[j]=_e->bounces[2*j+sub];
//
//			modseed*=2;
//		}
//		score=SCORE_FUNC(tuple, _weight, _dimension);
//		//update the max_score if necessary-----------------------
//		if (score>max_score)
//			max_score=score;
//		//--------------------------------------------------------
//	}
//	//------------------------------------------------------------
//
//	delete [] tuple;
//	return max_score;
//}
//
///*****************************************************************
//this function gets the score range of an entry with respect to a linear
//preference function
//para:
//e: the leaf entry
//weight: the prefrence vector
//dimension: dimensionality
//range: the range returned
//*****************************************************************/
//
//float RTree::get_score_range (Entry *_e, float *_weight, int _dimension, float *_range)
//{
//	float max_score=-MAXREAL; //this is the score to be returned
//	float min_score=MAXREAL;
//	//	float *tuple=new float[dimension]; //this is the instantiation 
//	//of a tuple (from one of the cornor points of the MBR)
//
//	//decide the final counter when the for loop ends-------------
//	float end_cnt=1;
//	for (int k=0; k<_dimension; k++)
//		end_cnt*=2;
//	//now we implement a for loop with dimension layers-----------
//	for (int i=0; i<end_cnt; i++)
//	{
//		//now instantiate the tuple and get its score-------------
//		int modseed=1;
//		float score=0;
//		for (int j=0; j<_dimension; j++)
//		{
//			int sub=i/modseed%2;
//			//			tuple[j]=_e[2*j+sub];
//			score+=_e->bounces[2*j+sub]*_weight[j];
//
//			modseed*=2;
//		}
//		//update the max/min_score if necessary-------------------
//		if (score>max_score)
//			max_score=score;
//		if (score<min_score)
//			min_score=score;
//		//--------------------------------------------------------
//	}
//	//------------------------------------------------------------
//
//	//	delete [] tuple;
//	_range[0]=min_score; _range[1]=max_score;
//	return max_score;
//}
//
///*****************************************************************
//this function performs a rank-inquiry query with linear function using
//the breadth-first search.
//para:
//weight: an array (with size equal to the dimensionality) storing
//the query vector
//qscore: the query score
//rslt: the link list storing the ids of the returned tuples (disabled
//for the time being)
//
//Coded by Yufei Tao 05/04/02
//*****************************************************************/
//
//void RTree::rank_qry_inquiry(float *_weight, float _qscore, int *_rslt)
//{
//	//clear the node access summary-----------------------------------
//	for (int i=0; i<10; i++)
//		na[i]=0;
//	//----------------------------------------------------------------
//	load_root();
//	root_ptr->rank_qry_inquiry(_weight, _qscore, _rslt);
//	delete root_ptr;
//	root_ptr=NULL;
//}
//
///*****************************************************************
//this function performs a constrain ranked query with linear function using
//the breadth-first search.
//para:
//weight: an array (with size equal to the dimensionality) storing
//the query vector
//qmbr: the query mbr
//k: top "k" query
//hp: the heap for the breadth-first search
//rslt: the link list storing the ids of the returned tuples
//
//Coded by Yufei Tao 05/04/02
//*****************************************************************/
//
//void RTree::rank_qry_constr_linear(float *_weight, float *_qmbr, int _k, Heap *_hp, int *_rslt)
//{
//	int found_cnt=0;
//
//	//clear the node access summary-----------------------------------
//	for (int i=0; i<10; i++)
//		na[i]=0;
//	//----------------------------------------------------------------
//	//first insert the root into the heap-------------------------
//	HeapEntry *he=new HeapEntry();
//	he->son1=root;
//	he->key=0;
//	he->level=1; //this is not important but make sure it is greather than 0
//	_hp->insert(he);
//	delete he;
//	//------------------------------------------------------------
//
//	while(_hp->used>0)
//	{
//		//remove the top heap-------------------------------------
//		HeapEntry *he=new HeapEntry();
//		_hp->remove(he);
//		int son=he->son1;
//		int level=he->level;
//		delete he;
//		//--------------------------------------------------------
//
//		if (level==0)
//		{
//			_rslt[found_cnt]=son;
//			found_cnt++;
//			if (found_cnt==_k)
//				return;
//		}
//		else  //non-leaf entry
//		{
//			RTNode *child=new RTNode(this, son);
//			for (int i=0; i<child->num_entries; i++)
//			{
//				float *ovrp=overlapRect(dimension, child->entries[i].bounces, _qmbr);
//				if (ovrp)
//				{
//					//evaluate the score of this entry----------------
//					//using the intersection area
//					float score=get_score_linear(ovrp, _weight, dimension);
//					delete []ovrp;
//					//now init a new heap entry-----------------------
//					HeapEntry *he=new HeapEntry();
//					he->son1=child->entries[i].son;
//					he->level=child->level;
//					he->key=-score; 
//					//since the heap sorts the entries in increasing order,
//					//and we need to visit the entries by decreasing order 
//					//of their scores, we take the negative.
//					_hp->insert(he);
//					delete he;
//					//------------------------------------------------
//				}
//			}
//			na[child->level]++;
//			delete child;
//		}
//	}
//}
//
///*****************************************************************
//this function performs a ranked query with linear function using
//the breadth-first search.
//para:
//weight: an array (with size equal to the dimensionality) storing
//the query vector
//k: top "k" query
//hp: the heap for the breadth-first search
//rslt: the link list storing the ids of the returned tuples
//
//Coded by Yufei Tao 20/03/02
//*****************************************************************/
//
//void RTree::rank_qry_linear(float *_weight, int _k, Heap *_hp, int *_rslt)
//{
//	int found_cnt=0;
//
//	//clear the node access summary-----------------------------------
//	for (int i=0; i<10; i++)
//		na[i]=0;
//	//----------------------------------------------------------------
//	//first insert the root into the heap-------------------------
//	HeapEntry *he=new HeapEntry();
//	he->son1=root;
//	he->key=0;
//	he->level=1; //this is not important but make sure it is greather than 0
//	_hp->insert(he);
//	delete he;
//	//------------------------------------------------------------
//
//	while(_hp->used>0)
//	{
//		//remove the top heap-------------------------------------
//		HeapEntry *he=new HeapEntry();
//		_hp->remove(he);
//		int son=he->son1;
//		int level=he->level;
//		delete he;
//		//--------------------------------------------------------
//
//		if (level==0)
//		{
//			_rslt[found_cnt]=son;
//			found_cnt++;
//			if (found_cnt==_k)
//				return;
//		}
//		else
//		{
//			RTNode *child=new RTNode(this, son);
//			for (int i=0; i<child->num_entries; i++)
//			{
//				//evaluate the score of this entry----------------
//				float score=get_score_linear(&(child->entries[i]), _weight, dimension);
//				//now init a new heap entry-----------------------
//				HeapEntry *he=new HeapEntry();
//				he->son1=child->entries[i].son;
//				he->level=child->level;
//				he->key=-score; 
//				//since the heap sorts the entries in increasing order,
//				//and we need to visit the entries by decreasing order 
//				//of their scores, we take the negative.
//				_hp->insert(he);
//				delete he;
//				//------------------------------------------------
//			}
//			na[child->level]++;
//			delete child;
//		}
//	}
//}
//
///*****************************************************************
//this function performs a ranked query with a monotone function using
//the breadth-first search.
//para:
//weight: an array (with size equal to the dimensionality) storing
//the query vector
//k: top "k" query
//hp: the heap for the breadth-first search
//rslt: the link list storing the ids of the returned tuples
//SCORE_FUNC: the score function used to evaluate a tuple
//
//Coded by Yufei Tao 20/03/02
//*****************************************************************/
//
//void RTree::rank_qry_monotone(float *_weight, int _k, Heap *_hp, int *_rslt, 
//	float SCORE_FUNC(const float *, const float *, int))
//{
//	int found_cnt=0;
//
//	//clear the node access summary-----------------------------------
//	for (int i=0; i<10; i++)
//		na[i]=0;
//	//----------------------------------------------------------------
//	//first insert the root into the heap-------------------------
//	HeapEntry *he=new HeapEntry();
//	he->son1=root;
//	he->key=0;
//	he->level=1; //this is not important but make sure it is greather than 0
//	_hp->insert(he);
//	delete he;
//	//------------------------------------------------------------
//
//	while(_hp->used>0)
//	{
//		//remove the top heap-------------------------------------
//		HeapEntry *he=new HeapEntry();
//		_hp->remove(he);
//		int son=he->son1;
//		int level=he->level;
//		delete he;
//		//--------------------------------------------------------
//
//		if (level==0)
//		{
//			_rslt[found_cnt]=son;
//			found_cnt++;
//			if (found_cnt==_k)
//				return;
//		}
//		else
//		{
//			RTNode *child=new RTNode(this, son);
//			for (int i=0; i<child->num_entries; i++)
//			{
//				//evaluate the score of this entry----------------
//				float score=get_score_linear(&(child->entries[i]), _weight, dimension, SCORE_FUNC);
//				//now init a new heap entry-----------------------
//				HeapEntry *he=new HeapEntry();
//				he->son1=child->entries[i].son;
//				he->level=child->level;
//				he->key=-score; 
//				//since the heap sorts the entries in increasing order,
//				//and we need to visit the entries by decreasing order 
//				//of their scores, we take the negative.
//				_hp->insert(he);
//				delete he;
//				//------------------------------------------------
//			}
//			na[child->level]++;
//			delete child;
//		}
//	}
//}
//
///*****************************************************************
//this function performs a kNN search using the best-first algo.
//para:
//q: the query point
//k: top "k" query
//hp: the heap for the best-first search
//fardist: the farthest distance from the query to the k-th NN
//
//Coded by Yufei Tao 16/12/02
//*****************************************************************/
///* Modified by Jeff, Dec 2006 */
///*
//void RTree::kNN(float *_q, int _k, Heap *_hp, float &fardist, vector<CInfoPoint>& _rslt,long query_type)
//{
//fardist=MAXREAL;
//int found_cnt=0;
//
////clear the node access summary-----------------------------------
//for (int i=0; i<10; i++)
//na[i]=0;
////----------------------------------------------------------------
//
////first insert the root into the heap-------------------------
//HeapEntry *he=new HeapEntry();
//he->son1=root;
//he->key=0;
//he->father = -1;
//he->level=1; //this is not important but make sure it is greather than 0
//he->init_HeapEntry(dimension);
//for (int j = 0; j < 2 * dimension; j ++) 
//he->bounces[j] = 0.0;	
//_hp->insert(he);
//// delete he;
//// ------------------------------------------------------------
//
//while(_hp->used>0)
//{
////remove the top heap-------------------------------------
//HeapEntry *he=new HeapEntry();
//he->init_HeapEntry(dimension);
////for (int j = 0; j < 2 * dimension; j ++) 
////	he->bounces[j] = 0.0;
//_hp->remove(he);
//int son=he->son1;
//int level=he->level;
//if (level==0)
//fardist=he->key;
//// delete he;
//// --------------------------------------------------------
//
//if (level==0)
//{
//RTNode *child=new RTNode(this, he->father);
//for(int i = 0; i < child->num_entries; i++)
//{
//if (child->entries[i].son == son)
//{
//if (query_type == -1)
//{
//CInfoPoint point;					
//point.id = son;
//point.x = child->entries[i].bounces[0];
//point.y = child->entries[i].bounces[dimension];
//point.type = child->entries[i].type;
//_rslt.push_back(point);
//found_cnt++;
//break;
//}
//else if (query_type == child->entries[i].type)
//{					
//CInfoPoint point;					
//point.id = son;
//point.x = child->entries[i].bounces[0];
//point.y = child->entries[i].bounces[dimension];
//point.type = child->entries[i].type;
//_rslt.push_back(point);
//found_cnt++;
//break;
//}
//}
//}
////_rslt.push_back(son);
////found_cnt++;
//delete child;
//delete he;
//if (found_cnt==_k)
//return;
//}
//else
//{
//delete he;
//RTNode *child=new RTNode(this, son);
//for (int i=0; i<child->num_entries; i++)
//{
////now init a new heap entry-----------------------
//HeapEntry *he=new HeapEntry();
//he->son1=child->entries[i].son;
//he->level=child->level;
//he->key=MINDIST(_q, child->entries[i].bounces, dimension); 
//he->father = son;
//he->init_HeapEntry(dimension);
//for (int j = 0; j < 2 * dimension; j ++) 
//he->bounces[j] = 0.0;
//_hp->insert(he);
//delete he;
////------------------------------------------------
//}
//na[child->level]++;
//delete child;
//}
//}
//}
//*/
//
//
////////////////////////////////////////////////////////////////////////////
//// functions added for skyline
////////////////////////////////////////////////////////////////////////////
//
////=============================================================
//// BBS for computing the skyline in a subspace
//
//// _hp: heap to be used
//// _rslt: result skyline list to be used
//// _active_dim: a boolean array whose size equals the dimensionality of the R-tree
////              the i-th element equals 'true' if the i-th dimension is involved in the subspace
//
//// this function was based on the bbs algorithm written by Greg, and later modified by yufei tao.
//// coded by yufei tao, may 2005
////=============================================================
//
//void RTree::bbs_subspace(Heap *_hp, float *_rslt, int &_rsltcnt, bool *_active_dim)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
//	for (j = 0; j < 2 * dimension; j ++) 
//		he->bounces[j] = 0.0;
//	_hp->insert(he);
//	//	delete he;
//	//------------------------------------------------------------
//
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//		//		HeapEntry *he = new HeapEntry();
//		//		he->init_HeapEntry(dimension); 
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces = new float[dimension];
//		for (int i = 0; i < dimension; i ++)
//		{
//			if (_active_dim[i])
//				bounces[i] = he->bounces[2 * i];
//			else
//				bounces[i] = 0;
//		}
//
//		//		delete he;
//		//--------------------------------------------------------
//
//		if (level==0)
//		{
//			if (!sky_dom(dimension, _rslt, _rsltcnt, bounces))
//			{
//				memcpy(&_rslt[dimension * _rsltcnt], bounces, sizeof(float) * dimension);
//				_rsltcnt ++;
//			}
//		}
//		else
//		{
//			if (!sky_dom(dimension, _rslt, _rsltcnt, bounces))
//			{
//				RTNode *child = new RTNode(this, son);
//
//				for (i = 0; i < child->num_entries; i ++)
//				{
//					/*
//					if (child -> level == 1)
//					{
//					int cnt = 0;
//					for (j = 0; j < dimension; j ++)
//					if (_active_dim[j])
//					{
//					printf("dim: %d: [%.1f, %.1f]\n", cnt, child->entries[i].bounces[2 * j], 
//					child->entries[i].bounces[2 * j + 1]);
//					cnt++;
//					}
//					printf("\n");
//					}
//					*/
//					float *b = new float[dimension];
//					for (j = 0; j < dimension; j ++) 
//					{
//						if (_active_dim[j])
//							b[j] = child->entries[i].bounces[2 * j];
//						else
//							b[j] = 0;
//					}
//
//					if (!sky_dom(dimension, _rslt, _rsltcnt, b))
//					{
//						//HeapEntry *he = new HeapEntry();
//						he->dim = dimension;
//						he->son1 = child->entries[i].son;
//						he->level = child->level;
//
//						for (int j = 0; j < dimension; j ++)
//							if (_active_dim[j])
//							{
//								he->key = b[j]; break;
//							}
//
//							//he->bounces = new float[2 * dimension];
//							for (j = 0; j < dimension; j++)
//								he->bounces[2 * j] = b[j]; 
//
//							_hp->insert(he);
//							//delete he;
//							//------------------------------------------------
//					}
//					delete [] b;
//				}
//				delete child;
//			}
//		}
//		delete [] bounces;
//	}
//	return;
//}
////=============================================================
////this function checks whether the new point is dominated by any 
////skyline point already found
//
////dim:
////rslt: an array with size dim * rsltcnt. the first dim numbers capture
////      the coordinates of the first point, the next dim numbers for the
////      the next point, and so on. 
////rsltcnt: number of skyline points so far 
////pt: the point being tested
//
////coded by yufei tao, may 2005
////=============================================================
//
//bool dominate(int _dim, float *_pt1, float *_pt2)
//{
//	bool ret = true;
//	for (int i = 0; i < _dim; i ++)
//	{
//		if (_pt1[i] > _pt2[i])
//		{
//			ret = false;
//			break;
//		}
//	}
//
//	return ret;
//}
//
//bool RTree::sky_dom(int _dim, float *_rslt, int _rsltcnt, float *_pt)
//{
//	bool ret = false; 
//	for (int i = 0; i < _rsltcnt; i ++)
//	{
//		float *s_pt = &_rslt[i * _dim];
//		if (dominate(_dim, s_pt, _pt))
//		{
//			ret = true; break;
//		}
//	}
//	return ret;
//}
//
////=============================================================
//// BBS for computing the skyline
//
//// _hp: heap to be used
//// _rslt: result skyline list to be used
//// _rsltcnt: number of skyline points
//
//// this function was modified from the bbs_subspace algorithm written by yufei tao.
//// modified by Junfeng Hu, Sep 27, 2007
////=============================================================
//void RTree::bbs(Heap *_hp, float *_rslt, int &_rsltcnt)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
//	for (j = 0; j < 2 * dimension; j ++) 
//		he->bounces[j] = 0.0;
//	_hp->insert(he);
//
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces = new float[dimension];
//		for (int i = 0; i < dimension; i ++)
//			bounces[i] = he->bounces[2 * i];
//
//		if (level==0)
//		{
//			if (!sky_dom(dimension, _rslt, _rsltcnt, bounces))
//			{
//				memcpy(&_rslt[dimension * _rsltcnt], bounces, sizeof(float) * dimension);
//				_rsltcnt ++;
//			}
//		}
//		else
//		{
//			if (!sky_dom(dimension, _rslt, _rsltcnt, bounces))
//			{
//				RTNode *child = new RTNode(this, son);
//
//				for (i = 0; i < child->num_entries; i ++)
//				{
//					float *b = new float[dimension];
//					for (j = 0; j < dimension; j ++) 
//						b[j] = child->entries[i].bounces[2 * j];
//
//					if (!sky_dom(dimension, _rslt, _rsltcnt, b))
//					{
//						he->dim = dimension;
//						he->son1 = child->entries[i].son;
//						he->level = child->level;
//						he->key = sky_mindist(dimension, b);
//
//						for (j = 0; j < 2*dimension; j++)
//							he->bounces[j] = child->entries[i].bounces[j]; 
//
//						_hp->insert(he);
//					}
//					delete [] b;
//				}
//				delete child;
//			}
//		}
//		delete [] bounces;
//	}
//	delete he;
//}
//
//float sky_mindist(int _dim, float *_bounces)
//{
//	float ret = 0.0;
//	for (int i = 0; i < _dim; i++)
//		ret += _bounces[i];
//
//	return ret;
//}
//
//void RTree::bbs_topk(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, float* weighting, int k)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
//	for (j = 0; j < 2 * dimension; j ++) 
//		he->bounces[j] = 0.0;
//	_hp->insert(he);
//
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces = new float[dimension*2];
//		for (int i = 0; i < 2 *dimension; i ++)
//			bounces[i] = he->bounces[i];
//
//		if (level==0)
//		{
//			if ( _rsltcnt < k )
//			{
//				memcpy(&_rslt[dimension * _rsltcnt], bounces, sizeof(float) * dimension);
//				_rsltcnt ++;
//			}
//
//			if (_rsltcnt == k)			
//				break;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//
//			for (i = 0; i < child->num_entries; i ++)
//			{
//				he->dim = dimension;
//				he->son1 = child->entries[i].son;
//				he->level = child->level;
//
//				memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//
//				for (j = 0; j < dimension; j ++)
//				{
//					he->key += weighting[j] * he->bounces[j];
//				}					
//
//				_hp->insert(he);
//			}
//			delete child;
//
//		}
//		delete [] bounces;
//	}
//	delete he;
//}
//
//
////void RTree::bbs_topkth(Heap *_hp, float *_rslt,  float* _qry_pt, float* weighting, int k)
////{
////	int i, j;
////	int _rsltcnt = 0;
////
////	//first insert the root into the heap-------------------------
////	_hp->used=0;
////	HeapEntry *he=new HeapEntry();
////	he->dim = dimension;
////	he->son1 = root;
////	he->key = 0;
////	he->level = 1; //this is not important but make sure it is greather than 0
////	he->bounces = new float[2 * dimension];
////	for (j = 0; j < 2 * dimension; j ++) 
////		he->bounces[j] = 0.0;
////	_hp->insert(he);
////
////	while(_hp->used > 0)
////	{
////		//remove the top heap-------------------------------------
////
////		_hp->remove(he);
////		int son = he->son1;
////		int level = he->level;
////		float *bounces = new float[dimension*2];
////		for (int i = 0; i < 2 *dimension; i ++)
////			bounces[i] = he->bounces[i];
////
////		if (level==0)
////		{
////			if ( _rsltcnt < k )
////			{
////				//memcpy(&_rslt[dimension * _rsltcnt], bounces, sizeof(float) * dimension);
////				_rsltcnt ++;
////			}
////
////			if (_rsltcnt == k)
////			{
////				memcpy(_rslt, bounces, sizeof(float) * dimension);
////				break;
////			}
////		}
////		else
////		{
////			RTNode *child = new RTNode(this, son);
////
////			for (i = 0; i < child->num_entries; i ++)
////			{
////				he->dim = dimension;
////				he->son1 = child->entries[i].son;
////				he->level = child->level;
////
////				memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
////
////				for (j = 0; j < dimension; j ++)
////				{
////					he->key += weighting[j] * he->bounces[j];
////				}					
////
////				_hp->insert(he);
////			}
////			delete child;
////
////		}
////		delete [] bounces;
////	}
////	delete he;
////}
//
//void RTree::bbs_topkth(Heap *_hp, float *_rslt,  float* _qry_pt, float* weighting, int k)
//{
//	int i, j;
//	int _rsltcnt = 0;
//
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
//	for (j = 0; j < 2 * dimension; j ++) 
//		he->bounces[j] = 0.0;
//	_hp->insert(he);
//
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces = new float[dimension*2];
//		for (int i = 0; i < 2 *dimension; i ++)
//			bounces[i] = he->bounces[i];
//
//		if (level==0)
//		{
//
//			_rsltcnt ++;
//
//			if (_rsltcnt == k)
//			{
//				for (i = 0; i < dimension; i++)
//					_rslt[i] = bounces[2*i];
//				break;
//			}
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//
//			for (i = 0; i < child->num_entries; i ++)
//			{
//				he->dim = dimension;
//				he->son1 = child->entries[i].son;
//				he->level = child->level;
//
//				he->key=0;
//
//				memcpy(he->bounces, child->entries[i].bounces, 2*(sizeof(float)* dimension));
//
//				for (j = 0; j < dimension; j ++)
//				{
//					he->key += weighting[j] * he->bounces[j * 2];
//				}					
//
//				_hp->insert(he);
//			}
//			delete child;
//
//		}
//		delete [] bounces;
//	}
//	delete he;
//}
//
//
///*
//_rslt：不可比点集合
//个数
//查询点
//一个矢量
//在不可比的点集中，查询点在此向量下的rank。   ranking_querypoint + bedominatnum + 1 = 查询点在整个数据集中于此矢量下的rank
//*/
//void RTree::bbs_ranking_querypoint( float *_rslt, int &_rsltcnt, float* _qry_pt, float* weighting, int &ranking_querypoint)
//{
//	int i, j;
//	ranking_querypoint = 0;
//	//---------------
//
//	float score_querypoint = 0.0,score_incompoint;
//	for (j = 0; j < dimension; j ++)
//	{
//		score_querypoint += weighting[j] * _qry_pt[j];
//	}	
//
//	for (j = 0; j < _rsltcnt; j ++)
//	{
//		score_incompoint=0;
//
//		for (i = 0; i < dimension; i ++)
//		{
//			score_incompoint += weighting[i] * _rslt[j*dimension+i];
//		}
//	
//		if(score_incompoint<score_querypoint)
//			ranking_querypoint++;
//
//
//	}
//
//
//
//}
//
//
//void RTree::traverse_points(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int &bedominatnum)
//{
//	int i, j;
//	_rsltcnt = 0;
//	bedominatnum=0;
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
//	for (j = 0; j < 2 * dimension; j ++) 
//		he->bounces[j] = 0.0;
//	_hp->insert(he);
//	ofstream outfile("allData.txt");
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces_left = new float[dimension];
//		float *bounces_right = new float[dimension];
//
//		if (level==0)
//		{
//			for (int i = 0; i <  dimension; i ++)
//			{
//				bounces_left[i] = he->bounces[ i*2];
//				bounces_right[i] = he->bounces[ i*2 +1];
//			}
//			outfile<<bounces_left[0]<<" "<<bounces_left[1]<<" "<<bounces_left[2]<<endl;
//		}
//		else
//		{
//			
//				RTNode *child = new RTNode(this, son);				
//				for (i = 0; i < child->num_entries; i ++)
//				{
//					
//						he->dim = dimension;
//						he->son1 = child->entries[i].son;
//						he->level = child->level;
//						memcpy(he->bounces, child->entries[i].bounces, 2*(sizeof(float)* dimension));
//						he->key =0;
//						_hp->insert(he);					
//				}
//				delete child;
//		}
//		delete [] bounces_left;
//		delete [] bounces_right;
//	}
//	delete he;
//}
//
//
///*
//hp:
//_rslt:不可比点的集合
//_rsltcnt：不可比点的个数
//_qry_pt：查询点
//bedominatnum：求得的支配q点的点的个数
//*/
//
//void RTree::incomparable_points1(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int &bedominatnum)
//{
//	int i, j;
//	_rsltcnt = 0;
//	bedominatnum=0;
//	int a;
//	float b;
//	float *aa= new float[dimension];
//	
//	ifstream in("Independent_2d_100k.txt");
//
//	for(i=0;i<100000;i++)
//	{
//		in>>a;
//		for(j=0;j<dimension;j++)
//		{
//			in>>b;
//			aa[j]=b;
//		}
//
//		if (dominate(dimension, aa, _qry_pt))//if (!dominate(dimension, bounces_right, _qry_pt)&&!dominate(dimension,_qry_pt , bounces_left))
//		{
//			bedominatnum++;
//		}
//		else if(!dominate(dimension, _qry_pt, aa))
//		{
//			memcpy(&_rslt[dimension * _rsltcnt], aa, sizeof(float) * dimension);
//			_rsltcnt ++;
//		}
//
//	}	
//}
//
//void RTree::incomparable_points(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int &bedominatnum)
//{
//	int i, j;
//	_rsltcnt = 0;
//	bedominatnum=0;
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
// 
//	//------------------------------------------------
//	memset(he->bounces, 0.0, 2 * dimension * sizeof(float) );
//	/*for (j = 0; j < 2 * dimension; j ++) 
//		he->bounces[j] = 0.0;*/
//	//------------------------------------------------
//	_hp->insert(he);
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces_left = new float[dimension];
//		float *bounces_right = new float[dimension];
//
//		for (i = 0; i <  dimension; i ++)
//		{
//			bounces_left[i] = he->bounces[ i*2];
//			bounces_right[i] = he->bounces[ i*2 +1];
//		}
//	
//		if (level==0)
//		{
//			if (dominate(dimension, bounces_right, _qry_pt))//if (!dominate(dimension, bounces_right, _qry_pt)&&!dominate(dimension,_qry_pt , bounces_left))
//			{
//				bedominatnum++;
//			}
//			else if(!dominate(dimension, _qry_pt, bounces_left))
//			{
//				memcpy(&_rslt[dimension * _rsltcnt], bounces_left, sizeof(float) * dimension);
//				_rsltcnt ++;
//			}
//		}
//		else
//		{
//			if (!dominate(dimension,_qry_pt , bounces_left))
//			{
//				RTNode *child = new RTNode(this, son);
//				float *left = new float[dimension];
//
//				for (i = 0; i < child->num_entries; i ++)
//				{
//
//
//						he->dim = dimension;
//						he->son1 = child->entries[i].son;
//						he->level = child->level;
//						memcpy(he->bounces, child->entries[i].bounces, 2*sizeof(float)* dimension);
//						he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//						_hp->insert(he);
//
//
//				}
//				delete child;
//				delete [] left;
//			}
//		}
//		delete [] bounces_left;
//		delete [] bounces_right;
//	}
//	delete he;
//}
//
//
//void RTree::incomparable_points_reuse(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int &bedominatnum)
//{
//	int i, j;
//	_rsltcnt = 0;
//	bedominatnum=0;
//
//	HeapEntry *he=new HeapEntry();
//	he->bounces = new float[2 * dimension];
//
//	//first time, insert the root into the heap-------------------------
//	if(_hp->used==0)
//	{		
//		he->dim = dimension;
//		he->son1 = root;
//		he->key = 0;
//		he->level = 1; //this is not important but make sure it is greather than 0		
//		//------------------------------------------------
//		memset(he->bounces, 0.0, 2 * dimension * sizeof(float) );
//		_hp->insert(he);
//	}
//
//	//the reuse heap which store the visired nodes -------------------------
//	Heap *hp_reuse = new Heap();
//	hp_reuse->init(dimension);
//
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces_left = new float[dimension];
//		float *bounces_right = new float[dimension];
//
//		for (i = 0; i <  dimension; i ++)
//		{
//			bounces_left[i] = he->bounces[ i*2];
//			bounces_right[i] = he->bounces[ i*2 +1];
//		}
//	
//		if (level==0)
//		{
//			if (dominate(dimension, bounces_right, _qry_pt))//if (!dominate(dimension, bounces_right, _qry_pt)&&!dominate(dimension,_qry_pt , bounces_left))
//			{
//				bedominatnum++;
//			}
//			else if(!dominate(dimension, _qry_pt, bounces_left))
//			{
//				memcpy(&_rslt[dimension * _rsltcnt], bounces_left, sizeof(float) * dimension);
//				_rsltcnt ++;
//			}
//			hp_reuse->insert(he);
//
//		}
//		else
//		{
//			if (!dominate(dimension,_qry_pt , bounces_left))
//			{
//				RTNode *child = new RTNode(this, son);
//				float *left = new float[dimension];
//
//				for (i = 0; i < child->num_entries; i ++)
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					memcpy(he->bounces, child->entries[i].bounces, 2*sizeof(float)* dimension);
//					he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//					_hp->insert(he);
//				}
//				delete child;
//				delete [] left;
//			}
//			else
//				hp_reuse->insert(he);
//		}
//		delete [] bounces_left;
//		delete [] bounces_right;
//	}
//
//	while(hp_reuse->used > 0)
//	{
//		hp_reuse->remove(he);
//		_hp->insert(he);
//	}
//
//	delete he;
//}
//
//
////对每个查询点找出why-not权重矢量，并把  notCount whynot矢量的个数 个矢量返回，要求查询点的rank 为ranks中的一个
///*
//* Description: 判断一个点是否满足条件：在权重向量中存在 count 个向量，使得该点是这些向量的第top-k' 
//*/
//float* RTree::checkQueryPoint(Heap *_hp, float* _qry_pt, float* weightSet, int weightCount, int dimension, int missCnt, int* ranks, char* dataSetName)
//{
//	int count = 0;
//	int incomNum, bedominatnum, rank;
//	float* incompPoints = new float[3000000];
//	float* whynWeightSet = new float[missCnt * dimension];
//	// 找出不可比的点集
//	int _rsltcnt;
//	incomparable_points( _hp, incompPoints, incomNum, _qry_pt, bedominatnum );
//	
//
//	char qFileName[100];
//	memset( qFileName, '\0', 100 );
//	sprintf(qFileName, ".\\%d_%d_%s", *ranks, missCnt, dataSetName);
//	ofstream outfile(qFileName, ios::app);
//	outfile<<"查询点: ";
//	for( int qry_pos = 0; qry_pos < dimension; qry_pos++)
//		outfile<<_qry_pt[qry_pos]<<'\t';
//	outfile<<endl;
//
//	//对一个矢量计算查询点的等级
//	for( int i = 0; i < weightCount; i++ )
//	{  
//		bbs_ranking_querypoint( incompPoints, incomNum, _qry_pt, &weightSet [ i * dimension ], rank);
//		outfile<<"weight: "<<setw(12);
//		for( int w_pos = 0; w_pos < dimension; w_pos++)
//			outfile<<weightSet [ i * dimension + w_pos]<<'\t';
//
//		outfile<< "\trank: "<<rank + bedominatnum + 1<<endl;
//		//如果查询点的rank 为所要求的，说明这个查询点满足
//		int drift = 2;
//		if( ranks[ count ] == 11 )
//			drift = 0;
//		if( rank + bedominatnum + 1 <= ranks[count] + drift &&  rank + bedominatnum + 1 >= ranks[count] - drift)
//		{
//			memcpy( &whynWeightSet[ count * dimension ], &weightSet[ i * dimension ], sizeof(float) * dimension );
//			count++;
//
//			// 若满足要求的向量个数和whynot的个数相同
//			if( count == missCnt )
//			{
//				outfile<<"-------------rank+bedominatnum+1:\t"<<rank+bedominatnum+1<<"---------------------"<<endl;
//				outfile<<"--------------------------ranks[count] :\t"<<ranks[count - 1] <<"---------------------"<<endl;
//				for( int cnt = 0; cnt < count; cnt++ )
//				{
//					outfile<<"weight:";
//					for( int dim = 0; dim < dimension; dim++ )
//						outfile<<" "<<whynWeightSet[ cnt * dimension + dim ];
//					outfile<<endl;
//				}
//				outfile.close();
//				delete incompPoints;
//				return whynWeightSet;
//			}
//		}
//	}
//	outfile<<"--------------------------------------------------------------"<<endl;
//	outfile.close();
//	delete incompPoints;
//	return NULL;
//}
//
//void RTree::SampleWeightsbyq(float quality_of_answer, float probability_guarantee, float *sample, float *quwey_lower, float *quwey_upper )
//{
//	long samplesize;
//	//int dimension = w_origin.weighting.size();
//	//srand ( time(NULL) );
//
//	samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer/100))+1;
//
//	for(int i = 0; i < samplesize; i++)
//	{
//		for(int j = 0; j < dimension; j++)
//		{	
//			int int_r = rand();
//			long base = RAND_MAX-1;
//			float f_r  = ((float) int_r) / base;
//			float temp = (quwey_upper[j]  -  quwey_lower[j]) * f_r +  quwey_lower[j];
//			
// 			//float temp = (float)(quwey_lower[j] + rand() % (int)(quwey_upper[j]  - quwey_lower[j]+1) );
//			sample[i*dimension+j] = temp;
//		}
//		//cout  << endl;
//	}	
//}
//
///*
//* Description: 在满足( 不可比的点 - 查询点 ) * 权重矢量 = 0 的边界上选取权重矢量样本
//* Author: zll
//* Parameters:
//*  sampleSize: 在此平面取样的大小
//*  sample: 存放取样后的矢量，个数为 sampleSize 个
//*  incomparable_point: 与查询点不可比的所有点
//*  _qry_pt: 查询点
//*/
//void RTree::SampleWeights( int sampleSize, float *sample, float* incomPoints, int incomnum, float* _qry_pt)
//{
//	int  samplePos = 0;
//	float* factors = new float[dimension];
//	float* variables = new float[dimension];// 一个等式各变量的取值
//	// count：表示在一个边界上要取的点的个数，最小为1
//	int count =  ceil( ( float ) sampleSize / incomnum );
//	bool flag = true;
//	// 不可比点的个数，即边界数，在每个边界上求得 count 个 向量
//	for( int boundIndex = 0; boundIndex < incomnum && flag; boundIndex++ )
//	{
//		// 计算系数，将两点的各维相减，并排序
//		for( int dim = 0; dim < dimension; dim++ )
//			factors[ dim ] = incomPoints[ boundIndex * dimension + dim ] - _qry_pt[ dim ];
//		sort(factors, factors + dimension );
//
//		if ( count * ( boundIndex + 1 ) > sampleSize )
//		{
//			count = sampleSize - boundIndex * count;
//			flag = false;
//		}
//		// 在边界上每次取一个向量
//		for( int vecIndex = 0; vecIndex < count; vecIndex++ )
//		{
//			float sum = 0.0; // 一个向量各维度的和
//			float x = 0.0;			
//			// 随机生成 变量的取值
//			SampleWeights(variables, 1, dimension - 1, false);
//
//			// 求最后一个变量
//			for( int dim = 0; dim < dimension - 1; dim++ )
//			{
//				// 先计算系数为负数相加的结果
//				if( factors[ dim ] < 0 )
//					x = x - factors[ dim ] * variables[ dim ];
//				else
//				{
//					// variables[ dim ] 的取值应该比 ( 当前所得的结果 / 当前变量对应的系数 ) 小，否则可能导致其它系数为正且未设置的变量只能取负值。
//					if( variables[ dim ] * factors[ dim ] > x )
//						variables[ dim ] = ((float) rand() / RAND_MAX ) * ( x / factors[ dim ] );
//					x = x - variables[ dim ] * factors[ dim ];
//				}
//				sum += variables[ dim ];
//			}
//			variables[ dimension - 1 ] = x / factors[ dimension - 1 ];
//			sum += variables[dimension - 1];
//
//			// 对取样的向量各维进行重新计算，使得各维的和为1
//			for( int dim = 0; dim < dimension; dim++ )
//			{
//				sample[samplePos] = variables[ dim ] /sum;
//				samplePos++;
//			}
//		}
//	}
//	delete factors;
//	delete variables;
//}
//
///*
//* Description: 随机生成一个权重矢量集合, 矢量各维度值的和为1
//* Author: zll
//* Parameters:
//*  sample: 随机生成的权重矢量集合
//*  dimension: 矢量的维度
//*  sampleSize: 矢量个数
//* Return:  
//*/
//void RTree::SampleWeights( float* sample, int dimension, int sampleSize, bool isWeight)
//{
//	for(int i = 0; i < sampleSize; i++)
//	{
//		float sum = 0;
//		for(int j = 0; j < dimension; j++)
//		{			
//			float temp = (float)(rand())/RAND_MAX;
//			sample[i*dimension+j] = temp;
//			sum += temp;
//		}
//		if( isWeight )
//			for(int j = 0; j < dimension; j++)
//				sample[i*dimension+j] = sample[i*dimension+j] /sum ;	
//	}
//}
//
//void RTree::modifyWandK_approximate(Heap *_hp, float *modifyweight, int &modifyk, float &penalty, float* _qry_pt, float* weightset,int &weightnum, int k, float quality_of_answer, float probability_guarantee, float& pruning)
//{//modifyweight, modifyk存放结果，
//	int bedomnum, i, j, incomnum;  //incomnum是incomparable点的个数，bedomnum是dominate查询点q的点的个数，
//	float *incomparable = new float[3000000];
//
//	incomparable_points(_hp,  incomparable, incomnum, _qry_pt, bedomnum);   //求incomparable点的集合，放在数组incomparable中，并求出incomnum和 bedomnum的值
//
//	int qminrank;
//	qminrank = bedomnum + 1;
//	//qmaxrank = incomnum + bedomnum + 1;  //求出查询点q可能的排序的范围
//
//	long samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer/100))+1;    //样本大小
//	float *sampleweight = new float[samplesize * dimension];
//		if( incomnum != 0 )
//		SampleWeights( samplesize, sampleweight, incomparable, incomnum, _qry_pt);
//	//( quality_of_answer,  probability_guarantee,  sampleweight);  //取样，放在数组sampleweight中；
//	else
//		SampleWeights( sampleweight, dimension, samplesize, true);
//	//( quality_of_answer,  probability_guarantee,  sampleweight);  //取样，放在数组sampleweight中；
//
//
//	int *rankofq = new int[samplesize];   //用来记录查询点q在每个样本权重下的排序位置
//	float *currentweight = new float[dimension];
//	int ranking;
//
//	for(i=0; i<samplesize ; i++)   //对于每个样本向量下查询点q的排序
//	{
//		//-------------------------------------------------------------
//		memcpy(currentweight, &sampleweight [ i * dimension ], dimension * sizeof(float) );
//		/*for (j = 0; j < dimension; j ++)
//		{
//			currentweight[j] = sampleweight[ i * dimension + j];
//		}*/
//		//--------------------------------------------------------------
//		bbs_ranking_querypoint( incomparable, incomnum,  _qry_pt,  currentweight, ranking);
//		rankofq[i]=ranking+bedomnum+1;
//	}
//
//	int qmaxrank_real = 0;
//	/*ofstream outfile(".\\parax.txt", ios::app);*/
//	for(i=0; i<weightnum ; i++)   //对于每个why-not向量下查询点q的排序
//	{
//		//-------------------------------------------------------
//	    memcpy( currentweight, &weightset [i * dimension], sizeof(float) * dimension );
//		/*	for (j = 0; j < dimension; j ++)
//		{
//			currentweight[j] = weightset[ i * dimension + j];
//		}*/
//		//-------------------------------------------------------
//		bbs_ranking_querypoint( incomparable, incomnum,  _qry_pt,  currentweight, ranking);
//		/*outfile<<"_qry_pt "<<endl;
//		for( int dd = 0; dd < dimension; dd++ )
//			outfile<<_qry_pt[dd]<<"  ";*/
//		//outfile<<endl;
//		//outfile<<"currentweight"<<currentweight[0]<<endl;
//		//outfile<<"ranking+bedomnum+1   "<<ranking+bedomnum+1<<endl;
//		if(qmaxrank_real<ranking+bedomnum+1)
//			qmaxrank_real = ranking+bedomnum+1;		
//	}
//
//	float parameter_k, parameter_w = 0, ww = 0;
//	parameter_k = 0.5/(qmaxrank_real-k);
//
//	for(i=0; i<weightnum ; i++)   
//	{
//		ww=0;
//		for (j = 0; j < dimension; j ++)
//			ww += weightset[ i * dimension + j] * weightset[ i * dimension + j];
//		
//		ww += 1;
//		ww = sqrt(ww); 
//		parameter_w += ww;			
//	}
//	parameter_w = 0.5/parameter_w;
//
//	/*outfile<<"parameter_w  "<<parameter_w<<endl;*/
//
//	int minrank, id, x;
//	float *temptweight = new float[samplesize * dimension];
//	int *temptrankofq = new int [samplesize]; //根据查询点q在样本权重下排序的位置对样本进行排序，从小到大，排序后的权重放在temptweight，对应的q的排序放在temptrankofq；看看是否可以优化
//
//	for(i=0; i<samplesize ; i++)
//	{
//		minrank=100000;
//		id=0;
//		for (j = 0; j < samplesize; j ++)
//		{
//			if(rankofq[j]!=0 && rankofq[j] < minrank)
//			{
//				id=j;
//				minrank=rankofq[j];
//			}
//		}
//		//--------------------------------------------------------------
//         memcpy(&temptweight [i * dimension], &sampleweight [id * dimension], sizeof(float) * dimension );
//		//for (x = 0; x < dimension; x ++)
//		//{
//		//	temptweight[i * dimension + x] = sampleweight[ id * dimension + x];
//		//}
//		//--------------------------------------------------------------
//
//		temptrankofq[i] = minrank;
//		rankofq[id] = 0;
//
//		//cout <<  temptrankofq[i] << ", ";
//	}
//	//float *modifyweight = new float[1000];
//	float *mindeltaweight = new float[weightnum * dimension];  //用来存放k小于等于某个值时，每个why-not向量delta值最小所对应的样本向量
//	float *mindelta = new float[weightnum];        //用来存放k小于等于某个值时，每个why-not向量最小的delta值
//	//int modifyk;
//	float deltaw, currentpenalty ;
//	bool tag;
//
//	penalty=0;
//
//	//-----------------------------------
//	memset( mindelta, 0, sizeof (float) * weightnum );
//	/*for(i = 0; i < weightnum; i ++)    
//	{
//		mindelta[i] = 0;
//	}*/
//	//-----------------------------------
//
//	//for(i = 0; i < weightnum; i ++)    //用temptrankofq最小的那个向量来初始化modifyweight，modifyk， mindeltaweight和mindelta
//	//{
//	//	//------------------------------------------------------
//	//	memcpy( &modifyweight [i * dimension], temptweight, sizeof( float ) * dimension );
//	//	memcpy( &mindeltaweight[i * dimension], temptweight, sizeof( float ) * dimension );
//	//	//------------------------------------------------------
//	//	for(x = 0; x < dimension; x ++)
//	//	{
//	//		//------------------------------------------------------------------
//	//		//modify0
//	//		//minde1ltaweight[i * dimension + x] = temptweight[ x ];
//	//		//-------------------------------------------------------------------
//	//		mindelta[i] = mindelta[i] + (modifyweight[i * dimension + x] - weightset[i * dimension + x])*(modifyweight[i * dimension + x] - weightset[i * dimension + x]);
//	//		//cout <<  modifyweight[i * dimension + x] << ", ";
//	//	}
//
//	//	//cout << endl;
//
//	//	mindelta[i] = sqrt(mindelta[i]); 
//
//	//	penalty += mindelta[i];
//
//	//}
//
//	//penalty = penalty * parameter_w; 
//
//	//========================================
//	// 以 why -not 向量初始化 penalty
//	modifyk = qmaxrank_real;	
//	if( modifyk > k )
//		penalty = (modifyk - k) * parameter_k; 
//	
//
//	//===================================
//	memcpy( modifyweight, weightset, sizeof( float ) * dimension * weightnum );
//	for( int mindeltaNum = 0; mindeltaNum < weightnum; mindeltaNum++)
//		mindelta[mindeltaNum] =  5.0;
//	//if(modifyk>k)
//	//	penalty += (modifyk - k) * parameter_k; //计算penalty
//
//	//cout << "penalty: " <<  penalty << endl;
//
//
//	//for(int wnum = 0; wnum < weightnum; wnum ++)    //用temptrankofq最小的那个向量来初始化modifyweight，modifyk， mindeltaweight和mindelta
//	//{
//	//	//------------------------------------------------------
//	//	//memcpy( &modifyweight [wnum * dimension], temptweight, sizeof( float ) * dimension );
//	//	//memcpy( &mindeltaweight[wnum * dimension], temptweight, sizeof( float ) * dimension );
//	//	//------------------------------------------------------
//	//	for(x = 0; x < dimension; x ++)
//	//	{
//	//		//------------------------------------------------------------------
//	//		//modify0
//	//		//minde1ltaweight[i * dimension + x] = temptweight[ x ];
//	//		//-------------------------------------------------------------------
//	//		mindelta[wnum] = mindelta[wnum] + (temptweight[wnum * dimension + x] - weightset[wnum * dimension + x])*(temptweight[wnum * dimension + x] - weightset[wnum * dimension + x]);
//	//		//cout <<  modifyweight[i * dimension + x] << ", ";
//	//	}
//	//}
//	for(i=0; i<samplesize ; i++)   //找delta最小的W和k
//	{
//		if(temptrankofq[i] > qmaxrank_real )  //提前终止
//		{
//			//outfile<<"samplesize - i"<<samplesize - i<<endl;
//			pruning = ( samplesize - i ) / (float)samplesize;
//			break;
//		}
//		tag=false;		
//
//		for(j = 0; j < weightnum; j ++)//更新mindeltaweight
//		{
//			deltaw= 0;
//			for(x = 0; x < dimension; x ++)
//			{
//				deltaw = deltaw + (temptweight[ i * dimension + x ] - weightset[j * dimension + x]) * (temptweight[ i * dimension + x ] - weightset[j * dimension + x]);
//			}
//
//			deltaw = sqrt(deltaw); 
//
//			//更新mindeltaweight
//
//			if(deltaw < mindelta[j])
//			{
//				tag=true;
//				mindelta[j] = deltaw;
//
//				//---------------------------------------------------------
//				memcpy( &mindeltaweight [j * dimension], &temptweight [i * dimension], sizeof( float ) * dimension );
//				/*for(x = 0; x < dimension; x ++)
//				{
//					mindeltaweight[j * dimension + x] = temptweight[ i * dimension + x ];					
//				}*/
//				//---------------------------------------------------------
//			}
//		}
//
//		if(tag==false)
//			continue;
//
//		deltaw = 0;
//
//		for(j = 0; j < weightnum; j ++) //如果更新后mindeltaweight包含k等于某个值时对应的向量，需计算penalty，验证其是否是解
//			deltaw = deltaw + mindelta[j];
//
//		float pp;
//		if(temptrankofq[i] > k)
//			pp = deltaw * parameter_w + (temptrankofq[i]-k) * parameter_k;
//		else
//			pp = deltaw * parameter_w;
//
//		if( penalty >  pp)
//		{
//			//-------------------------------------------------------
//			memcpy( modifyweight , mindeltaweight , sizeof( float ) * dimension * weightnum );
//			//for(j = 0; j < weightnum; j ++)
//			//{			
//			//	for(x = 0; x < dimension; x ++)
//			//	{
//			//		modifyweight[j * dimension + x] = mindeltaweight[j * dimension + x];
//			//		//cout <<  modifyweight[j * dimension + x] << ", ";
//			//	}
//			//	
//			//	//cout << endl;
//			//}
//            //------------------------------------------
//			modifyk = temptrankofq[i];
//
//			penalty= pp ;
//			//cout <<    "penalty" <<penalty<< endl;
//		}
//	}
//
//	delete incomparable;
//	delete sampleweight;
//	delete rankofq;
//	delete currentweight;
//	delete temptweight;
//	delete temptrankofq;
//	delete mindeltaweight;
//	delete mindelta;
//}
//
//void RTree::modifyWandK_approximate_reuse(Heap *_hp, float *modifyweight, int &modifyk, float &penalty, float* _qry_pt, float* weightset,int &weightnum, int k, float quality_of_answer, float probability_guarantee, float& pruning)
//{//modifyweight, modifyk存放结果，
//	int bedomnum, i, j, incomnum;  //incomnum是incomparable点的个数，bedomnum是dominate查询点q的点的个数，
//	float *incomparable = new float[3000000];
//
//	incomparable_points_reuse(_hp,  incomparable, incomnum, _qry_pt, bedomnum);   //求incomparable点的集合，放在数组incomparable中，并求出incomnum和 bedomnum的值
//
//	int qminrank;
//	qminrank = bedomnum + 1;
//	//qmaxrank = incomnum + bedomnum + 1;  //求出查询点q可能的排序的范围
//
//	long samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer/100))+1;    //样本大小
//	float *sampleweight = new float[samplesize * dimension];
//		if( incomnum != 0 )
//		SampleWeights( samplesize, sampleweight, incomparable, incomnum, _qry_pt);
//	//( quality_of_answer,  probability_guarantee,  sampleweight);  //取样，放在数组sampleweight中；
//	else
//		SampleWeights( sampleweight, dimension, samplesize, true);
//	//( quality_of_answer,  probability_guarantee,  sampleweight);  //取样，放在数组sampleweight中；
//
//
//	int *rankofq = new int[samplesize];   //用来记录查询点q在每个样本权重下的排序位置
//	float *currentweight = new float[dimension];
//	int ranking;
//
//	for(i=0; i<samplesize ; i++)   //对于每个样本向量下查询点q的排序
//	{
//		//-------------------------------------------------------------
//		memcpy(currentweight, &sampleweight [ i * dimension ], dimension * sizeof(float) );
//		/*for (j = 0; j < dimension; j ++)
//		{
//			currentweight[j] = sampleweight[ i * dimension + j];
//		}*/
//		//--------------------------------------------------------------
//		bbs_ranking_querypoint( incomparable, incomnum,  _qry_pt,  currentweight, ranking);
//		rankofq[i]=ranking+bedomnum+1;
//	}
//
//	int qmaxrank_real = 0;
//	/*ofstream outfile(".\\parax.txt", ios::app);*/
//	for(i=0; i<weightnum ; i++)   //对于每个why-not向量下查询点q的排序
//	{
//		//-------------------------------------------------------
//	    memcpy( currentweight, &weightset [i * dimension], sizeof(float) * dimension );
//		/*	for (j = 0; j < dimension; j ++)
//		{
//			currentweight[j] = weightset[ i * dimension + j];
//		}*/
//		//-------------------------------------------------------
//		bbs_ranking_querypoint( incomparable, incomnum,  _qry_pt,  currentweight, ranking);
//		/*outfile<<"_qry_pt "<<endl;
//		for( int dd = 0; dd < dimension; dd++ )
//			outfile<<_qry_pt[dd]<<"  ";*/
//		//outfile<<endl;
//		//outfile<<"currentweight"<<currentweight[0]<<endl;
//		//outfile<<"ranking+bedomnum+1   "<<ranking+bedomnum+1<<endl;
//		if(qmaxrank_real<ranking+bedomnum+1)
//			qmaxrank_real = ranking+bedomnum+1;		
//	}
//
//	float parameter_k, parameter_w = 0, ww = 0;
//	parameter_k = 0.5/(qmaxrank_real-k);
//
//	for(i=0; i<weightnum ; i++)   
//	{
//		ww=0;
//		for (j = 0; j < dimension; j ++)
//			ww += weightset[ i * dimension + j] * weightset[ i * dimension + j];
//		
//		ww += 1;
//		ww = sqrt(ww); 
//		parameter_w += ww;			
//	}
//	parameter_w = 0.5/parameter_w;
//
//	/*outfile<<"parameter_w  "<<parameter_w<<endl;*/
//
//	int minrank, id, x;
//	float *temptweight = new float[samplesize * dimension];
//	int *temptrankofq = new int [samplesize]; //根据查询点q在样本权重下排序的位置对样本进行排序，从小到大，排序后的权重放在temptweight，对应的q的排序放在temptrankofq；看看是否可以优化
//
//	for(i=0; i<samplesize ; i++)
//	{
//		minrank=100000;
//		id=0;
//		for (j = 0; j < samplesize; j ++)
//		{
//			if(rankofq[j]!=0 && rankofq[j] < minrank)
//			{
//				id=j;
//				minrank=rankofq[j];
//			}
//		}
//		//--------------------------------------------------------------
//         memcpy(&temptweight [i * dimension], &sampleweight [id * dimension], sizeof(float) * dimension );
//		//for (x = 0; x < dimension; x ++)
//		//{
//		//	temptweight[i * dimension + x] = sampleweight[ id * dimension + x];
//		//}
//		//--------------------------------------------------------------
//
//		temptrankofq[i] = minrank;
//		rankofq[id] = 0;
//
//		//cout <<  temptrankofq[i] << ", ";
//	}
//	//float *modifyweight = new float[1000];
//	float *mindeltaweight = new float[weightnum * dimension];  //用来存放k小于等于某个值时，每个why-not向量delta值最小所对应的样本向量
//	float *mindelta = new float[weightnum];        //用来存放k小于等于某个值时，每个why-not向量最小的delta值
//	//int modifyk;
//	float deltaw, currentpenalty ;
//	bool tag;
//
//	penalty=0;
//
//	//-----------------------------------
//	memset( mindelta, 0, sizeof (float) * weightnum );
//	/*for(i = 0; i < weightnum; i ++)    
//	{
//		mindelta[i] = 0;
//	}*/
//	//-----------------------------------
//
//	//for(i = 0; i < weightnum; i ++)    //用temptrankofq最小的那个向量来初始化modifyweight，modifyk， mindeltaweight和mindelta
//	//{
//	//	//------------------------------------------------------
//	//	memcpy( &modifyweight [i * dimension], temptweight, sizeof( float ) * dimension );
//	//	memcpy( &mindeltaweight[i * dimension], temptweight, sizeof( float ) * dimension );
//	//	//------------------------------------------------------
//	//	for(x = 0; x < dimension; x ++)
//	//	{
//	//		//------------------------------------------------------------------
//	//		//modify0
//	//		//minde1ltaweight[i * dimension + x] = temptweight[ x ];
//	//		//-------------------------------------------------------------------
//	//		mindelta[i] = mindelta[i] + (modifyweight[i * dimension + x] - weightset[i * dimension + x])*(modifyweight[i * dimension + x] - weightset[i * dimension + x]);
//	//		//cout <<  modifyweight[i * dimension + x] << ", ";
//	//	}
//
//	//	//cout << endl;
//
//	//	mindelta[i] = sqrt(mindelta[i]); 
//
//	//	penalty += mindelta[i];
//
//	//}
//
//	//penalty = penalty * parameter_w; 
//
//	//========================================
//	// 以 why -not 向量初始化 penalty
//	modifyk = qmaxrank_real;	
//	if( modifyk > k )
//		penalty = (modifyk - k) * parameter_k; 
//	
//
//	//===================================
//	memcpy( modifyweight, weightset, sizeof( float ) * dimension * weightnum );
//	for( int mindeltaNum = 0; mindeltaNum < weightnum; mindeltaNum++)
//		mindelta[mindeltaNum] =  5.0;
//	//if(modifyk>k)
//	//	penalty += (modifyk - k) * parameter_k; //计算penalty
//
//	//cout << "penalty: " <<  penalty << endl;
//
//
//	//for(int wnum = 0; wnum < weightnum; wnum ++)    //用temptrankofq最小的那个向量来初始化modifyweight，modifyk， mindeltaweight和mindelta
//	//{
//	//	//------------------------------------------------------
//	//	//memcpy( &modifyweight [wnum * dimension], temptweight, sizeof( float ) * dimension );
//	//	//memcpy( &mindeltaweight[wnum * dimension], temptweight, sizeof( float ) * dimension );
//	//	//------------------------------------------------------
//	//	for(x = 0; x < dimension; x ++)
//	//	{
//	//		//------------------------------------------------------------------
//	//		//modify0
//	//		//minde1ltaweight[i * dimension + x] = temptweight[ x ];
//	//		//-------------------------------------------------------------------
//	//		mindelta[wnum] = mindelta[wnum] + (temptweight[wnum * dimension + x] - weightset[wnum * dimension + x])*(temptweight[wnum * dimension + x] - weightset[wnum * dimension + x]);
//	//		//cout <<  modifyweight[i * dimension + x] << ", ";
//	//	}
//	//}
//	for(i=0; i<samplesize ; i++)   //找delta最小的W和k
//	{
//		if(temptrankofq[i] > qmaxrank_real )  //提前终止
//		{
//			//outfile<<"samplesize - i"<<samplesize - i<<endl;
//			pruning = ( samplesize - i ) / (float)samplesize;
//			break;
//		}
//		tag=false;		
//
//		for(j = 0; j < weightnum; j ++)//更新mindeltaweight
//		{
//			deltaw= 0;
//			for(x = 0; x < dimension; x ++)
//			{
//				deltaw = deltaw + (temptweight[ i * dimension + x ] - weightset[j * dimension + x]) * (temptweight[ i * dimension + x ] - weightset[j * dimension + x]);
//			}
//
//			deltaw = sqrt(deltaw); 
//
//			//更新mindeltaweight
//
//			if(deltaw < mindelta[j])
//			{
//				tag=true;
//				mindelta[j] = deltaw;
//
//				//---------------------------------------------------------
//				memcpy( &mindeltaweight [j * dimension], &temptweight [i * dimension], sizeof( float ) * dimension );
//				/*for(x = 0; x < dimension; x ++)
//				{
//					mindeltaweight[j * dimension + x] = temptweight[ i * dimension + x ];					
//				}*/
//				//---------------------------------------------------------
//			}
//		}
//
//		if(tag==false)
//			continue;
//
//		deltaw = 0;
//
//		for(j = 0; j < weightnum; j ++) //如果更新后mindeltaweight包含k等于某个值时对应的向量，需计算penalty，验证其是否是解
//			deltaw = deltaw + mindelta[j];
//
//		float pp;
//		if(temptrankofq[i] > k)
//			pp = deltaw * parameter_w + (temptrankofq[i]-k) * parameter_k;
//		else
//			pp = deltaw * parameter_w;
//
//		if( penalty >  pp)
//		{
//			//-------------------------------------------------------
//			memcpy( modifyweight , mindeltaweight , sizeof( float ) * dimension * weightnum );
//			//for(j = 0; j < weightnum; j ++)
//			//{			
//			//	for(x = 0; x < dimension; x ++)
//			//	{
//			//		modifyweight[j * dimension + x] = mindeltaweight[j * dimension + x];
//			//		//cout <<  modifyweight[j * dimension + x] << ", ";
//			//	}
//			//	
//			//	//cout << endl;
//			//}
//            //------------------------------------------
//			modifyk = temptrankofq[i];
//
//			penalty= pp ;
//			//cout <<    "penalty" <<penalty<< endl;
//		}
//	}
//
//	delete incomparable;
//	delete sampleweight;
//	delete rankofq;
//	delete currentweight;
//	delete temptweight;
//	delete temptrankofq;
//	delete mindeltaweight;
//	delete mindelta;
//}
//
//
//
//
//void RTree::modifyQ_accurate(Heap *_hp, float* _qry_pt, float* weightset,int weightnum, int k, float *&qpmodified)
//{
//	int i;
//	float *currentweight = new float[dimension];
//	float *kthpoint = new float[dimension];
//	float *kthpoints = new float[weightnum * dimension]; //存放每个why-not向量下ranking是k的那个点
//
//	for (i = 0; i < weightnum; i ++) //求每个why-not向量下ranking是k的那个点
//	{
//		//-----------------------------------------------
//		memcpy( currentweight, &weightset[i * dimension], dimension * sizeof(float));
//		/*for (j = 0; j < dimension; j ++)
//		{
//			currentweight[j] = weightset[ i * dimension + j];
//		}*/
//		//-----------------------------------------------
//		bbs_topkth(_hp, kthpoint,  _qry_pt,  currentweight,  k);
//		memcpy(&kthpoints[dimension * i], kthpoint, sizeof(float) * dimension);
//	}
////	qpmodified = q_Quadprog(_qry_pt, weightset, kthpoints, weightnum, dimension);
//
//	delete [] currentweight;
//	delete [] kthpoint;
//}
//
//
//void RTree::modifyQ_approximate(Heap *_hp,  float* _qry_pt, float* weightset,int &weightnum, int k, float *qpmodified)
//{
//	int i, j;
//	float *currentweight = new float[dimension];
//	float *kthpoint = new float[dimension];
//	float *kthnearest = new float[dimension];
//
//	//---------------------------------------
//	memcpy( qpmodified, _qry_pt, dimension * sizeof(float));
//	/*for (j = 0; j < dimension; j ++)
//	{			
//		qpmodified[j] = _qry_pt[j];
//	}*/
//	//----------------------------------------
//
//	for (i = 0; i < weightnum; i ++) 
//	{
//		memcpy( currentweight, &weightset[ i * dimension], sizeof(float) * dimension );
//		//----------------------------------------------------
//		/*for (j = 0; j < dimension; j ++)
//		{
//			currentweight[j] = weightset[ i * dimension + j];
//		}*/
//		//------------------------------------------------------
//		bbs_topkth(_hp, kthpoint,  _qry_pt,  currentweight,  k);	
//
//		point_nearest_q( kthnearest,  _qry_pt, currentweight,  kthpoint);//求平面上到查询点q最近的点
//
//		for (j = 0; j < dimension; j ++)
//		{
//			if(kthnearest[j] < qpmodified[j])
//				qpmodified[j] = kthnearest[j];
//		}
//	}
//
//	delete [] currentweight;
//	delete [] kthpoint;
//	delete [] kthnearest;
//}
//
//
////void RTree::point_nearest_q( float *_rslt, float* _qry_pt, float* weight, float *kthpoint)
////{
////
////	int i, j;
////	float para, para1,para2;
////
////	para=0.0;
////	para1=0.0;
////	para2=0.0;
////
////	for(i = 0; i <dimension; i ++)
////	{
////		para1 += weight[i] * (_qry_pt[i] - kthpoint[i]);
////		para2 += weight[i] * weight[i];
////	}
////
////	para = para1/para2;
////
////	for(i = 0; i <dimension; i ++)
////	{
////		_rslt[i] = _qry_pt[i] - weight[i] * para;
////		_rslt[i] = ( _rslt[i] < 0 ? 0 : _rslt[i] );
////	}
////}
//
//void RTree::point_nearest_q( float *_rslt, float* _qry_pt, float* weight, float *kthpoint)
//{
//
//	int i, j;
//	bool tag;
//	float para, para1,para2, para3;;
//	float *point_on_plane = new float[dimension];
//	int *tag1 = new int[dimension];
//
//	tag=false;
//	para=0.0;
//	para1=0.0;
//	para2=0.0;
//	para3=0.0;
//
//	for(i = 0; i <dimension; i ++)
//	{
//		para1 += weight[i] * (_qry_pt[i] - kthpoint[i]);
//		para2 += weight[i] * weight[i];
//
//		tag1[i] =1;
//		para3 += weight[i] *  kthpoint[i];
//	}
//
//	para = para1/para2;
//
//	for(i = 0; i <dimension; i ++)
//	{
//		_rslt[i] = _qry_pt[i] - weight[i] * para;
//		if(_rslt[i] < 0)
//		{
//			_rslt[i]=0;
//			tag1[i]=0;
//			tag = true;
//		}
//	}
//
//	while(tag==true)
//	{
//		para=0.0;
//		para1=0.0;
//		para2=0.0;
//
//		tag=false;
//
//		for(i = 0; i <dimension; i ++)
//			point_on_plane[i]=0.0;
//
//		for(i = 0; i <dimension; i ++)
//		{
//			if(tag1[i]==1)
//			{
//				point_on_plane[i]=para3/weight[i];
//				break;
//			}
//		}
//
//		for(i = 0; i <dimension; i ++)
//		{
//			if(tag1[i]==1)
//			{
//				para1 += weight[i] * (_qry_pt[i] - point_on_plane[i]);
//				para2 += weight[i] * weight[i];
//			}
//		}
//
//		para = para1/para2;
//		for(i = 0; i <dimension; i ++)
//		{
//			if(tag1[i]==1)
//			{
//				_rslt[i] = _qry_pt[i] - weight[i] * para;
//				if(_rslt[i] < 0)
//				{
//					_rslt[i]=0;
//					tag1[i]=0;
//					tag = true;
//				}
//
//			}
//		}
//
//	}
//
//
//
//
//	/*	for(i = 0; i <dimension; i ++)
//	{
//	para1 += weight[i] * (_qry_pt[i] - kthpoint[i]);
//	para2 += weight[i] * weight[i];
//
//	tag1[i] =1;
//	point_on_plane[i] =0.0;
//	para3 += weight[i] *  kthpoint[i];
//	}
//
//	para = para1/para2;
//
//	for(i = 0; i <dimension; i ++)
//	{
//	_rslt[i] = _qry_pt[i] - weight[i] * para;
//	if(_rslt[i] < 0)
//	{
//	_rslt[i]=0;
//	tag1[i]=0;
//	tag = true;
//	}
//	}
//
//
//
//	if(tag==true)
//	{
//	para=0.0;
//	para1=0.0;
//	para2=0.0;
//
//	for(i = 0; i <dimension; i ++)
//	{
//	if(tag1[i]==1)
//	{
//	point_on_plane[i]=para3/weight[i];
//	break;
//	}
//	}
//
//	for(i = 0; i <dimension; i ++)
//	{
//	if(tag1[i]==1)
//	{
//	para1 += weight[i] * (_qry_pt[i] - point_on_plane[i]);
//	para2 += weight[i] * weight[i];
//	}
//	}
//
//	para = para1/para2;
//	for(i = 0; i <dimension; i ++)
//	{
//	if(tag1[i]==1)
//	_rslt[i] = _qry_pt[i] - weight[i] * para;
//	}
//
//	}
//
//	*/
//
//
//	delete [] point_on_plane;
//	delete [] tag1;
//}
//
//void RTree::modify_k_W_Q(Heap *_hp,  float* _qry_pt, float* _qry_pt1, float* weightset, int &weightnum, int k, float *qpmodified, float* modifyweight,int &modifyk, float quality_of_answer, float probability_guarantee, float &penalty, float& pruning)
//{
//	penalty=100000.0;
//
//	float samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer/100))+1;    //样本大小
//
//	float *samplepoint = new float[samplesize * dimension];
//	//SampleWeights( quality_of_answer,  probability_guarantee,  samplepoint);  //取样，放在数组samplepoint中；
//	SampleWeightsbyq( quality_of_answer, probability_guarantee, samplepoint,  _qry_pt1,  _qry_pt );
//
//	float *tempmodifyweight = new float[weightnum * dimension];
//	int tempmodifyk ;
//	float temppenalty;
//	float *tempt_qry_pt = new float[dimension];
//	float pruning_tmp = 0.0;
//	float qqmodify = 0.0;
//	float qqmin = 0.0;
//	for(int i=0;i<samplesize;i++)
//	{
//		//-------------------------------------------------------------------------
//		memcpy(tempt_qry_pt, &samplepoint[i*dimension], dimension * sizeof(float));	
//		//for(j=0;j <dimension; j ++)
//		//{
//		//	tempt_qry_pt[j]=samplepoint[i*dimension+j];
//		//}
//		//-----------------------------------------------------------------------------
//		modifyWandK_approximate(_hp, tempmodifyweight, tempmodifyk, temppenalty,  tempt_qry_pt, weightset,weightnum,  k,  quality_of_answer,  probability_guarantee, pruning);
//		pruning_tmp += pruning;
//
//		for( int q_pos = 0; q_pos < dimension; q_pos++ )
//		{
//			qqmodify += ( _qry_pt[q_pos] - tempt_qry_pt[q_pos] ) * ( _qry_pt[q_pos] - tempt_qry_pt[q_pos] );
//			qqmin += ( _qry_pt[q_pos] - _qry_pt1[q_pos] ) * ( _qry_pt[q_pos] - _qry_pt1[q_pos] );
//		}
//		qqmodify = sqrt( qqmodify );
//		qqmin = sqrt( qqmin );
//		temppenalty = temppenalty * 0.5; 
//		temppenalty += 0.5 * qqmodify / qqmin ;
//
//		if(temppenalty < penalty)
//		{
//			modifyk=tempmodifyk;
//			penalty = temppenalty;
//
//			//-------------------------------------------------------------------------------------------------
//			memcpy(modifyweight,tempmodifyweight, sizeof(float) * dimension*weightnum);
//			//for(j = 0; j < weightnum; j ++)//能否用memcpy(modifyweight,tempmodifyweight, sizeof(float) * dimension*weightnum);
//			//{
//			//	for(x = 0; x < dimension; x ++)
//			//	{
//			//		modifyweight[j * dimension + x] = tempmodifyweight[j * dimension + x];
//			//	}
//
//			//}
//			//-----------------------------------------------------------------------------------------------------
//	
//			//--------------------------------------------------------------------------------------------------------
//			memcpy(qpmodified,tempt_qry_pt, sizeof(float) * dimension);
//			//for(x = 0; x < dimension; x ++)//能否用memcpy(qpmodified,tempt_qry_pt, sizeof(float) * dimension);
//			//{
//			//	qpmodified[x] =tempt_qry_pt[ x];
//			//}
//			//---------------------------------------------------------------------------------------------------
//		}
//	}
//	delete samplepoint;
//	delete tempmodifyweight;
//	delete tempt_qry_pt;
//	pruning = pruning_tmp / samplesize;
//}
//
//void RTree::modify_k_W_Q_reuse(Heap *_hp,  float* _qry_pt, float* _qry_pt1, float* weightset, int &weightnum, int k, float *qpmodified, float* modifyweight,int &modifyk, float quality_of_answer, float probability_guarantee, float &penalty, float& pruning)
//{
//	penalty=100000.0;
//
//	float samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer/100))+1;    //样本大小
//
//	float *samplepoint = new float[samplesize * dimension];
//	//SampleWeights( quality_of_answer,  probability_guarantee,  samplepoint);  //取样，放在数组samplepoint中；
//	SampleWeightsbyq( quality_of_answer, probability_guarantee, samplepoint,  _qry_pt1,  _qry_pt );
//
//	float *tempmodifyweight = new float[weightnum * dimension];
//	int tempmodifyk ;
//	float temppenalty;
//	float *tempt_qry_pt = new float[dimension];
//	float pruning_tmp = 0.0;
//	float qqmodify = 0.0;
//	float qqmin = 0.0;
//	for(int i=0;i<samplesize;i++)
//	{
//		//-------------------------------------------------------------------------
//		memcpy(tempt_qry_pt, &samplepoint[i*dimension], dimension * sizeof(float));	
//		//for(j=0;j <dimension; j ++)
//		//{
//		//	tempt_qry_pt[j]=samplepoint[i*dimension+j];
//		//}
//		//-----------------------------------------------------------------------------
//		modifyWandK_approximate_reuse(_hp, tempmodifyweight, tempmodifyk, temppenalty,  tempt_qry_pt, weightset,weightnum,  k,  quality_of_answer,  probability_guarantee, pruning);
//		pruning_tmp += pruning;
//
//		for( int q_pos = 0; q_pos < dimension; q_pos++ )
//		{
//			qqmodify += ( _qry_pt[q_pos] - tempt_qry_pt[q_pos] ) * ( _qry_pt[q_pos] - tempt_qry_pt[q_pos] );
//			qqmin += ( _qry_pt[q_pos] - _qry_pt1[q_pos] ) * ( _qry_pt[q_pos] - _qry_pt1[q_pos] );
//		}
//		qqmodify = sqrt( qqmodify );
//		qqmin = sqrt( qqmin );
//		temppenalty = temppenalty * 0.5; 
//		temppenalty += 0.5 * qqmodify / qqmin ;
//
//		if(temppenalty < penalty)
//		{
//			modifyk=tempmodifyk;
//			penalty = temppenalty;
//
//			//-------------------------------------------------------------------------------------------------
//			memcpy(modifyweight,tempmodifyweight, sizeof(float) * dimension*weightnum);
//			//for(j = 0; j < weightnum; j ++)//能否用memcpy(modifyweight,tempmodifyweight, sizeof(float) * dimension*weightnum);
//			//{
//			//	for(x = 0; x < dimension; x ++)
//			//	{
//			//		modifyweight[j * dimension + x] = tempmodifyweight[j * dimension + x];
//			//	}
//
//			//}
//			//-----------------------------------------------------------------------------------------------------
//	
//			//--------------------------------------------------------------------------------------------------------
//			memcpy(qpmodified,tempt_qry_pt, sizeof(float) * dimension);
//			//for(x = 0; x < dimension; x ++)//能否用memcpy(qpmodified,tempt_qry_pt, sizeof(float) * dimension);
//			//{
//			//	qpmodified[x] =tempt_qry_pt[ x];
//			//}
//			//---------------------------------------------------------------------------------------------------
//		}
//	}
//	delete samplepoint;
//	delete tempmodifyweight;
//	delete tempt_qry_pt;
//	pruning = pruning_tmp / samplesize;
//}
//
//
//
//
////=============================================================
//// BBS for computing the skyline in a constrained space
//
//// _hp: heap to be used
//// _rslt: result skyline list to be used
//// _rsltcnt: number of skyline points
//
//// this function was modified from the bbs_subspace algorithm written by yufei tao.
//// modified by Junfeng Hu, Sep 27, 2007
////=============================================================
//void RTree::bbs_constrained(Heap *_hp, float *_rslt, int &_rsltcnt, float *_bounces)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
//	for (j = 0; j < 2 * dimension; j ++) 
//		he->bounces[j] = 0.0;
//	_hp->insert(he);
//	//	delete he;
//	//------------------------------------------------------------
//
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//		//		HeapEntry *he = new HeapEntry();
//		//		he->init_HeapEntry(dimension); 
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//		float *bounces = new float[dimension];
//		for (int i = 0; i < dimension; i ++)
//			bounces[i] = he->bounces[2 * i];
//
//		//		delete he;
//		//--------------------------------------------------------
//
//		if (level==0)
//		{
//			//checked if this point is in the constrained space
//			if(MINDIST(bounces, _bounces, dimension) > FLOATZERO)
//			{
//				delete bounces;
//				continue;
//			}
//
//			if (!sky_dom(dimension, _rslt, _rsltcnt, bounces))
//			{
//				memcpy(&_rslt[dimension * _rsltcnt], bounces, sizeof(float) * dimension);
//				_rsltcnt ++;
//			}
//		}
//		else
//		{
//			if (!sky_dom(dimension, _rslt, _rsltcnt, bounces))
//			{
//				RTNode *child = new RTNode(this, son);
//
//				for (i = 0; i < child->num_entries; i ++)
//				{
//					float *b = new float[dimension];
//					for (j = 0; j < dimension; j ++) 
//						b[j] = child->entries[i].bounces[2 * j];
//
//					//checked if this entry is in the constrained subspace
//					if(MbrMINDIST(child->entries[i].bounces, _bounces, dimension) > FLOATZERO)
//					{
//						delete b;
//						continue;
//					}
//
//					if (!sky_dom(dimension, _rslt, _rsltcnt, b))
//					{
//						//HeapEntry *he = new HeapEntry();
//						he->dim = dimension;
//						he->son1 = child->entries[i].son;
//						he->level = child->level;
//						he->key = sky_mindist(dimension, b);
//
//						//he->bounces = new float[2 * dimension];
//						for (j = 0; j < 2*dimension; j++)
//							he->bounces[j] = child->entries[i].bounces[j]; 
//
//						_hp->insert(he);
//						//delete he;
//						//------------------------------------------------
//					}
//					delete [] b;
//				}
//				delete child;
//			}
//		}
//		delete [] bounces;
//	}
//	return;	
//}
//
////////////////////////////////////////////////////////////////////////////
//// functions added for reverse skyline
////////////////////////////////////////////////////////////////////////////
//
////=============================================================
//// BBRS algorithm for reverse skyline queries
//
//// _hp:      heap to be used
//// _rslt:    reverse skyline points list
//// _rsltcnt: number of reverse skyline points  
//// _qry_pt:  query point
//
//// coded by Junfeng Hu, on Nov. 3, 2007
////=============================================================
//void RTree::bbrs(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	Heap *brq_hp = new Heap();	// heap for bool range query
//	brq_hp->init(dimension, 1000);
//
//	float *global_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of global skyline points
//	int global_skyline_cnt = 0; // number of global skyline points
//
//	// insert all entries of the root in the Heap _hp sorted by distance from query point
//	_hp->used = 0;
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; // not important, but make sure it's greater than 0
//	he->bounces = new float[2*dimension];
//	for (i = 0; i < 2 * dimension; i++)
//		he->bounces[i] = 0.0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove top entry
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		// if the top entry is globally dominated by some global skyline points found yet
//		// then we discard it.
//		if (global_dom(dimension, global_skyline_pts, global_skyline_cnt, he->bounces, _qry_pt))
//			continue;
//
//		if (level == 0) 
//		{
//			// a new global skyline point is found
//
//			// (1) insert it into global skyline list
//			float *gpt = new float[dimension]; // the new global skyline point
//			for (i = 0; i < dimension; i++)
//				gpt[i] = he->bounces[2*i];
//
//			memcpy(&global_skyline_pts[dimension*global_skyline_cnt], gpt, sizeof(float)*dimension);
//			global_skyline_cnt++;
//
//			// (2) execute a window query based on gpt and the query point _qry_pt
//			float *range = new float[2*dimension];
//			for (i = 0; i < dimension; i++)
//			{
//				if (gpt[i] < _qry_pt[i])
//				{
//					range[2*i] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//					range[2*i + 1] = _qry_pt[i];
//				}
//				else
//				{
//					range[2*i] = _qry_pt[i];
//					range[2*i + 1] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//				}
//			}
//
//			brq_hp->used = 0;
//			if (!bool_range_query(brq_hp, range, gpt)) // if bool range query return false
//			{
//				// a new reverse skyline point is found
//				memcpy(&_rslt[dimension*_rsltcnt], gpt, sizeof(float)*dimension);
//				_rsltcnt++;
//			}
//
//			delete[] range;
//			delete[] gpt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				if (!global_dom(dimension, global_skyline_pts, global_skyline_cnt
//					, child->entries[i].bounces, _qry_pt))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//					memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}		
//	}
//
//	delete brq_hp;
//	delete he;
//	delete[] global_skyline_pts;
//}
//
////=============================================================
//// supporting function used to test whether a mbr (including a
//// point) is globally dominated by a list of points
//
//// _dim:        number of dimensions
//// _skyline:    list of points, which usually are skyline points
//// _skylinecnt: number of points in _skyline
//// _mbr:        the mbr to be tested
//// _qry_pt:	    query point of reverse skyline
//
//// coded by Junfeng Hu, on Nov. 3, 2007
////=============================================================
//bool RTree::global_dom(int _dim, float *_skyline, int _skylinecnt, float *_mbr, float *_qry_pt)
//{
//	int i, j;
//	bool ret = false;
//
//	for (i = 0; i < _skylinecnt; i++)
//	{
//		bool flag = true;
//		for (j = 0; j < _dim; j++)
//		{
//			if ( (_skyline[i*_dim + j] - _qry_pt[j])*(_mbr[2*j] - _qry_pt[j]) < 0 
//				|| (_skyline[i*_dim + j] - _qry_pt[j])*(_mbr[2*j + 1] - _qry_pt[j]) < 0 )
//			{
//				flag = false;
//				break;
//			}
//		}
//
//		if (flag)
//		{
//			// this point, with its coordinates transformed, is to be checked if it is globally dominated 
//			float *pt = new float[_dim]; 
//			for (j = 0; j < _dim; j++)
//			{
//				if ( (_mbr[2*j] - _qry_pt[j]) < (_mbr[2*j + 1] - _qry_pt[j]) )
//					pt[j] = _mbr[2*j] - _qry_pt[j];
//				else
//					pt[j] = _mbr[2*j + 1] - _qry_pt[j];
//
//				if (pt[j] < 0) pt[j] = (-1)*pt[j]; // make sure it's positive
//			}
//
//			// the global skyline point whose coordinates are transformed
//			float *gpt = new float[_dim]; 
//			for (j = 0; j < _dim; j++)
//			{
//				gpt[j] = _skyline[i*_dim + j] - _qry_pt[j];
//				if (gpt[j] < 0) gpt[j] = (-1)*gpt[j]; // make sure it's positive
//			}
//
//			if (dominate(_dim, gpt, pt))
//			{
//				ret = true; 
//				break;
//			}
//
//			delete[] gpt;
//			delete[] pt;
//		}
//	}
//
//	return ret;
//}
//
////=============================================================
//// Algorithm for boolean range query
//
//// _hp:     heap to be used
//// _mbr:    range to be tested
//// _gpt:	center of _mbr, which should be ignored in bbrs alogrithm
//
//// coded by Junfeng Hu, on Nov. 9, 2007
////=============================================================
//bool RTree::bool_range_query(Heap *_hp, float *_mbr, float *_gpt)
//{
//	int i, j;
//
//	// construct the root node 
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //make sure it is greater than 0
//	he->bounces = new float[2*dimension]; // it will be released in the destructor of he	
//	for (i = 0; i < 2*dimension; i++)
//		he->bounces[i] = 0.0;
//
//	// insert the root node into the heap
//	_hp->used = 0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove the top heap
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		if (level == 0) // it is a leaf node
//		{
//			// pruning the global skyline itself
//			bool is_center = true;
//			for (j = 0; j < dimension; j++)
//			{
//				if ((he->bounces[2*j] - _gpt[j]) > FLOATZERO 
//					|| (_gpt[j] - he->bounces[2*j]) > FLOATZERO) 
//				{
//					is_center = false;
//					break;
//				}
//			}
//			if (is_center)
//				continue;
//
//			float *pt = new float[dimension]; // the point in this leaf node
//			for (i = 0; i < dimension; i++)
//				pt[i] = he->bounces[2*i];
//
//			// check if this point is in the range _mbr
//			if (MINDIST(pt, _mbr, dimension) < FLOATZERO)
//			{
//				// then return true, since a point is found in the range 
//				delete[] pt;
//				delete he;
//				return true;
//			}
//			delete[] pt;
//		}
//		else // it is an internal node
//		{
//			float *bottom_left = new float[dimension]; // the bottom-left corner of mbr
//			float *top_right = new float[dimension];   // the top-right corner of mbr
//
//			// check if at least one edge of bounces is in the range
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				// pruning all the entries that do not intersect with our range
//				if (MbrMINDIST(child->entries[i].bounces, _mbr, dimension) > FLOATZERO)
//					continue;
//
//				// pruning the global skyline itself
//				if (child->level == 0)
//				{
//					bool is_center = true;
//					for (j = 0; j < dimension; j++)
//					{
//						if ((child->entries[i].bounces[2*j] - _gpt[j]) > FLOATZERO 
//							|| (_gpt[j] - child->entries[i].bounces[2*j]) > FLOATZERO) 
//						{
//							is_center = false;
//							break;
//						}
//					}
//					if (is_center)
//						continue;
//				}
//
//				bool flag = false; // this flag is true if _gpt is at one edge of current mbr
//				for (j = 0; j < dimension; j++)
//				{
//					bottom_left[j] = child->entries[i].bounces[2*j];
//					top_right[j] = child->entries[i].bounces[2*j + 1];
//					if ( ((bottom_left[j] - _gpt[j]) < FLOATZERO && (_gpt[j] - bottom_left[j]) < FLOATZERO) 
//						|| ((top_right[j] - _gpt[j]) < FLOATZERO && (_gpt[j] - top_right[j]) < FLOATZERO))
//					{
//						flag = true;
//						break;
//					}
//				}
//
//				bool rslt = false;
//				if (!flag)
//				{
//					// Both the bottom-left corner and the top-right corner is in the range
//					if (MINDIST(bottom_left, _mbr, dimension) < FLOATZERO 
//						&& MINDIST(top_right, _mbr, dimension) < FLOATZERO)
//					{
//						rslt = true; // There's at least one point in the range
//					}
//					// Only the bottom-left corner is in the range
//					else if (MINDIST(bottom_left, _mbr, dimension) < FLOATZERO)
//					{
//						float *pt = new float[dimension];
//						for (j = 0; j < dimension; j++)
//						{
//							memcpy(pt, top_right, sizeof(float) * dimension);
//							pt[j] = bottom_left[j];
//
//							if (MINDIST(pt, _mbr, dimension) < FLOATZERO)
//							{
//								rslt = true; // There's at least one point in the range
//								break;
//							}
//						}
//						delete[] pt;
//					}
//					// Only the top-right corner is in the range
//					else if (MINDIST(top_right, _mbr, dimension) < FLOATZERO)
//					{
//						float *pt = new float[dimension];
//						for (j = 0; j < dimension; j++)
//						{
//							memcpy(pt, bottom_left, sizeof(float) * dimension);
//							pt[j] = top_right[j];
//
//							if (MINDIST(pt, _mbr, dimension) < FLOATZERO)
//							{
//								rslt = true; // There's at least one point in the range
//								break;
//							}
//						}
//						delete[] pt;
//					}
//				}
//
//
//				if (rslt) // There's at least one point in the range
//				{
//					delete[] bottom_left;
//					delete[] top_right;
//					delete he;
//					return true;
//				}
//				else
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = 0; // it dosen't matter
//					for (j = 0; j < 2*dimension; j++)
//						he->bounces[j] = child->entries[i].bounces[j];
//
//					_hp->insert(he);
//				}
//			}
//
//			delete[] bottom_left;
//			delete[] top_right;
//		}
//	}
//
//	delete he;
//	return false;
//}
//
////=============================================================
//// Algorithm for dynamic skyline query
//
//// _hp:			heap to be used
//// _rslt:		list of dynamic skyline points
//// _rsltcnt:	number of dynamic skyline points
//// _qry_pt:		query point
//
//// coded by Junfeng Hu, on Nov. 14, 2007
////=============================================================
//void RTree::dynamic_bbs(Heap *_hp, float *_rslt, int &_rsltcnt, float *_qry_pt)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	//first insert the root into the heap-------------------------
//	_hp->used=0;
//	HeapEntry *he=new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //this is not important but make sure it is greather than 0
//	he->bounces = new float[2 * dimension];
//	for (j = 0; j < dimension; j ++)
//	{
//		he->bounces[2*j] = 0.0;
//		he->bounces[2*j+1] = 10000.0;
//	}
//	_hp->insert(he);
//
//	while(_hp->used > 0)
//	{
//		//remove the top heap-------------------------------------
//
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		if (dynamic_dom(dimension, _rslt, _rsltcnt, he->bounces, _qry_pt))
//			continue;
//
//		if (level==0)
//		{
//			float *pt = new float[dimension];
//			for (i = 0; i < dimension; i++)
//				pt[i] = he->bounces[2*i];
//
//			bool is_querypoint = true;
//			for (i = 0; i < dimension; i++)
//			{
//				if ((_qry_pt[i] - pt[i]) > FLOATZERO
//					|| (pt[i] - _qry_pt[i]) > FLOATZERO)
//				{
//					is_querypoint = false;
//					break;
//				}
//			}
//			if (!is_querypoint) // not query point
//			{
//				memcpy(&_rslt[dimension * _rsltcnt], pt, sizeof(float) * dimension);
//				_rsltcnt ++;
//			}
//
//			delete[] pt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//
//			for (i = 0; i < child->num_entries; i ++)
//			{
//				if (!dynamic_dom(dimension, _rslt, _rsltcnt, child->entries[i].bounces, _qry_pt))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//
//					he->key = dynamic_mindist(dimension, child->entries[i].bounces, _qry_pt);
//
//					for (j = 0; j < 2*dimension; j++)
//						he->bounces[j] = child->entries[i].bounces[j]; 
//
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}
//	}
//	delete he;
//
//	/*
//	int i, j;
//
//	// (1) transforming space
//	ofstream outfile("transformed_space.tmp");
//	if (!outfile)
//	{
//	cerr << "creating transformed_space.tmp failed!" << endl;
//	exit(1);
//	}
//	int data_count = 0;
//
//	HeapEntry *he = new HeapEntry();
//	he->init_HeapEntry(dimension); // it initials dim and bounces
//	he->son1 = root;
//	he->level = 1;
//	he->key = 0;
//	for (i = 0; i < 2*dimension; i++)
//	he->bounces[i] = 0.0;
//
//	_hp->used = 0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//	_hp->remove(he);
//	int son = he->son1;
//	int level = he->level;
//
//	if (level == 0) // it is a leaf node
//	{
//	float *pt = new float[dimension];
//	for (i = 0; i < dimension; i++)
//	{
//	pt[i] = he->bounces[2*i] - _qry_pt[i];
//	if (pt[i] < 0) pt[i] = (-1)*pt[i]; // ensure that it's positive
//	}
//
//	outfile << data_count++; 
//	for (i = 0; i < dimension; i++)
//	outfile << " " << fixed << setprecision(6) << pt[i];
//	outfile << "\n";
//
//	delete[] pt;
//	}
//	else // it is an internal node
//	{
//	RTNode *child = new RTNode(this, son);
//
//	for (i = 0; i < child->num_entries; i++)
//	{
//	he->dim = dimension;
//	he->son1 = child->entries[i].son;
//	he->level = child->level;
//	he->key = 0;
//	memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//	_hp->insert(he);
//	}
//
//	delete child;
//	}
//	}
//
//	outfile.close();
//	delete he;
//	*/
//
//	// (2) building a new R-Tree
//
//	// (3) using bbs to compute skyline based on the new R-Tree
//}
//
//float dynamic_mindist(int _dim, float *_bounces, float *_qry_pt)
//{
//	float ret = 0.0;
//	int i;
//	for (i = 0; i < _dim; i++)
//	{
//		float d1 = _bounces[2*i] - _qry_pt[i];		
//		float d2 = _bounces[2*i+1] - _qry_pt[i];
//		if (d1 * d2 < 0) // since _bounces[2*i] >= _bounces[2*i+1], here d2 > 0 and d1 < 0
//			ret += 0.0;
//		else // d1 >= 0 and d2 >= 0, or d1 <= 0 and d2 <= 0
//		{
//			if (d1 < 0) d1 *= (-1);
//			if (d2 < 0) d2 *= (-1);
//			ret += d1 < d2 ? d1 : d2; // the smaller of d1 or d2
//		}
//	}
//	return ret;
//}
//
//bool dynamic_dom(int _dim, float *_rslt, int _rsltcnt, float *_bounces, float *_qry_pt)
//{
//	int i, j;
//
//	float *tfpt = new float[_dim]; // lowest corner of mbr with coordinates transformed
//	for (i = 0; i < _dim; i++)
//	{
//		float d1 = _bounces[2*i] - _qry_pt[i];		
//		float d2 = _bounces[2*i+1] - _qry_pt[i];
//		if (d1 * d2 < 0) // since _bounces[2*i] >= _bounces[2*i+1], here d2 > 0 and d1 < 0
//			tfpt[i] = 0.0;
//		else // d1 >= 0 and d2 >= 0, or d1 <= 0 and d2 <= 0
//		{
//			if (d1 < 0) d1 *= (-1);
//			if (d2 < 0) d2 *= (-1);
//			tfpt[i] = d1 < d2 ? d1 : d2; // the smaller of d1 or d2
//		}
//	}
//
//	bool ret = false;
//
//	float *s_pt = new float[_dim]; // one skyline point
//	for (i = 0; i < _rsltcnt; i++)
//	{
//		for (j = 0; j < _dim; j++)
//		{
//			s_pt[j] = _rslt[i*_dim + j] - _qry_pt[j];
//			if (s_pt[j] < 0) s_pt[j] *= (-1); // ensure is is positive
//		}
//		if (dominate(_dim, s_pt, tfpt))
//		{
//			ret = true;
//			break;
//		}
//	}
//
//	delete[] s_pt;
//	delete[] tfpt;
//
//	return ret;
//}
//
//void RTree::reverse_skyline(Heap *_hp, float *_rslt, int &_rsltcnt, float *_qry_pt)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	float *dynamic_skyline = new float[MAX_SKYLINE_RESULT]; // set of dynamic skyline points
//	Heap *dynamic_hp = new Heap();	// heap for dynamic skyline query
//	dynamic_hp->init(dimension, 1000);
//
//	float *global_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of global skyline points
//	int global_skyline_cnt = 0; // number of global skyline points
//
//	// insert all entries of the root in the Heap _hp sorted by distance from query point
//	_hp->used = 0;
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; // not important, but make sure it's greater than 0
//	he->bounces = new float[2*dimension];
//	for (i = 0; i < 2 * dimension; i++)
//		he->bounces[i] = 0.0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove top entry
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		// if the top entry is globally dominated by some global skyline points found yet
//		// then we discard it.
//		if (global_dom(dimension, global_skyline_pts, global_skyline_cnt, he->bounces, _qry_pt))
//			continue;
//
//		if (level == 0) 
//		{
//			// a new global skyline point is found
//
//			// (1) insert it into global skyline list
//			float *gpt = new float[dimension]; // the new global skyline point
//			for (i = 0; i < dimension; i++)
//				gpt[i] = he->bounces[2*i];
//
//			memcpy(&global_skyline_pts[dimension*global_skyline_cnt], gpt, sizeof(float)*dimension);
//			global_skyline_cnt++;
//
//			int dynamic_skyline_cnt = 0; // number of dynamic skyline points
//			// compute the dynamic skyline, as gpt is the query point
//			dynamic_hp->used = 0;
//			dynamic_bbs(dynamic_hp, dynamic_skyline, dynamic_skyline_cnt, gpt);
//
//			float *b = new float[2*dimension];
//			for (i = 0; i < dimension; i++)
//			{
//				b[2*i] = _qry_pt[i];
//				b[2*i + 1] = _qry_pt[i];
//			}
//			if (!dynamic_dom(dimension, dynamic_skyline, dynamic_skyline_cnt, b, gpt)) 
//			{
//				// a new reverse skyline is found
//				memcpy(&_rslt[dimension*_rsltcnt], gpt, sizeof(float)*dimension);
//				_rsltcnt++;
//			}
//			delete[] b;
//
//			delete[] gpt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				if (!global_dom(dimension, global_skyline_pts, global_skyline_cnt
//					, child->entries[i].bounces, _qry_pt))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//					memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}		
//	}
//
//	delete dynamic_hp;
//	delete he;
//	delete[] global_skyline_pts;	
//	delete[] dynamic_skyline;
//}
//
////=============================================================
//// RSSA algorithm for reverse skyline queries
//
//// _hp:      heap to be used
//// _rslt:    reverse skyline points list
//// _rsltcnt: number of reverse skyline points  
//// _qry_pt:  query point
//
//// coded by Junfeng Hu, on Nov. 27, 2007
////=============================================================
//void RTree::rssa(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	Heap *brq_hp = new Heap();	// heap for bool range query
//	brq_hp->init(dimension, 1000);
//
//	float *global_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of global skyline points
//	int global_skyline_cnt = 0; // number of global skyline points
//
//	// insert all entries of the root in the Heap _hp sorted by distance from query point
//	_hp->used = 0;
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; // not important, but make sure it's greater than 0
//	he->bounces = new float[2*dimension];
//	for (i = 0; i < 2 * dimension; i++)
//		he->bounces[i] = 0.0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove top entry
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		// if the top entry is globally dominated by some global skyline points found yet
//		// then we discard it.
//		if (global_dom(dimension, global_skyline_pts, global_skyline_cnt, he->bounces, _qry_pt))
//			continue;
//
//		if (level == 0) 
//		{
//			// a new global skyline point is found
//
//			// (1) insert it into global skyline list
//			float *gpt = new float[dimension]; // the new global skyline point
//			for (i = 0; i < dimension; i++)
//				gpt[i] = he->bounces[2*i];
//
//			memcpy(&global_skyline_pts[dimension*global_skyline_cnt], gpt, sizeof(float)*dimension);
//			global_skyline_cnt++;
//
//			// (2) filter step
//			if (DDR(gpt, _qry_pt))
//			{
//				//TODO...
//			}
//			else if (DADR(gpt, _qry_pt))
//			{
//				//TODO...
//			}
//			else
//			{
//				// (3) refinement step
//				//TODO...				
//			}
//
//			delete[] gpt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				if (!global_dom(dimension, global_skyline_pts, global_skyline_cnt
//					, child->entries[i].bounces, _qry_pt))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//					memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}		
//	}
//
//	delete brq_hp;
//	delete he;
//	delete[] global_skyline_pts;
//}
//
//bool RTree::DDR(float* _p, float *_q)
//{
//	// TODO...
//	return true;
//}
//
//bool RTree::DADR(float* _p, float *_q)
//{
//	// TODO...
//	return true;
//}
//
////////////////////////////////////////////////////////////////////////////
//// functions added for variations of reverse skyline
////////////////////////////////////////////////////////////////////////////
//
////=============================================================
//// BBRS algorithm for reverse skyline queries in a constrained 
//// space
//
//// _hp:      heap to be used
//// _rslt:    reverse skyline points list
//// _rsltcnt: number of reverse skyline points  
//// _qry_pt:  query point
//// _bounces: constrained region
//
//// coded by Junfeng Hu, on Nov. 28, 2007
////=============================================================
//void RTree::bbrs_constrained(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, float *_bounces)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	Heap *brq_hp = new Heap();	// heap for bool range query
//	brq_hp->init(dimension, 1000);
//
//	float *global_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of global skyline points
//	int global_skyline_cnt = 0; // number of global skyline points
//
//	// insert all entries of the root in the Heap _hp sorted by distance from query point
//	_hp->used = 0;
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; // not important, but make sure it's greater than 0
//	he->bounces = new float[2*dimension];
//	for (i = 0; i < 2 * dimension; i++)
//		he->bounces[i] = 0.0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove top entry
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		// if the top entry is globally dominated by some global skyline points found yet
//		// then we discard it.
//		if (global_dom(dimension, global_skyline_pts, global_skyline_cnt, he->bounces, _qry_pt))
//			continue;
//
//		if (level == 0) 
//		{
//			// a new global skyline point is found
//
//			// (1) insert it into global skyline list
//			float *gpt = new float[dimension]; // the new global skyline point
//			for (i = 0; i < dimension; i++)
//				gpt[i] = he->bounces[2*i];
//
//			// check if this point is in the constrained space
//			if(MINDIST(gpt, _bounces, dimension) > FLOATZERO)
//			{
//				delete gpt;
//				continue;
//			}
//
//			memcpy(&global_skyline_pts[dimension*global_skyline_cnt], gpt, sizeof(float)*dimension);
//			global_skyline_cnt++;
//
//			// (2) execute a window query based on gpt and the query point _qry_pt
//			float *range = new float[2*dimension];
//			for (i = 0; i < dimension; i++)
//			{
//				if (gpt[i] < _qry_pt[i])
//				{
//					range[2*i] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//					range[2*i + 1] = _qry_pt[i];
//				}
//				else
//				{
//					range[2*i] = _qry_pt[i];
//					range[2*i + 1] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//				}
//			}
//
//			brq_hp->used = 0;
//			if (!bool_range_query(brq_hp, range, gpt)) // if bool range query return false
//			{
//				// a new reverse skyline point is found
//				memcpy(&_rslt[dimension*_rsltcnt], gpt, sizeof(float)*dimension);
//				_rsltcnt++;
//			}
//
//			delete[] range;
//			delete[] gpt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				//checked if this entry is in the constrained subspace
//				if(MbrMINDIST(child->entries[i].bounces, _bounces, dimension) > FLOATZERO)
//				{
//					continue;
//				}
//
//				if (!global_dom(dimension, global_skyline_pts, global_skyline_cnt
//					, child->entries[i].bounces, _qry_pt))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//					memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}		
//	}
//
//	delete brq_hp;
//	delete he;
//	delete[] global_skyline_pts;
//}
//
//void RTree::rssa_constrained(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, float *_bounces)
//{
//	// TODO...
//}
//
//
//bool cmp_value_p(const struct value_skyline_point &x, const struct value_skyline_point &y)
//{
//	return x.num > y.num;
//}
//
//void RTree::bbrs_enum(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k)
//{
//	int i, j;
//	_rsltcnt = 0;
//	float *reverse_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of reverse skyline points
//	int reverse_skyline_cnt = 0;
//
//	Heap *brq_hp = new Heap();	// heap for bool range query
//	brq_hp->init(dimension, 1000);
//
//	float *global_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of global skyline points
//	int global_skyline_cnt = 0; // number of global skyline points
//
//	// insert all entries of the root in the Heap _hp sorted by distance from query point
//	_hp->used = 0;
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; // not important, but make sure it's greater than 0
//	he->bounces = new float[2*dimension];
//	for (i = 0; i < 2 * dimension; i++)
//		he->bounces[i] = 0.0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove top entry
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		// if the top entry is globally dominated by some global skyline points found yet
//		// then we discard it.
//		if (global_dom(dimension, global_skyline_pts, global_skyline_cnt, he->bounces, _qry_pt))
//			continue;
//
//		if (level == 0) 
//		{
//			// a new global skyline point is found
//
//			// (1) insert it into global skyline list
//			float *gpt = new float[dimension]; // the new global skyline point
//			for (i = 0; i < dimension; i++)
//				gpt[i] = he->bounces[2*i];
//
//			memcpy(&global_skyline_pts[dimension*global_skyline_cnt], gpt, sizeof(float)*dimension);
//			global_skyline_cnt++;
//
//			// (2) execute a window query based on gpt and the query point _qry_pt
//			float *range = new float[2*dimension];
//			for (i = 0; i < dimension; i++)
//			{
//				if (gpt[i] < _qry_pt[i])
//				{
//					range[2*i] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//					range[2*i + 1] = _qry_pt[i];
//				}
//				else
//				{
//					range[2*i] = _qry_pt[i];
//					range[2*i + 1] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//				}
//			}
//
//			brq_hp->used = 0;
//			if (!bool_range_query(brq_hp, range, gpt)) // if bool range query return false
//			{
//				// a new reverse skyline point is found
//				memcpy(&reverse_skyline_pts[dimension*reverse_skyline_cnt], gpt, sizeof(float)*dimension);
//				reverse_skyline_cnt++;
//			}
//
//			delete[] range;
//			delete[] gpt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				if (!global_dom(dimension, global_skyline_pts, global_skyline_cnt
//					, child->entries[i].bounces, _qry_pt))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//					memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}		
//	}
//
//	vector<struct value_skyline_point> enum_reverse_skyline;
//	float *range = new float[2*dimension];
//	float *reverse_spt = new float[dimension];
//	for (i = 0; i < reverse_skyline_cnt; i++)
//	{
//		struct value_skyline_point p;
//		p.i = i;
//		brq_hp->used = 0;
//		for (j = 0; j < dimension; j++)
//		{
//			range[2*j] = reverse_skyline_pts[i*dimension+j];
//			range[2*j+1] = 10000.0;
//
//			reverse_spt[2*j] = reverse_skyline_pts[i*dimension+j];
//		}
//		p.num = enum_range_query(brq_hp, range, reverse_spt);
//		enum_reverse_skyline.push_back(p);
//	}
//	delete[] range;
//	delete[] reverse_spt;
//
//	// sort skyline points using the compare function cmp_value_p.
//	sort(enum_reverse_skyline.begin(), enum_reverse_skyline.end(), cmp_value_p);
//
//	for (i = 0; i < _k; i++)
//		for (j = 0; j < dimension; j++)
//			_rslt[dimension*i+j] = reverse_skyline_pts[dimension*(enum_reverse_skyline[i].i)+j];
//	_rsltcnt = reverse_skyline_cnt;
//
//	delete brq_hp;
//	delete he;
//	delete[] global_skyline_pts;
//	delete[] reverse_skyline_pts;
//}
//
//void RTree::rssa_enum(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k)
//{
//	// TODO...
//}
//
////=============================================================
//// Algorithm for enumeration range query, which return number 
//// of points in the range.
//
//// _hp:     heap to be used
//// _mbr:    range to be tested
//// _gpt:	center of _mbr, which should be ignored in bbrs alogrithm
//
//// coded by Junfeng Hu, on Nov. 29, 2007
////=============================================================
//int RTree::enum_range_query(Heap *_hp, float *_mbr, float *_gpt)
//{
//	int i, j;
//	int num = 0;
//
//	// construct the root node 
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; //make sure it is greater than 0
//	he->bounces = new float[2*dimension]; // it will be released in the destructor of he	
//	for (i = 0; i < 2*dimension; i++)
//		he->bounces[i] = 0.0;
//
//	// insert the root node into the heap
//	_hp->used = 0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove the top heap
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		if (level == 0) // it is a leaf node
//		{
//			// pruning the global skyline itself
//			bool is_center = true;
//			for (j = 0; j < dimension; j++)
//			{
//				if ((he->bounces[2*j] - _gpt[j]) > FLOATZERO 
//					|| (_gpt[j] - he->bounces[2*j]) > FLOATZERO) 
//				{
//					is_center = false;
//					break;
//				}
//			}
//			if (is_center)
//				continue;
//
//			float *pt = new float[dimension]; // the point in this leaf node
//			for (i = 0; i < dimension; i++)
//				pt[i] = he->bounces[2*i];
//
//			// check if this point is in the range _mbr
//			if (MINDIST(pt, _mbr, dimension) < FLOATZERO)
//			{
//				// then return true, since a point is found in the range 
//				// 				delete[] pt;
//				// 				delete he;
//				// 				return true;
//				num++;
//			}
//			delete[] pt;
//		}
//		else // it is an internal node
//		{
//			RTNode *child = new RTNode(this, son);
//
//			float *bounces = new float[dimension]; // the point in this leaf node
//			for (i = 0; i < 2*dimension; i++)
//				bounces[i] = he->bounces[i];
//
//			// check if the internal node is completely in _mbr
//			bool flag = true;
//			for (i = 0; i < dimension; i++)
//			{
//				if (bounces[2*i] < _mbr[2*i] || bounces[2*i+1] > _mbr[2*i+1])
//				{
//					flag = false;
//					break;
//				}
//			}
//			if (flag)
//			{
//				num += child->num_entries;
//				continue;
//			}
//
//			// else then the internal node intersect with _mbr
//			for (i = 0; i < child->num_entries; i++)
//			{
//				// pruning all the entries that do not intersect with our range
//				if (MbrMINDIST(child->entries[i].bounces, _mbr, dimension) > FLOATZERO)
//					continue;
//
//				// pruning the global skyline itself
//				if (child->level == 0)
//				{
//					bool is_center = true;
//					for (j = 0; j < dimension; j++)
//					{
//						if ((child->entries[i].bounces[2*j] - _gpt[j]) > FLOATZERO 
//							|| (_gpt[j] - child->entries[i].bounces[2*j]) > FLOATZERO) 
//						{
//							is_center = false;
//							break;
//						}
//					}
//					if (is_center)
//						continue;
//				}
//
//				he->dim = dimension;
//				he->son1 = child->entries[i].son;
//				he->level = child->level;
//				he->key = 0; // it dosen't matter
//				for (j = 0; j < 2*dimension; j++)
//					he->bounces[j] = child->entries[i].bounces[j];
//
//				_hp->insert(he);
//			}
//		}
//	}
//
//	delete he;
//	return num;
//}
//
//void RTree::bbrs_ranked(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k, float (*_mindist)(int , float*, float* ))
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	Heap *brq_hp = new Heap();	// heap for bool range query
//	brq_hp->init(dimension, 1000);
//
//	float *global_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of global skyline points
//	int global_skyline_cnt = 0; // number of global skyline points
//
//	// insert all entries of the root in the Heap _hp sorted by distance from query point
//	_hp->used = 0;
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; // not important, but make sure it's greater than 0
//	he->bounces = new float[2*dimension];
//	for (i = 0; i < 2 * dimension; i++)
//		he->bounces[i] = 0.0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove top entry
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		// if the top entry is globally dominated by some global skyline points found yet
//		// then we discard it.
//		if (global_dom(dimension, global_skyline_pts, global_skyline_cnt, he->bounces, _qry_pt))
//			continue;
//
//		if (level == 0) 
//		{
//			// a new global skyline point is found
//
//			// (1) insert it into global skyline list
//			float *gpt = new float[dimension]; // the new global skyline point
//			for (i = 0; i < dimension; i++)
//				gpt[i] = he->bounces[2*i];
//
//			memcpy(&global_skyline_pts[dimension*global_skyline_cnt], gpt, sizeof(float)*dimension);
//			global_skyline_cnt++;
//
//			// (2) execute a window query based on gpt and the query point _qry_pt
//			float *range = new float[2*dimension];
//			for (i = 0; i < dimension; i++)
//			{
//				if (gpt[i] < _qry_pt[i])
//				{
//					range[2*i] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//					range[2*i + 1] = _qry_pt[i];
//				}
//				else
//				{
//					range[2*i] = _qry_pt[i];
//					range[2*i + 1] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//				}
//			}
//
//			brq_hp->used = 0;
//			if (!bool_range_query(brq_hp, range, gpt)) // if bool range query return false
//			{
//				// a new reverse skyline point is found
//				memcpy(&_rslt[dimension*_rsltcnt], gpt, sizeof(float)*dimension);
//				_rsltcnt++;
//			}
//
//			delete[] range;
//			delete[] gpt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				if (!global_dom(dimension, global_skyline_pts, global_skyline_cnt
//					, child->entries[i].bounces, _qry_pt))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = _mindist(dimension, child->entries[i].bounces, _qry_pt);
//					memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}		
//	}
//
//	delete brq_hp;
//	delete he;
//	delete[] global_skyline_pts;
//}
//
//void RTree::rssa_ranked(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k, float (*_mindist)(int, float*, float*))
//{
//	// TODO...
//}
//
//bool RTree::global_dom_skyband(int _dim, float *_skyline, int _skylinecnt, float *_mbr, float *_qry_pt, int _k)
//{
//	int i, j;
//	bool ret = true;
//	int num = 0;
//
//	for (i = 0; i < _skylinecnt; i++)
//	{
//		bool flag = true;
//		for (j = 0; j < _dim; j++)
//		{
//			if ( (_skyline[i*_dim + j] - _qry_pt[j])*(_mbr[2*j] - _qry_pt[j]) < 0 
//				|| (_skyline[i*_dim + j] - _qry_pt[j])*(_mbr[2*j + 1] - _qry_pt[j]) < 0 )
//			{
//				flag = false;
//				break;
//			}
//		}
//
//		if (flag)
//		{
//			// this point, with its coordinates transformed, is to be checked if it is globally dominated 
//			float *pt = new float[_dim]; 
//			for (j = 0; j < _dim; j++)
//			{
//				if ( (_mbr[2*j] - _qry_pt[j]) < (_mbr[2*j + 1] - _qry_pt[j]) )
//					pt[j] = _mbr[2*j] - _qry_pt[j];
//				else
//					pt[j] = _mbr[2*j + 1] - _qry_pt[j];
//
//				if (pt[j] < 0) pt[j] = (-1)*pt[j]; // make sure it's positive
//			}
//
//			// the global skyline point whose coordinates are transformed
//			float *gpt = new float[_dim]; 
//			for (j = 0; j < _dim; j++)
//			{
//				gpt[j] = _skyline[i*_dim + j] - _qry_pt[j];
//				if (gpt[j] < 0) gpt[j] = (-1)*gpt[j]; // make sure it's positive
//			}
//
//			if (dominate(_dim, gpt, pt))
//			{
//				num++;
//				if (num > _k)
//				{
//					ret = false; 
//					break;
//				}
//			}
//
//			delete[] gpt;
//			delete[] pt;
//		}
//	}
//
//	return ret;
//}
//
//void RTree::bbrs_skyband(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k)
//{
//	int i, j;
//	_rsltcnt = 0;
//
//	Heap *brq_hp = new Heap();	// heap for bool range query
//	brq_hp->init(dimension, 1000);
//
//	float *global_skyline_pts = new float[MAX_SKYLINE_RESULT]; // set of global skyline points
//	int global_skyline_cnt = 0; // number of global skyline points
//
//	// insert all entries of the root in the Heap _hp sorted by distance from query point
//	_hp->used = 0;
//	HeapEntry *he = new HeapEntry();
//	he->dim = dimension;
//	he->son1 = root;
//	he->key = 0;
//	he->level = 1; // not important, but make sure it's greater than 0
//	he->bounces = new float[2*dimension];
//	for (i = 0; i < 2 * dimension; i++)
//		he->bounces[i] = 0.0;
//	_hp->insert(he);
//
//	while (_hp->used > 0)
//	{
//		// remove top entry
//		_hp->remove(he);
//		int son = he->son1;
//		int level = he->level;
//
//		// if the top entry is globally dominated by some global skyline points found yet
//		// then we discard it.
//		if (global_dom_skyband(dimension, global_skyline_pts, global_skyline_cnt, he->bounces, _qry_pt, _k))
//			continue;
//
//		if (level == 0) 
//		{
//			// a new global skyline point is found
//
//			// (1) insert it into global skyline list
//			float *gpt = new float[dimension]; // the new global skyline point
//			for (i = 0; i < dimension; i++)
//				gpt[i] = he->bounces[2*i];
//
//			memcpy(&global_skyline_pts[dimension*global_skyline_cnt], gpt, sizeof(float)*dimension);
//			global_skyline_cnt++;
//
//			// (2) execute a window query based on gpt and the query point _qry_pt
//			float *range = new float[2*dimension];
//			for (i = 0; i < dimension; i++)
//			{
//				if (gpt[i] < _qry_pt[i])
//				{
//					range[2*i] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//					range[2*i + 1] = _qry_pt[i];
//				}
//				else
//				{
//					range[2*i] = _qry_pt[i];
//					range[2*i + 1] = _qry_pt[i] + 2*(gpt[i] - _qry_pt[i]);
//				}
//			}
//
//			brq_hp->used = 0;
//			// gpt is pruned only if a windows query for this point returns more than _k points.
//			if (enum_range_query(brq_hp, range, gpt) <= _k) 
//			{
//				// a new reverse skyline point is found
//				memcpy(&_rslt[dimension*_rsltcnt], gpt, sizeof(float)*dimension);
//				_rsltcnt++;
//			}
//
//			delete[] range;
//			delete[] gpt;
//		}
//		else
//		{
//			RTNode *child = new RTNode(this, son);
//			for (i = 0; i < child->num_entries; i++)
//			{
//				if (!global_dom_skyband(dimension, global_skyline_pts, global_skyline_cnt
//					, child->entries[i].bounces, _qry_pt, _k))
//				{
//					he->dim = dimension;
//					he->son1 = child->entries[i].son;
//					he->level = child->level;
//					he->key = MINDIST(_qry_pt, child->entries[i].bounces, dimension);
//					memcpy(he->bounces, child->entries[i].bounces, sizeof(float)*2*dimension);
//					_hp->insert(he);
//				}
//			}
//			delete child;
//		}		
//	}
//
//	delete brq_hp;
//	delete he;
//	delete[] global_skyline_pts;	
//}
//
//void RTree::rssa_skyband(Heap *_hp, float *_rslt, int &_rsltcnt, float* _qry_pt, int _k)
//{
//	// TODO...
//}