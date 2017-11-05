/* entry.cpp
   implementation of class Entry */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "entry.h"
#include "rtnode.h"
#include "../linlist/linlist.h"
//------------------------------------------------------------
Entry::Entry()
  //this ructor does nothing.  remember you should call init_entry
  //to initialize if you use this ructor
{
	son_ptr = NULL;
	bounces = NULL;
}
//------------------------------------------------------------
Entry::Entry(int _dimension, RTree *rt)
{
    dimension = _dimension;
    my_tree = rt;
    bounces = new float[2*dimension];
    son_ptr = NULL;
    son = 0;
	level = 0;
	type=10;
}


//------------------------------------------------------------
Entry::~Entry()
{
    if (bounces)
		delete [] bounces;
    if (son_ptr != NULL)
	    delete son_ptr;
}
//------------------------------------------------------------
void Entry::del_son()
{
	if (son_ptr != NULL)
	{
		delete son_ptr;
		son_ptr = NULL;
	}
}	
//------------------------------------------------------------
Linkable* Entry::gen_Linkable()
{
	Linkable *new_link = new Linkable(dimension);
	new_link -> son = son;
	new_link->type=type;
	//memcpy(new_link -> bounces, bounces, 2 * dimension * sizeof(float));
	for (int i = 0; i < 2 * dimension; i ++)
		new_link -> bounces[i] = bounces[i];
	new_link -> level = level;

	return new_link;
}
//------------------------------------------------------------
int Entry::get_size()
{
    return 2 * dimension * sizeof(float) + sizeof(int)+sizeof(long);
	  //for bounces and son
}
//------------------------------------------------------------
RTNode* Entry::get_son()
{
    if (son_ptr == NULL)
	    son_ptr = new RTNode(my_tree, son);

    return son_ptr;
}
//------------------------------------------------------------
void Entry::init_entry(int _dimension, RTree *_rt)
{
	dimension = _dimension;
    my_tree = _rt;
    bounces = new float[2 * dimension];
    son_ptr = NULL;
    son = 0;
	level = 0;
	type=10;
}
//------------------------------------------------------------
void Entry::read_from_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(bounces, buffer, i);

    memcpy(&son, &buffer[i], sizeof(int));
    i += sizeof(int);
	memcpy(&type,&buffer[i],sizeof(long));
	i+=sizeof(long);
}
//------------------------------------------------------------
SECTION Entry::section(float *mbr)
{
    bool inside;
    bool overlap;

    overlap = TRUE;
    inside = TRUE;

    for (int i = 0; i < dimension; i++)
    {
		if (mbr[2 * i] > bounces[2 * i + 1] ||  mbr[2 * i + 1] < bounces[2 * i])
			overlap = FALSE;
		if (mbr[2 * i] < bounces[2 * i] ||
			mbr[2 * i + 1] > bounces[2 * i + 1])
			inside = FALSE;
    }
    if (inside)
		return INSIDE;
    else if (overlap)
		return OVERLAP;
    else
		return S_NONE;
}
//------------------------------------------------------------
void Entry::set_from_Linkable(Linkable *link)
{
	son = link -> son;
	dimension = link -> dimension;
	type=link->type;
	memcpy(bounces, link -> bounces, 2 * dimension * sizeof(float));
	level = link -> level;

	my_tree = NULL;
	son_ptr = NULL;
}
//------------------------------------------------------------
void Entry::write_to_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(buffer, bounces, i);

    memcpy(&buffer[i], &son, sizeof(int));
    i += sizeof(int);
	//memcpy()
	memcpy(&buffer[i],&type,sizeof(long));
	i+=sizeof(long);
}
//------------------------------------------------------------
bool Entry::operator == (Entry &_d)
  //this function compares two entries based on (1)son (2)dimension (3)extents
{
	if (son != _d.son) return false;
	if (type!=_d.type) return false;
	if (dimension != _d.dimension) return false;
	for (int i = 0; i < 2 * dimension; i++)
		if (fabs(bounces[i] - _d.bounces[i]) > FLOATZERO) return false;
	return true;
}
//------------------------------------------------------------
Entry& Entry::operator = (Entry &_d)
  //this function assigns all fieds of _d with the same values of this entry
{
    dimension = _d.dimension;
    son = _d.son;
	type=_d.type;
    son_ptr = _d.son_ptr;
    memcpy(bounces, _d.bounces, sizeof(float) * 2 * dimension);
    my_tree = _d.my_tree;
	level = _d.level;

    return *this;
}
//------------------------------------------------------------