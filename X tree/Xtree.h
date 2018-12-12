#ifndef XTREE_H
#define XTREE_H
#include <stack>
#include <ctime>
#include "node.h"
#define maxoverlap 100

class Xtree
{
public:
	unsigned int m;
	unsigned int M;
	Xtree();
	dimtype n;
	Xtree(unsigned int max_children, dimtype nn);
	void xtinsert(point d);
	bool onrec(point p, point gmin, point gmax);
	void multidivide(vector<node*> lower, int unsigned numth,int np, int unsigned inc, node*& l, node*& r, typecor& S, bool& w);
	typecor distance_points(point p1, point p2);
	typecor distance_rec_point(point d, point cmin, point cmax);
	void bigger_Rec(point Amin, point Amax, point Bmin, point Bmax, point& Cmin, point& Cmax);
	void adjusttree(node**& p, stack<node**>& path, point d);
	void splitnode(node**& p, stack<node**>& path);
	typecor margin(node* g1, node* g2);
	bool multiexecute_calcS(vector<node*> lower, int unsigned start, int unsigned end, node*& l, node*& r, typecor& S);
	bool overlapaxis(point p1min, point p1max, point p2min, point p2max, dimtype x, typecor& men, typecor& may);
	void neighbors(node* p, vector <point>& points, point center, typecor radio);
	typecor area_overlap(point p1min, point p1max, point p2min, point p2max);
	bool stardivide(node**& p, node*& l, node*& r,dimtype x);
	bool minimaldivide(node**& p, node*& l, node*& r, dimtype x);
	bool outlimits(point limsl, point limsr, node* c);
	typecor area(point cmin, point cmax);
	node* root;
	void chooseleaf(point d, node**& p, stack<node**>& path);
	virtual ~Xtree();

protected:

private:
};

#endif // XTREE_H
