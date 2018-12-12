#include "node.h"

int unsigned ith;

bool byI(point a, point b)
{
	return a[ith]<b[ith];
}

node::node()
{
	data = 0;
	supernode = 0;
	splitaxis = -1;
}

node::node(vector <node*> c)
{
	data = 0;
	supernode = 0;
	children = c;
	splitaxis = -1;
}

node::node(point d)
{
	data = new point(d);
	supernode = 0;
	splitaxis = -1;
}

node::node(point mi, point ma)
{
	splitaxis = -1;
	cmax = ma;
	cmin = mi;
	data = 0;
	supernode = 0;
}

void node::update_rec(point mi, point ma)
{
	cmax = ma;
	cmin = mi;
}

void multiexecute_calcmin(vector<point> v, int unsigned start, int unsigned numdim, vector<typecor>&result)
{

	for (int unsigned i = start; i < numdim+start; i++)
	{
		priority_queue<typecor, vector< typecor >, greater< typecor > > PQ;
		//cout << " v size " << v.size() << endl;
		for (int unsigned j = 0; j < v.size(); j++)
		{
			PQ.push(v[j][i]);
			//cout << v[j][i] << ",";
		}
		
		result[i]=(PQ.top()-dis);
		//cout << "MENOR " << PQ.top() << endl;
	}
	//cout << endl;
}

void multiexecute_calcmax(vector<point> v, int unsigned start, int unsigned numdim,vector<typecor>&result)
{
	
	for (int unsigned i = start; i < numdim + start; i++)
	{
		priority_queue<typecor, vector< typecor >, less< typecor > > PQ;
		for (int unsigned j = 0; j < v.size(); j++)
		{
			PQ.push(v[j][i]);
		}
		result[i] = (PQ.top()+dis);
	}
}

void node::set_rec()
{

	cmin.resize(0);
	cmax.resize(0);
	

	if (data != 0)
	{

		for (unsigned int i = 0; i < data->size(); i++)
		{
			cmin.push_back((*data)[i] - dis);
			cmax.push_back((*data)[i] + dis);
		}
		
	}
	else
	{	
		cmin.resize(children[0]->cmin.size());
		cmax.resize(children[0]->cmin.size());
		vector <point> mins;
		vector <point> maxs;
		for (int unsigned i = 0; i < children.size(); i++)
		{
			mins.push_back(children[i]->cmin);
			maxs.push_back(children[i]->cmax);
		}
		vector <thread> threads;
		int np = 13;//thread::hardware_concurrency();
		int inc = mins[0].size() / np;

		
		for (int i = 0, j = 0; j < inc; i += np, j++) {
			threads.push_back(std::thread(&multiexecute_calcmin, mins,i,np, ref(cmin)));

		}
		for (int i = 0, j = 0; j < inc; i += np, j++) {

			threads.push_back(std::thread(&multiexecute_calcmax, maxs, i, np, ref(cmax)));

		}
		for (auto& th : threads) th.join();
		
		
	}
}

node::~node()
{
	//dtor
}
