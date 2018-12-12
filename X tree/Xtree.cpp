#include "Xtree.h"


int unsigned it;

bool bylower(node* a, node* b)
{
	return a->cmin[it] < b->cmin[it];
}

bool byupper(node* a, node* b)
{
	return a->cmax[it] < b->cmax[it];
}

Xtree::Xtree()
{
	root = new node;
	M = m = n=0;
}

Xtree::Xtree(unsigned int max_children, dimtype nn)
{
	M = max_children;
	m = M * 40 / 100;
	n = nn;
	root = new node;
	root->splitaxis = rand() % n;
}

typecor Xtree::area_overlap(point p1min, point p1max, point p2min, point p2max)
{
	typecor men, may;
	point newmin(n), newmax(n);
	for (dimtype i = 0; i < n; i++)
	{
		if (!overlapaxis(p1min, p1max, p2min, p2max, i, men, may))
			return 0;
		newmin[i] = men;
		newmax[i] = may;
	}
	return area(newmin, newmax);
}

typecor Xtree::area(point cmin, point cmax)
{
	typecor a = 1;
	for (unsigned int i = 0; i < n; i++)
	{
	
		a *= (cmax[i] - cmin[i]);
	
	}

	return fabs(a);
}

bool Xtree::overlapaxis(point p1min, point p1max, point p2min, point p2max, dimtype x, typecor& men, typecor& may)
{
	may = min(p1max[x], p2max[x]);
	men = max(p1min[x], p2min[x]);
	if (men < may)
		return 1;
	return 0;
}

bool Xtree::onrec(point p, point gmin, point gmax)
{
	if (gmin.size() > 0) {
		for (unsigned int i = 0; i < n; i++)
			if (p[i] > gmin[i] && gmax[i] > p[i])
				return 1;
	}
	return 0;
}

void Xtree::neighbors(node* p, vector <point>& points, point center, typecor radio)
{
	point limitsr(n);
	point limitsl(n);
	unsigned int k, i;
	for (k = 0; k < n; k++)
	{
		limitsr[k] = center[k] + radio;
		limitsl[k] = center[k] - radio;
	}

	vector <node*> rectangles;
	if (p->children[0]->data != 0)
	{
		for (i = 0; i < p->children.size(); i++)
		{
			typecor diss = distance_points(center, *(p->children[i]->data));
			//cout << "diss" << diss << endl;;
			if (diss <= radio)
				points.push_back(*(p->children[i]->data));
		}
		return;
	}

	for (i = 0; i < p->children.size(); i++)
	{
		if (outlimits(limitsl, limitsr, p->children[i]))
			continue;
		rectangles.push_back(p->children[i]);
	}

	for (int unsigned i = 0; i < rectangles.size(); i++)
	{
		neighbors(rectangles[i], points, center, radio);
	}
	
}

bool Xtree::outlimits(point limsl, point limsr, node* c)
{
	for (unsigned int k = 0; n; k++)
		if (c->cmax[k] < limsl[k] || c->cmin[k] > limsr[k])
			return 0;
	return 1;
}

void Xtree::chooseleaf(point d, node**& p, stack<node**>& path)
{
	path.push(p);

	if ((*p)->children.size() == 0 || (*p)->children[0]->data != 0) {
		return;
	}
	node** F = &((*p)->children[0]);
	for (int unsigned i = 1; i < (*p)->children.size(); i++)
	{
		if (distance_rec_point(d, (*p)->children[i]->cmin, (*p)->children[i]->cmax) < distance_rec_point(d, (*F)->cmin, (*F)->cmax))
			F = &((*p)->children[i]);
		else if (distance_rec_point(d, (*p)->children[i]->cmin, (*p)->children[i]->cmax) == distance_rec_point(d, (*F)->cmin, (*F)->cmax))
		{
			if (area((*F)->cmin, (*F)->cmax) > area((*p)->children[i]->cmin, (*p)->children[i]->cmax))
				F = &((*p)->children[i]);
		}
	}
	p = F;
	chooseleaf(d, p, path);
}

void Xtree::xtinsert(point d)
{
	stack<node**> path;
	node** p = &root;
	chooseleaf(d, p, path);
	(*p)->children.push_back(new node(d));
	(*p)->children[(*p)->children.size() - 1]->set_rec();
	if ((*p)->children.size() > M && (*p)->supernode==0)
	{
		clock_t t2 = clock();
		splitnode(p, path);
		clock_t t3 = clock();
		double time = (double(t3 - t2) / CLOCKS_PER_SEC);
		cout << "time total splitting: " << time << endl;
	}
	
	adjusttree(p, path,d);
	
	///cout << "time adjusting: " << time << endl;

}



void Xtree::adjusttree(node**& p, stack<node**>& path, point d)
{

	if(!onrec(d,(*path.top())->cmin, (*path.top())->cmax))
		(*path.top())->set_rec();
	
	if ((*path.top())->children.size() > M && (*path.top())->supernode==0) {
		splitnode(path.top(), path);
		(*path.top())->set_rec();
	}
	if (p == &root)
		return;
	path.pop();
	p = path.top();
	adjusttree(p, path,d);

}

bool Xtree::multiexecute_calcS(vector<node*> lower, int unsigned start, int unsigned end, node*& l, node*& r, typecor& S)
{
	for (unsigned int k = start; k <= end; k++)
	{
		vector<node*> lg1, lg2;
		lg1.insert(lg1.begin(), lower.begin(), lower.begin() + k);
		lg2.insert(lg2.begin(), lower.begin() + k, lower.end());
		node* temp = new node(lg1);
		temp->set_rec();
		node* temp1 = new node(lg2);
		temp1->set_rec();
		S = area_overlap(temp->cmin, temp->cmax, temp1->cmin, temp1->cmax);
		if (S == 0)
		{
			l = temp;
			r = temp1;
			return 1;
		}

	}
	return 0;
}

void Xtree::multidivide(vector<node*> lower, int unsigned numth,int np,int unsigned inc, node*& l, node*& r,typecor& S,bool& w)
{
	w = multiexecute_calcS(lower, m + (numth*inc), m + ((numth + 1)*inc), l, r,S);

}


bool Xtree::minimaldivide(node**& p, node*& l, node*& r, dimtype x)
{
	
	vector<node*> lower;
	lower  = (*p)->children;

	it = x;
	sort(lower.begin(), lower.end(), bylower);
	
	for (unsigned int k = m; k <= M - m ; k++)
	{
		vector<node*> lg1, lg2;
		lg1.insert(lg1.begin(), lower.begin(), lower.begin() + k);
		lg2.insert(lg2.begin(), lower.begin() + k, lower.end());
		node* temp = new node(lg1);
		temp->set_rec();
		node* temp1 = new node(lg2);
		temp1->set_rec();
		typecor S2;
		S2 = area_overlap(temp->cmin, temp->cmax, temp1->cmin, temp1->cmax);
		if (S2 ==0)
		{
			l = temp;
			r = temp1;
			return 1;
		}
	}
	return 0;
}
bool Xtree::stardivide(node**& p, node*& l, node*& r, dimtype x)
{
	//cout << "entra" << endl;
	vector<node*> lower;
	lower = (*p)->children;
	typecor S1 = maxoverlap;

	//cout << "empieza a calcular" << endl;
	//for (it = 0; it < (*p)->cmin.size(); it++)
	//{
	it = x;
		sort(lower.begin(), lower.end(), bylower);
	
		//cout << "hizo sort" << endl;
		for (unsigned int k = 1; k <= M - 2 * m + 2; k++)
		{
			vector<node*> lg1, lg2;
			lg1.insert(lg1.begin(), lower.begin(), lower.begin() + m - 1 + k);
			
			lg2.insert(lg2.begin(), lower.begin() + m - 1 + k, lower.end());
			
	
		
			node* temp = new node(lg1);
			temp->set_rec();
			node* temp1 = new node(lg2);
			temp1->set_rec();
			typecor S2;
			//S2 = margin(temp, temp1);
			S2 = area_overlap(temp->cmin, temp->cmax, temp1->cmin, temp1->cmax);
			if (S2 < S1)
			{
				S1 = S2;
				l = temp;
				r = temp1;
			}
			
		}

	//}

		if (S1 == maxoverlap)
			return 0;
		else
			return 1;
}

typecor	Xtree::margin(node* g1, node* g2)
{
	typecor p1, p2;
	p1 = p2 = 0;
	for (int unsigned i = 0; i < n; i++)
	{
		p1 += fabs(g1->cmin[i] - g1->cmax[i]);
		p2 += fabs(g2->cmin[i] - g2->cmax[i]);
	}

	return p1 + p2;
}




void Xtree::splitnode(node**& p, stack<node**>& path)
{

	node *l, *r;
	dimtype x=rand()%n;
	x = (*p)->splitaxis;
	bool sn=minimaldivide(p, l, r, x);
	if (sn == 0)
	{
		(*p)->supernode = 1;
			return;
	}
	l->splitaxis = x;
	r->splitaxis = x;
	if (p == &root)
	{
		node* aux = new node;
		aux->children.push_back(l);
		aux->children.push_back(r);
		root = aux;
		root->splitaxis = x;
	}
	else
	{
		path.pop();
		for (int unsigned e = 0; e < (*path.top())->children.size(); e++) {
			if ((*path.top())->children[e] == *p)
			{
				(*path.top())->children.erase((*path.top())->children.begin() + e);
			}
		}
		(*path.top())->children.push_back(l);
		(*path.top())->children.push_back(r);
		(*path.top())->set_rec();
		p = path.top();
		if ((*path.top())->children.size() > M)
		{
			splitnode(p, path);
		}
	}
}

void Xtree::bigger_Rec(point Amin, point Amax, point Bmin, point Bmax, point& Cmin, point& Cmax)
{
	for (int unsigned i = 0; i < n; i++)
	{
		Cmin.push_back(min(Amin[i], Bmin[i]));
		Cmax.push_back(max(Amax[i], Bmax[i]));
	}
}



/*clock_t t0 = clock();
clock_t t1 = clock();
double time = (double(t1 - t0) / CLOCKS_PER_SEC);*/


typecor Xtree::distance_points(point p1, point p2)
{
	typecor distance = 0;
	for (unsigned int i = 0; i <n; i++)
	{
		distance += (p1[i] - p2[i])*(p1[i] - p2[i]);

	}

	distance = sqrt(distance);
	return distance;
}


typecor Xtree::distance_rec_point(point d, point cmin, point cmax)
{
	point center;
	for (unsigned int i = 0; i < n; i++)
	{
		center.push_back(cmax[i] - ((cmax[i] - cmin[i]) / 2));
	}
	return distance_points(center, d);
}


Xtree::~Xtree()
{
	//dtor
}
