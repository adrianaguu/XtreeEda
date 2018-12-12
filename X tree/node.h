#ifndef NODE_H
#define NODE_H
#include <utility>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <thread>  
#include <queue>
using namespace std;
typedef double typecor;
typedef vector <typecor> point;
#define dis 0.25
typedef int dimtype;


class node
{
    public:
        node();
        virtual ~node();
        point cmin, cmax;
		dimtype splitaxis;
        vector <node*> children;
		vector<dimtype> split_history;
		bool supernode;
        point* data;
        void update_rec (point mi, point ma);
        void set_rec ();
        node (point mi, point ma);
		node(vector <node*> c);
        node (point d);

    protected:

    private:
};



#endif // NODE_H
