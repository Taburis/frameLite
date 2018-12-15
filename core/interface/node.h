
#ifndef XFLOW_NODE_H
#define XFLOW_NODE_H
#include <vector>
#include <iostream>
#include <string>
#include "ParaSet.h"

class nodeBase {
		public : nodeBase(){};
				 ~nodeBase(){};
				 //return 1 : sccuess full evaluted the node, other wise return -1 for error
				 virtual int evaluate(){ return 0;}
				 void connect(nodeBase &, int );
				 void reset(){ is_done_ = 0;};
				 bool isReady();
				 void _done(){is_done_ = 1;}
				 int Ndependency(){ return up_stream.size();}
				 std::vector<nodeBase*>& dependency(){ return up_stream;}
				 std::vector<nodeBase*>& downStream(){ return down_stream;}

				 std::string name_;
				 bool is_done_ = 0;
		public:
				 // the up & down streams stands for the nodes connect to / from this node.
				 // the beginning node has up_stream.size()==0
				 // the endding node has down_stream.size()==0
				 std::vector<nodeBase *> up_stream;
				 std::vector<nodeBase *> down_stream;
};

bool nodeBase::isReady(){
		if(up_stream.size()==0) return 1;
		for(auto it = up_stream.begin(); it!=up_stream.end(); ++it){
				if((*it)->is_done_ ==0) return 0;
		}
		return 1;
}

void nodeBase::connect(nodeBase & nd, int type){
		// type = 1 connected out 
		//      =-1 connected in
		if( type ==1){
			   	down_stream.push_back(&nd);
			   	nd.up_stream.push_back(this);
		}
		else if(type ==-1){
			   	up_stream.push_back(&nd);
				nd.down_stream.push_back(this);
		}
		else std::cout<<" type error: "<<type<<std::endl;
}

class node : public nodeBase{
		public : node(){};
				 ~node(){};
				 void add_cfg(ParaSet & ps) {pset_ = &ps; }
				 node & operator+(node& rhs){
						 cout<<this->name_<<" + "<<rhs.name_<<endl;
						 connect(rhs, 1);
						 return rhs;
				 }
		public : 
				 ParaSet *pset_ ; 
};

#endif
