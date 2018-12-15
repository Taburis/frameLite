
#include "node.h"

class graphBase {
		public : graphBase(){};
				 virtual ~graphBase(){};
				 void forwardPropagate(nodeBase &nd, int nstep = 0){
						 if(nd.isReady()) {
								 nd.is_done_=1;
								 forward_chain.push_back(&nd);
						 }
						 else return;
						 //		cout<<nd.name_<<": "<<nd.is_done_<<", ";
						 if(nd.down_stream.size() == 0){
								 //				cout<<endl;
								 return;
						 }
						 for(auto it= nd.downStream().begin(); it!= nd.downStream().end(); ++it){
								 forwardPropagate((**it), nstep+1);
						 }
				 }
				 size_t makeChain(nodeBase & nd){
						 forward_chain.clear();
						 forwardPropagate(nd);
						 resetChain();
						 return forward_chain.size();
				 }
				 void resetChain(){
						 for(auto it= forward_chain.begin(); it!= forward_chain.end(); ++it){
								 (*it)->reset();
						 }
				 }
				 void printChain(){
						 for(auto it= forward_chain.begin(); it!= forward_chain.end(); ++it){
								 std::cout<<(*it)->name_<<", ";
						 }
						 std::cout<<"END"<<std::endl;
				 }
				 void addNode(const char *key, nodeBase &nd){ node_list[key] = &nd; }
				 void removeNode(const char* key){ node_list.erase(key); }
				 nodeBase *getGraphHeader(){
						 nodeBase *ptr = NULL;
						 for(auto it = node_list.begin(); it!= node_list.end(); ++it){
								 if((*it).second->Ndependency()!=0)continue;
								 if(ptr !=0) {
										 std::cout<<"not unique header in graphBase"<<std::endl;
										 return NULL;
								 }
								 else ptr = (*it).second;
						 }
						 return ptr;
				 }
				 void append(const char * key1, const char * key2){
						 node_list[key1]->connect(*(node_list[key2]), 1);
				 }
				 void run_forward_chain (){
						 auto header_ = getGraphHeader();
						 //auto ndList = forwardPropagate(*header_);
						 for(auto it = forward_chain.begin(); it!= forward_chain.end(); it++){
								 (*it)->evaluate();
						 }
				 }
				 int evaluate(){ run_forward_chain(); return 0;}

		public : 
				 nodeBase* graphBase_header;
				 std::unordered_map<std::string, nodeBase *> node_list;
				 std::vector<nodeBase *> forward_chain;
				 std::unordered_map<long , nodeBase *> id_map;
};

class graph: public graphBase{
		public : graph(){};
				 virtual ~graph(){};
				 nodeBase & find_header(nodeBase & np){
						 //find out the header by input part of the node
						 if(np.up_stream.size() == 0) return np;
						 else return find_header( *(np.up_stream.at(0)));
				 }
				 graph( node & nd){ makeChain(find_header(nd));}
};
