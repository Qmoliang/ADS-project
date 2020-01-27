#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <math.h>
using namespace std;

struct Building{
	int build_No; 
	int exec_time; 
	int tot_time; 
	int minheap_pos; 
	Building(int build_No, int tot_time); 
};

Building::Building(int number, int time){
	build_No = number;
	exec_time = 0;
	tot_time = time;
	minheap_pos = -1;
}

//minheap
struct MinHeap{
	vector<Building *> build_queue; 
	void insert(Building *);
	void remove(Building *);
	void increasekey(Building *, int addvalue);
};

int minheapleft(int par){return 2*par + 1;}
int minheapright(int par){return 2*(par + 1);}
int minheapparent(int child){return floor((child - 1)/2);}

//topdown operation for delete and increasekey
void topdown(vector<Building *> & build_queue, Building * cur_Building){
	int cur_pos, left_pos, right_pos;
	cur_pos = cur_Building->minheap_pos;
	while(1){ 
		left_pos = minheapleft(cur_pos);
		right_pos = minheapright(cur_pos);
		if(left_pos < build_queue.size()){ 
			if(build_queue.at(left_pos)->exec_time < build_queue.at(cur_pos)->exec_time)
				cur_pos = left_pos;
			if(right_pos < build_queue.size() && build_queue.at(right_pos)->exec_time < build_queue.at(cur_pos)->exec_time)
				cur_pos = right_pos;
			if(cur_pos != cur_Building->minheap_pos){
				swap(build_queue.at(cur_pos), build_queue.at(cur_Building->minheap_pos));
				swap(build_queue.at(cur_pos)->minheap_pos, build_queue.at(cur_Building->minheap_pos)->minheap_pos);
			}
			else break;
		}
		else break; 
	}
}

void MinHeap::insert(Building * newtask){
	build_queue.push_back(newtask);
	if(build_queue.size() > 2000)
		throw length_error("Exceeding Maximum Load");

	newtask->minheap_pos = build_queue.size() - 1; 
	int pt_pos;
	while(newtask->minheap_pos > 0){ 
		pt_pos = minheapparent(newtask->minheap_pos);
		Building * parentBuilding = build_queue.at(pt_pos);
		if(newtask->exec_time < parentBuilding->exec_time){ 
			swap(build_queue.at(newtask->minheap_pos), build_queue.at(pt_pos));
			swap(newtask->minheap_pos, parentBuilding->minheap_pos);
		}
		else break;
	}
}

void MinHeap::remove(Building * oldtask){
	Building * cur_Building = build_queue.at(build_queue.size() - 1); 
	swap(build_queue.at(oldtask->minheap_pos), build_queue.at(build_queue.size() - 1));
	swap(build_queue.at(oldtask->minheap_pos)->minheap_pos, build_queue.at(build_queue.size() - 1)->minheap_pos);
	build_queue.pop_back();
	topdown(build_queue, cur_Building);
}

void MinHeap::increasekey(Building * toincrease, int addvalue){
	toincrease->exec_time += addvalue;
	topdown(build_queue, toincrease);
}


// RBT
enum colortype {red, black};
struct RBNode{
	enum colortype color;
	Building* building;
	RBNode* parent;
	RBNode* left;
	RBNode* right;

	RBNode* sibling();
	RBNode(Building *);//store the information about building
};

RBNode* RBNode::sibling(){
	if(parent==nullptr) return nullptr;
	if(this == parent->left) return parent->right;
	else return parent->left;
}

RBNode::RBNode(Building * build){
	color = red;
	building = build;
	parent = nullptr;
	left = nullptr;
	right = nullptr;
}

Building* nullBuilding = new Building(0, 0);
RBNode* nullNode = new RBNode(nullBuilding);

struct RBTree{
	RBNode* root;

	void insert(Building*);
	void remove(Building*);
	RBTree();
	RBNode* search_No(int build_No);
	vector<Building*> search_between(int build_No1, int build_No2);	
};

RBTree::RBTree(){
	root = nullptr;
}

RBNode* RBTree::search_No(int build_No){// search key
	RBNode* getnode = root;
	while(1){
		if(getnode == nullptr){ 
			return nullNode;
		}
		if(build_No < getnode->building->build_No) 
			getnode = getnode->left;
		else if(build_No > getnode->building->build_No) 
			getnode = getnode->right;
		else
			return getnode;
	}
}

void search_between_loop(int build_No1, int build_No2, RBNode* root, vector<Building*> & build_queue){
	if(root == nullptr) return;
	if(root->building->build_No > build_No1) 
		search_between_loop(build_No1, build_No2, root->left, build_queue);
	if(root->building->build_No >= build_No1 && root->building->build_No <= build_No2) 
		build_queue.push_back(root->building);
	if(root->building->build_No < build_No2) 
		search_between_loop(build_No1, build_No2, root->right, build_queue);
}

vector<Building*> RBTree::search_between(int build_No1, int build_No2){
	vector<Building*> build_queue;
	search_between_loop(build_No1, build_No2, root, build_queue);
	if(build_queue.size() == 0) build_queue.push_back(nullBuilding);
	return build_queue;
}

RBNode* removeDegree1(RBNode* oldtask, RBNode* toraise, RBTree* tree){
	RBNode* v;
	if(oldtask->parent != nullptr){ 
		if(oldtask->parent->left == oldtask){
			oldtask->parent->left = toraise;
			v = oldtask->parent->right;
		}
		else{
			oldtask->parent->right = toraise;
			v = oldtask->parent->left;
		}
	}
	else{
		tree->root = toraise; 
		v = nullptr;
	}
	if(toraise != nullptr){
		toraise->parent = oldtask->parent;
		toraise->color = black;
	}
	return v;
}
void raiseup(RBNode * toraise, RBNode * gpnode, RBTree * tree){
	RBNode * sib = removeDegree1(gpnode, toraise, tree);
	gpnode->parent = toraise;
	gpnode->color = red;
}
void LLb(RBNode * ppnode, RBNode * gpnode, RBTree * tree){
	raiseup(ppnode, gpnode, tree);
	RBNode* c = ppnode->right;
	gpnode->left = c;
	ppnode->right = gpnode;
	if(c != nullptr) c->parent = gpnode;
}
void RRb(RBNode * ppnode, RBNode * gpnode, RBTree * tree){
	raiseup(ppnode, gpnode, tree);
	RBNode * c = ppnode->left;
	gpnode->right = c;
	ppnode->left = gpnode;
	if(c != nullptr) c->parent = gpnode;
}
void LRb(RBNode * pnode, RBNode * ppnode, RBNode * gpnode, RBTree * tree){
	raiseup(pnode, gpnode, tree);
	ppnode->parent = pnode;
	RBNode * b = pnode->left;
	RBNode * c = pnode->right;
	ppnode->right = b;
	gpnode->left = c;
	pnode->left = ppnode;
	pnode->right = gpnode;
	if(b != nullptr) b->parent = ppnode;
	if(c != nullptr) c->parent = gpnode;
}
void RLb(RBNode * pnode, RBNode * ppnode, RBNode * gpnode, RBTree * tree){
	raiseup(pnode, gpnode, tree);
	ppnode->parent = pnode;
	RBNode * b = pnode->right;
	RBNode * c = pnode->left;
	ppnode->left = b;
	gpnode->right = c;
	pnode->right = ppnode;
	pnode->left = gpnode;
	if(b != nullptr) b->parent = ppnode;
	if(c != nullptr) c->parent = gpnode;
}

void insertBinarySearchTree(RBNode *& root, RBNode * newnode){
	if(root == nullptr){ 
		root = newnode;
		return;
	}
	if(newnode->building->build_No == root->building->build_No) 
		throw invalid_argument("building number duplicated: "+to_string(newnode->building->build_No));
	newnode->parent = root;
	if(newnode->building->build_No < root->building->build_No) 
		insertBinarySearchTree(root->left, newnode);
	if(newnode->building->build_No > root->building->build_No) 
		insertBinarySearchTree(root->right, newnode);
}
void four_type_totation(RBNode * pnode, RBNode * ppnode, RBNode * gpnode, RBTree * tree){
	if(ppnode == gpnode->left && pnode == ppnode->left) LLb(ppnode, gpnode, tree);
	if(ppnode == gpnode->right && pnode == ppnode->right) RRb(ppnode, gpnode, tree);
	if(ppnode == gpnode->left && pnode == ppnode->right) LRb(pnode, ppnode, gpnode, tree);
	if(ppnode == gpnode->right && pnode == ppnode->left) RLb(pnode, ppnode, gpnode, tree);
}
void adjust_bothred(RBNode * pnode, RBTree * tree){
	RBNode * ppnode = pnode->parent;
	if(ppnode == nullptr) { 
		pnode->color = black;
		return;
	}
	RBNode * gpnode = ppnode->parent;
	if(gpnode == nullptr) { 
		ppnode->color = black;
		return;
	}
	if(pnode->color == red && ppnode->color == red){
		if(ppnode->sibling() != nullptr){
			if(ppnode->sibling()->color == red){ 
				ppnode->color = black;
				gpnode->color = red;
				ppnode->sibling()->color = black;
				adjust_bothred(gpnode, tree);
			}
			else four_type_totation(pnode, ppnode, gpnode, tree);
		}
		else four_type_totation(pnode, ppnode, gpnode, tree);
	}
}
void RBTree::insert(Building * newtask){
	RBNode * newnode = new RBNode(newtask);
	insertBinarySearchTree(root, newnode);
	adjust_bothred(newnode, this);
}


void adjust_deficient(RBNode * py, RBNode * v, RBTree * tree){
	if(py == nullptr) return; // y is root
	RBNode *w, *a, *b, *c, *x;
	enum colortype pycolor = py->color;
	bool Rxx = (v == py->left); // 1 for RXx, 0 for LXx
	w = Rxx? v->right : v->left;
	if(v->color == black){
		a = Rxx? v->left : v->right;
		if(w != nullptr && w->color == red){ // Lb2, Rb2 or Lb1, Rb1 case 2
			if(Rxx) LRb(w, v, py, tree);
			else RLb(w, v, py, tree);
			w->color = pycolor;
			py->color = black;
			return;
		}
		else{
			if(a != nullptr && a->color == red){ // Lb1, Rb1 case 1
				if(Rxx) LLb(v, py, tree);
				else RRb(v, py, tree);
				v->color = pycolor;
				py->color = black;
				a->color = black;
				return;
			}
			else{ // Lb0, Rb0
				if(pycolor == black){ // Lb0, Rb0 case 1
					v->color = red;
					if(py->parent == nullptr) return;
					else adjust_deficient(py->parent, py->sibling(), tree);
				}
				else{ // Lb0, Rb0 case 2
					v->color = red;
					py->color = black;
					return;
				}
			}
		}
	}
	else{ // v red
		a = Rxx? w->left : w->right;
		x = Rxx? w->right : w->left;
		if(x != nullptr && x->color == red){ // Lr(2), Rr(2) or Lr(1), Rr(1) case 2
			if(Rxx) LRb(x, v, py, tree);
			else RLb(x, v, py, tree);
			if(Rxx) {
				w->right = v->right;
				if(v->right != nullptr) v->right->parent = w;
				v->right = w;
			}
			else{
				w->left = v->left;
				if(v->left != nullptr) v->left->parent = w;
				v->left = w;
			}
			w->parent = v;
			py->color = black;
			return;
		}
		else{
			if(a != nullptr && a->color == red){ // Lr(1), Rr(1) case 1
				if(Rxx) LRb(w, v, py, tree);
				else RLb(w, v, py, tree);
				py->color = black;
				a->color = black;
				return;
			}
			else{ // Lr(0), Rr(0)
				if(Rxx) LLb(v, py, tree);
				else RRb(v, py, tree);
				py->color = black;
				w->color = red;
				return;
			}
		}

	}
}
RBNode * replace(RBNode * alternative){
	RBNode * thisnode = alternative->left;
	while(1){
		if(thisnode->right == nullptr){ // swap thisnode and alternative
			swap(alternative->building, thisnode->building);
			return thisnode;
		}
		else thisnode = thisnode->right;
	}
}
void RBTree::remove(Building * oldtask){
	RBNode * removenode = search_No(oldtask->build_No);
	RBNode * toraise;
	bool adjust = false; // whether need to continue adjusting
	if((removenode->left == nullptr) ^ (removenode->right == nullptr)){ // degree 1 node
		if(removenode->left != nullptr) toraise = removenode->left; // degree 1 node with left child
		else if(removenode->right != nullptr) toraise = removenode->right; // degree 1 node with right child
	}
	else{ // leaf or degree 2
		if(removenode->left != nullptr) removenode = replace(removenode); // degree 2 node
		toraise = removenode->left; // can be null
	}
	if(removenode->color == black){ // might neet to rebalance
		if(toraise == nullptr) adjust = true;
		else if(toraise->color == black) adjust = true;
	}
	RBNode * v = removeDegree1(removenode, toraise, this);
	if(adjust) adjust_deficient(removenode->parent, v, this);
}


enum operation{Insert, Print};
struct IOinterface{
	int gtime; // global time
	enum operation InorPr; // Insert or Print
	vector<int> params; // {build_No, tot_time} for Insert; {build_No1(, build_No2)} for Print
};


void get_task(string & line, IOinterface * thisterm){// convert a line of input command into an IOinterface 
	//   InorPr gtime
	thisterm->gtime = stoi(line.substr(0, line.find_first_of(":")));
	//   InorPr name
	if(line.find("Insert")!=string::npos) thisterm->InorPr = Insert;
	if(line.find("Print")!=string::npos) thisterm->InorPr = Print;
	//   InorPr parameters
	int left_pos = line.find_first_of("(") + 1;
	int right_pos = line.find_last_of(")");
	string parastr(line, left_pos, right_pos - left_pos); // substring with parameters
	int start = 0; 
	int end;
	for (size_t i = 0; i < parastr.size(); i++){
		if(parastr.at(i) == ','){
			end = i;
			thisterm->params.push_back(stoi(parastr.substr(start,end)));
			start = i+1;
		}
		if(i == parastr.size()-1)
			thisterm->params.push_back(stoi(parastr.substr(start,i+1)));
	}
}
void readinput(string filename, vector<IOinterface *> & sechedules){
	string line;
	ifstream infile(filename);
	if(infile.is_open()){
		while(getline(infile,line)){
			IOinterface * thisterm = new IOinterface;
			get_task(line, thisterm);
			sechedules.push_back(thisterm);
		}
		infile.close();
	}
}
void printbuilding(Building * building, ofstream & outfile){
	if(outfile.is_open()){
		outfile << "(" << building->build_No << ","
		                 << building->exec_time << ","
		                 << building->tot_time << ")";
	}
}
void printbuildings(vector<Building *> & build_queue, ofstream & outfile){
	if(outfile.is_open()){
		for (size_t i = 0; i < build_queue.size(); i++){
			printbuilding(build_queue.at(i), outfile);
			if(i < build_queue.size() - 1) outfile << ",";
		}
		outfile << endl;
	}
}
void printcomplete(Building * building, int global_time, ofstream & outfile){
	if(outfile.is_open()){
		outfile << "(" << building->build_No << ","
		                 << global_time << ")" << endl;
	}
}

// the city
struct City{
	vector<IOinterface *> sechedules; // the list of input operation
	MinHeap * executedTimes;
	RBTree * buildingNums;
	City();
	void ptbd(int build_No, ofstream &); // prints the triplet (build_No, exec_time, tot_time)
	void ptbds(int build_No1, int build_No2, ofstream &); // prints all triplets (build_No, exec_time, tot_time) for build_No in [build_No1, build_No2]
	void itbd(int build_No, int tot_time); // where build_No is different from existing building numbers and exec_time = zero
	Building * selectbuilding(); // select the building with the lowest exec_time to work on
	void take_one_plan(IOinterface *, ofstream &); // perform one InorPr of IOinterface
	void takeschedule(int & global_time, int end_time, Building * cur_building, ofstream &); // perform operation on the sechedules with gtime in (start_time, end_time]
};

City::City(){
	sechedules.clear();
	buildingNums = new RBTree();
	executedTimes = new MinHeap;
}
void City::ptbd(int build_No, ofstream & outfile){
	RBNode * toprint = buildingNums->search_No(build_No);
	printbuilding(toprint->building, outfile);
	outfile << endl;
}
void City::ptbds(int build_No1, int build_No2, ofstream & outfile){
	vector<Building *> toprint = buildingNums->search_between(build_No1, build_No2);
	printbuildings(toprint, outfile);
}
void City::itbd(int build_No, int tot_time){
	Building * newtask = new Building(build_No, tot_time);
	buildingNums->insert(newtask);
	executedTimes->insert(newtask);
}

vector<int> breaktie(vector<Building *> & build_queue, int root, int min_time){
	vector<int> leftmin = {};
	vector<int> rightmin = {};
	vector<int> globalmin = {root, build_queue.at(root)->build_No};
	int left_pos = minheapleft(root);
	int right_pos = minheapright(root);
	if(left_pos < build_queue.size())
		if(build_queue.at(left_pos)->exec_time == min_time){ 
			leftmin = breaktie(build_queue, left_pos, min_time);
			if(leftmin.at(1) < globalmin.at(1)) globalmin = leftmin;
		}
	if(right_pos < build_queue.size())
		if(build_queue.at(right_pos)->exec_time == min_time){ 
			rightmin = breaktie(build_queue, right_pos, min_time);
			if(rightmin.at(1) < globalmin.at(1)) globalmin = rightmin;
		}
	return globalmin;
}
Building * City::selectbuilding(){
	if(executedTimes->build_queue.size() == 0) return nullptr;
	if(executedTimes->build_queue.size() == 1) return executedTimes->build_queue.at(0);
	
	vector<int> selected; 
	selected = breaktie(executedTimes->build_queue, 0, executedTimes->build_queue.at(0)->exec_time);
	return executedTimes->build_queue.at(selected.at(0));
}
void City::take_one_plan(IOinterface * term, ofstream & outfile){
	if(term->InorPr == Insert){
		itbd(term->params.at(0), term->params.at(1));
	}
	if(term->InorPr == Print){
		if(term->params.size() == 1)
			ptbd(term->params.at(0), outfile);
		if(term->params.size() == 2)
			ptbds(term->params.at(0), term->params.at(1), outfile);
	}
}
void City::takeschedule(int & global_time, int end_time, Building * cur_building, ofstream & outfile){
	while(sechedules.size() != 0 && sechedules.front()->gtime <= end_time){
		if(cur_building != nullptr) executedTimes->increasekey(cur_building, sechedules.front()->gtime - global_time); // update cur_building executed time
		take_one_plan(sechedules.front(), outfile);
		global_time = sechedules.front()->gtime;
		sechedules.erase(sechedules.begin());
	}
	if(cur_building != nullptr) executedTimes->increasekey(cur_building, end_time - global_time); // update executed time of cur_building
	global_time = end_time;
}


int main(int argc, char* argv[]){
	int global_time = 0;
	int one_step_time = 5;
	int remaining_time, work_time;
	ofstream outputfile("output_file.txt");
	City * risingCity = new City();
	readinput(argv[1], risingCity->sechedules);
	while(1){
		Building * cur_building = risingCity->selectbuilding();
		if(cur_building == nullptr){ 
			if(risingCity->sechedules.size() == 0) break;
			risingCity->takeschedule(global_time, risingCity->sechedules.at(0)->gtime, cur_building, outputfile); 
		}
		else{
			remaining_time = cur_building->tot_time - cur_building->exec_time;
			work_time = (remaining_time <= one_step_time)? remaining_time : one_step_time;
			risingCity->takeschedule(global_time, global_time + work_time, cur_building, outputfile); 
			if(work_time == remaining_time){ 
				printcomplete(cur_building, global_time, outputfile);
				risingCity->executedTimes->remove(cur_building);
				risingCity->buildingNums->remove(cur_building);
			}
		}
	}
	outputfile.close();
	return 0;
}
