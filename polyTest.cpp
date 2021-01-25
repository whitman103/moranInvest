#include <vector>
#include <iostream>
#include <memory>

using namespace std;

class Base{
	public:
	Base();
	Base(int in);
	virtual int returnType();
	int type;
};

int Base::returnType(){
	return 0;
}

Base::Base(){
	this->type=0;
}

Base::Base(int in){
	this->type=in;
}

class Derived: public Base{
	public:
	Derived(int in);
	int returnType();
	int type;
};

Derived::Derived(int in){
	this->type=in+1;
}

int Derived::returnType(){
	return this->type;
}

int main(){

	vector<unique_ptr<Base> > objects;

	objects.emplace_back(new Base(2));
	objects.emplace_back(new Derived(3));

	for(int i=0;i<(int)objects.size();i++){
		cout<<objects[i]->returnType()<<endl;
	}

	return 0;
}