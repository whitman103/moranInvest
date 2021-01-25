#include <vector>
#include <iostream>
#include <chrono>
#include <fstream>
#include <memory>

using namespace std;

#include <boost/random/mersenne_twister.hpp>
boost::mt19937 generator;

#include "MoranFunctions.h"

int main(){
	generator.seed(time(NULL));


	const int simDim(50);
	volParams initParams;
	initParams.filled=false;
	initParams.cellType=0;

	vector<vector<vector<unique_ptr<simVolume> > > > simSpace(simDim,vector<vector<unique_ptr<simVolume> > >(simDim, vector<unique_ptr<simVolume> > (simDim,make_unique<simVolume>(initParams))));

	vector<vector<vector<int> > > proteinField(simDim,vector<vector<int> > (simDim,vector<int> (simDim,0)));
	proteinField[25][25][25]=2000;
	for(int i=0;i<8;i++){
		diffuseProteins(proteinField,generator);
	}

	simSpace[25][25][25]=(make_unique<cancerCell>(initParams));

	cout<<simSpace[25][25][25]->returnCellType()<<endl;



	return 0;
}