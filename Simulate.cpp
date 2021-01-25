#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <memory>

using namespace std;

#include <boost/random/mersenne_twister.hpp>
boost::mt19937 generator;

#include "MoranFunctions.h"

double randPull(){
	return (double)generator()/(double)generator.max();
}

int main(){
	generator.seed(time(NULL));


	const int simDim(50);
	volParams initParams;
	initParams.filled=false;
	initParams.cellType=4;

	vector<vector<vector<unique_ptr<simVolume> > > > simSpace;

	simSpace.resize(simDim);
	for(int i=0;i<simDim;i++){
		simSpace[i].resize(simDim);
		for(int j=0;j<simDim;j++){
			for(int k=0;k<simDim;k++){
				simSpace[i][j].push_back(make_unique<simVolume>());
			}
		}
	}

	vector<vector<vector<int> > > proteinField(simDim,vector<vector<int> > (simDim,vector<int> (simDim,0)));
	

	const int seedPatchSize(10);
	const int patchOffset(20);
	for(int i=patchOffset;i<patchOffset+seedPatchSize;i++){
		for(int j=patchOffset;j<patchOffset+seedPatchSize;j++){
			simSpace[i][j][0]=make_unique<cancerCell>(initParams);
		}
	}

	double growthProb(0.025);

	for(int growthRounds=0;growthRounds<50;growthRounds++){
		growthRound(simSpace,initParams);
	}

	ofstream testOut("cancerCells.txt");

	for(int i=0;i<simDim;i++){
		for(int j=0;j<simDim;j++){
			for(int k=0;k<simDim;k++){
				if(simSpace[i][j][k]->returnCellType()==4){
					testOut<<i<<" "<<j<<" "<<k<<" "<<k<<endl;
				}
			}
		}
	}

	testOut.close();

	proteinField[25][25][25]=2000;
	for(int i=0;i<8;i++){
		diffuseProteins(proteinField,generator);
	}



	return 0;
}