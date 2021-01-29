#include <vector>
#include <tuple>
#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <memory>

using namespace std;

#include "MoranFunctions.h"

simVolume::simVolume(){
	this->cellFilled=0;
	this->boundReceptors=0;
}

simVolume::simVolume(volParams setupParams){
	this->cellFilled=setupParams.filled;
	this->boundReceptors=0;
	this->cellAlive=setupParams.aliveCell;
}

int simVolume::returnCellType(){
	return 0;
}

cancerCell::cancerCell(volParams setupParams):simVolume{setupParams}{
	this->cellType=setupParams.cellType;
	this->apopThreshold=setupParams.boundReceptorThresh;
}

int cancerCell::returnCellType(){
	return this->cellType;
}

bool simVolume::checkBoundReceptors(){
	return false;
}

bool cancerCell::checkBoundReceptors(){
	if(this->boundReceptors>apopThreshold){
		return true;
	}
	else{
		return false;
	}
}

void diffuseProteins(vector<vector<vector<int> > >& proteinField, boost::mt19937& randGenerator){
	vector<vector<vector<int> > > holdField(proteinField.size(),vector<vector<int> > (proteinField[0].size(),vector<int> (proteinField[0][0].size(),0)));
	int outOfTop(0);
	for(int i=0;i<(int)proteinField.size();i++){
		for(int j=0;j<(int)proteinField[i].size();j++){
			for(int k=0;k<(int)proteinField[i][j].size();k++){
				diffuseRoutine(holdField,proteinField[i][j][k],make_tuple(i,j,k),randGenerator);
			}
		}
	}
	proteinField=holdField;
}

void diffuseRoutine(vector<vector<vector<int> > >& inField, int numToDiffuse, tuple<int,int,int> currentPlace, boost::mt19937& randGenerator){
	vector<tuple<int,int,int> > curNeighbors=generateCubicNeighbors(currentPlace,inField.size());
	for(int particle=0;particle<numToDiffuse;particle++){
		double rand1((double)randGenerator()/(double)randGenerator.max());
		rand1=floor(rand1*2*curNeighbors.size());
		if(rand1<curNeighbors.size()){
			auto[xTar,yTar,zTar]=curNeighbors[rand1];
			inField[xTar][yTar][zTar]+=1;
		}
		else{
			auto[xCur,yCur,zCur]=currentPlace;
			inField[xCur][yCur][zCur]+=1;
		}
	}
}

vector<tuple<int,int,int> > generateCubicNeighbors(tuple<int,int,int> currentPosition, int cubeSize){
	vector<tuple<int,int,int> > outNeighbors(6);
	auto[xCur,yCur,zCur]=currentPosition;
	outNeighbors[0]=make_tuple(iMinus(xCur),yCur,zCur);
	outNeighbors[1]=make_tuple(iPlus(xCur,cubeSize),yCur,zCur);
	outNeighbors[2]=make_tuple(xCur,iMinus(yCur),zCur);
	outNeighbors[3]=make_tuple(xCur,iPlus(yCur,cubeSize),zCur);
	outNeighbors[4]=make_tuple(xCur,yCur,iMinus(zCur));
	outNeighbors[5]=make_tuple(xCur,yCur,iPlus(zCur,cubeSize));
	return outNeighbors;
}

int iPlus(int currentPlace, int maxSize){
	if(currentPlace==maxSize-1) return (currentPlace);
	else return (currentPlace+1);
}

int iMinus(int currentPlace){
	if (currentPlace==0) return currentPlace;
	else return (currentPlace-1);
}

void growthRound(vector<vector<vector<unique_ptr<simVolume> > > >& inVolume, volParams& newCancerParams){
	double growthProb(0.025);
	const int simDim(inVolume.size());

	for(int i=0;i<simDim;i++){
		for(int j=0;j<simDim;j++){
			for(int k=0;k<simDim;k++){
				if(inVolume[i][j][k]->cellAlive){
					if(inVolume[i][j][k]->returnCellType()==4){
						double growthPull(randPull());
						if(inVolume[iPlus(i,simDim)][j][k]->returnCellType()==0&&growthPull<growthProb){
							inVolume[iPlus(i,simDim)][j][k]=make_unique<cancerCell>(newCancerParams);
						}
						if(inVolume[iMinus(i)][j][k]->returnCellType()==0&&growthPull<growthProb){
							inVolume[iMinus(i)][j][k]=make_unique<cancerCell>(newCancerParams);
						}
						if(inVolume[i][iPlus(j,simDim)][k]->returnCellType()==0&&growthPull<growthProb){
							inVolume[i][iPlus(j,simDim)][k]=make_unique<cancerCell>(newCancerParams);
						}
						if(inVolume[i][iMinus(j)][k]->returnCellType()==0&&growthPull<growthProb){
							inVolume[i][iMinus(j)][k]=make_unique<cancerCell>(newCancerParams);
						}
						if(inVolume[i][j][iPlus(k,simDim)]->returnCellType()==0&&growthPull<growthProb){
							inVolume[i][j][iPlus(k,simDim)]=make_unique<cancerCell>(newCancerParams);
						}
						if(inVolume[i][j][iMinus(k)]->returnCellType()==0&&growthPull<growthProb){
							inVolume[i][j][iMinus(k)]=make_unique<cancerCell>(newCancerParams);
						}
					}
				}
			}
		}
	}
}

void bindingRound(vector<vector<vector<unique_ptr<simVolume> > > >& inVolume, vector<vector<vector<int> > >& cytoField, double deltaT){
	for(int i=0;i<(int)inVolume.size();i++){
		for(int j=0;j<(int)inVolume[i].size();j++){
			for(int k=0;k<(int)inVolume[i][j].size();k++){
				if(inVolume[i][j][k]->cellAlive&&inVolume[i][j][k]->returnCellType()!=0){
					cytoField[i][j][k]-=inVolume[i][j][k]->bindIL(cytoField[i][j][k],deltaT);
				}
				if(inVolume[i][j][k]->returnCellType()==4&&inVolume[i][j][k]->checkBoundReceptors()==true){
					inVolume[i][j][k]->cellAlive=false;
					inVolume[i][j][k]->cellType=5;
				}
			}
		}
	}
}

int simVolume::bindIL(int cytoIn, double deltaT){
	double timeDif(0);
	double bindProp(0.1);
	int boundCytos(0);
	do{
		if(cytoIn>0){
			timeDif+=1/(cytoIn*bindProp)*log(1/randPull());
			boundCytos+=1;
			this->boundReceptors+=1;
			cytoIn-=1;
		}
		else{
			timeDif=deltaT+1;
		}
	}while(timeDif<deltaT);
	return boundCytos;
}