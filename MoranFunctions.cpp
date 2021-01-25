#include <vector>
#include <tuple>
#include <boost/random/mersenne_twister.hpp>

using namespace std;

#include "MoranFunctions.h"



simVolume::simVolume(volParams setupParams){
	this->cellFilled=setupParams.filled;
}

cancerCell::cancerCell(volParams setupParams):simVolume(setupParams){
	this->hold=setupParams.cellType;
}

int cancerCell::returnCellType(){
	return this->hold;
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
	if(currentPlace==maxSize) return (currentPlace);
	else return (currentPlace+1);
}

int iMinus(int currentPlace){
	if (currentPlace==0) return currentPlace;
	else return (currentPlace-1);
}