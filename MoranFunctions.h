#ifndef MORAN_H
#define MORAN_H

#include <vector>
#include <tuple>

extern boost::mt19937 generator;

extern double randPull();

struct volParams{
	bool filled;
	int cellType;
	int boundReceptorThresh;
	bool aliveCell;
};

class simVolume{
	public:
	bool cellFilled;
	bool cellAlive;
	int boundReceptors;
	simVolume();
	simVolume(volParams setupParams);
	int bindIL(int cytoIn, double deltaT);
	int cellType;
	virtual int returnCellType();
	virtual bool checkBoundReceptors();
	private:
};

class cancerCell: public simVolume{
	public:
	cancerCell(volParams setupParams);
	int apopThreshold;
	int returnCellType();
	bool checkBoundReceptors();
	private:
};

void diffuseProteins(vector<vector<vector<int> > >& proteinField, boost::mt19937& randGenerator);

void diffuseRoutine(vector<vector<vector<int> > >& inField, int numToDiffuse, tuple<int,int,int> currentPlace, boost::mt19937& randGenerator);

vector<tuple<int,int,int> > generateCubicNeighbors(tuple<int,int,int> currentPosition, int cubeSize);

int iPlus(int currentPlace, int maxSize);

int iMinus(int currentPlace);

void growthRound(vector<vector<vector<unique_ptr<simVolume> > > >& inVolume, volParams& newCancerParams);

void bindingRound(vector<vector<vector<unique_ptr<simVolume> > > >& inVolume, vector<vector<vector<int> > >& cytoField, double deltaT);

#endif