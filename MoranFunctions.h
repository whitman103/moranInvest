#ifndef MORAN_H
#define MORAN_H

#include <vector>
#include <tuple>


struct volParams{
	bool filled;
	int cellType;

};

class simVolume{
	public:
	bool cellFilled;
	simVolume(volParams setupParams);
	virtual int returnCellType();

	private:
};

class cancerCell: public simVolume{
	public:
	cancerCell(volParams setupParams);
	int returnCellType();
	int hold;
	private:
};

void diffuseProteins(vector<vector<vector<int> > >& proteinField, boost::mt19937& randGenerator);

void diffuseRoutine(vector<vector<vector<int> > >& inField, int numToDiffuse, tuple<int,int,int> currentPlace, boost::mt19937& randGenerator);

vector<tuple<int,int,int> > generateCubicNeighbors(tuple<int,int,int> currentPosition, int cubeSize);

int iPlus(int currentPlace, int maxSize);

int iMinus(int currentPlace);

#endif