/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/


/*

This file contains a class for a matrix for adjacency history

Created the: 30-11-2015
by: Wandrille Duchemin


Last modified the: 27-10-2016
by: Wandrille Duchemin

*/

#include "AdjMatrix.h"


void AdjMatrix::backtrackAux( vector< AdjTree *> * AdjacencyTrees, bool stochastic, bool alwaysGainAtTop , double c1proba , int Root1, int Root2)
{

	double c1 = getC1(Root1,Root2);
	double c0 = getC0(Root1,Root2);

	if(verbose)
		cout << "AdjMatrix::backtrackAux " << Root1 << "-" << Root2 << "-> c1: "<< c1 << " , c0: " << c0 <<endl;
	//verbose=true;

	bool c1best;

	if ( (c1 == c0) && ( c0 == worstScore) ) // both score are equal AND worst --> overflow?
	{
		
		cerr << "Potential overflow: only impossible solutions at the root of the backtrack."<<endl;

		double r = ((double) rand()/ RAND_MAX); //random result

		if(r <  c1proba)
			c1best = true;
		else
			c1best = false;
	}
	else if(c0 == worstScore)
	{
		c1best = true;
	}
	else if(c1 == worstScore)
	{
		c1best = false;
	}
	else
	{ // neither is the worst -> choose better or randomly
		if(alwaysGainAtTop) // we have to add a gain to c1
		{
			if(useBoltzmann) // for Boltzmann DeCo
			{
				c1 *= getBoltzmannGainBreakCost(1,0); 
			}
			else
			{
				c1 += getGainCost();
			}
		}

	
		if (c1 == c0) // both score are equal --> randomly determine which to choose
		{
			
			double r = ((double) rand()/ RAND_MAX); //random result

			if(r <  c1proba)
				c1best = true;
			else
				c1best = false;
		}
		else if(stochastic)//big score has a bigger proba to be chosen. Designed to be used with Boltzmann versio only.
		{
			double r = ((double) rand()/ RAND_MAX); //random result

			if(r < (c1 / (c1 + c0) ) )
				c1best = true;
			else
				c1best = false;
		}
		else // best score is chosen ( best means minimum in classical DeCo, and maximal in BoltzmannDeCo)
		{

			c1best = false;
			if(c0 > c1)
				c1best = true;
	
			if(useBoltzmann)//reverse -> we want the greater score to win
				c1best = !c1best;
		}
	}



	if(c1best) // creation of a new AdjTree
	{
		AdjTree * ATree = new AdjTree(alwaysGainAtTop); //new adjacency tree with a gain at its root
		AdjacencyTrees->push_back(ATree);
		ATree->setRootNode(backtrackAuxC1(Root1, Root2, AdjacencyTrees, stochastic));
		//cout << "plop5" << endl;
	}
	else // C0 -> continue backtrack. if another c1 is encountered it will be part of another AdjTree anyway
		backtrackAuxC0(Root1, Root2, AdjacencyTrees, stochastic);


	//cout << "plop6" << endl;
	for(vector< AdjTree *>::iterator it = AdjacencyTrees->begin(); it != AdjacencyTrees->end(); ++it)
	{
		(*it)->resetNodesId();
	}
	//cout << "plop7  "<< "nb of trees " << AdjacencyTrees->size() << endl;
}

bool AdjMatrix::isImpossibleCase(vector<AdjSolution> aSolC1, vector<AdjSolution> aSolC0)
{
	if(aSolC1[0].score == worstScore) //  && (aSolC0[0].score == bestScore))
		return true;

	return false;
	/*
	if(aSolC1.size() !=1)
		return false;

	if((aSolC1[0].score != worstScore) || (aSolC0[0].score != bestScore))
		return false;

	if((aSolC1[0].components.size() != 0) || (aSolC0[0].components.size() != 0))
		return false;

	return true;*/
}

void AdjMatrix::initRootMatrix()
{
	RootMatrix = MyMatrix();
	int d1 = TreeToMatrixId1.size();
	int d2 = TreeToMatrixId2.size();



	RootMatrix.setDim(d1, d2);




	rootComputed = false;


}

void AdjMatrix::scanSolutionForRooting(int id1, int id2)
{

	vector<AdjSolution> C1solutions = getSolutionC1( id1, id2);
	vector<AdjSolution> C0solutions = getSolutionC0( id1, id2);

	//2. all components of the solutions can't be root.

	for(unsigned i = 0; i < C1solutions.size(); i++)
	{
		if(C1solutions[i].components.size() > 0)
		{
			for (vector<AdjScore>::iterator it = C1solutions[i].components.begin() ; it != C1solutions[i].components.end(); ++it)//iterates over the AdjScore of the chosen solution
			{
				setRootMatrixNotPossible(it->id1,it->id2);//the solution uses that id pair -> is not a root

			}
		}
	}

	for(unsigned i = 0; i < C0solutions.size(); i++)
	{
		if(C0solutions[i].components.size() > 0)
		{
			for (vector<AdjScore>::iterator it = C0solutions[i].components.begin() ; it != C0solutions[i].components.end(); ++it)//iterates over the AdjScore of the chosen solution
			{
				setRootMatrixNotPossible(it->id1,it->id2);//the solution uses that id pair -> is not a root

			}
		}
	}

}

void AdjMatrix::scanForRoots()
{
	if(!isComputed())
	{
		throw Exception("AdjMatrix::scanForRoots : scanning for root when the matrix is not computed!!");
	}
	
	int d1 = TreeToMatrixId1.size();
	int d2 = TreeToMatrixId2.size();

	//1. checking if self could be a root
	for(unsigned i = 0; i < d1; i++)
	{
		for(unsigned j = 0; j < d2 ; j++)
		{
//			vector<AdjSolution> C1solutions = getSolutionC1( MatrixToTreeId1[i], MatrixToTreeId2[j]);
//			vector<AdjSolution> C0solutions = getSolutionC0( MatrixToTreeId1[i], MatrixToTreeId2[j]);
//
//			if(isImpossibleCase(C1solutions, C0solutions))
//				setRootMatrixNotPossible(MatrixToTreeId1[i], MatrixToTreeId2[j]);
//			else
				setRootMatrixPossible(MatrixToTreeId1[i], MatrixToTreeId2[j]);
		}
	}

	for(unsigned i = 0; i < d1; i++)
	{
		for(unsigned j = 0; j < d2 ; j++)
		{
			if(getRootMatrix(MatrixToTreeId1[i] , MatrixToTreeId2[j]))
				scanSolutionForRooting(MatrixToTreeId1[i] , MatrixToTreeId2[j]);
		}
	}

	rootComputed = true;	

}


vector < pair<int,int> > AdjMatrix::getMatrixRoots()
{
	if(!rootComputed)
	{
		//cout << "AdjMatrix::getMatrixRoots : tried to get the roots while they aren't set : setting them" << endl;
		scanForRoots();
	}

	int d1 = TreeToMatrixId1.size();
	int d2 = TreeToMatrixId2.size();

	vector < pair<int,int> > roots;

	for(unsigned i = 0; i < d1; i++)
	{
		for(unsigned j = 0; j < d2 ; j++)
		{
			int id1 = MatrixToTreeId1[i];
			int id2 = MatrixToTreeId2[j];
			if(getRootMatrix(id1,id2))
			{
				roots.push_back( pair <int,int> (id1 , id2)  );
			}
		}
	}
	return roots;
}

////////////WMODIF



//// methods that will be used by the score algebra
double AdjMatrix::addition(double const& a, double const& b)
{
	return a + b;
}

double AdjMatrix::multiplication(double const& a, double const& b)
{
	return a * b;
}

double AdjMatrix::getminimum(vector <double> const& v)
{
	double m = v[0];
	for(unsigned i = 1 ; i < v.size(); i++)
	{
		if(v[i] < m)
			m = v[i];
	}
	return m;
}

double AdjMatrix::getsum(vector <double> const& v)
{
	double s = 0;
	for(unsigned i = 0 ; i < v.size(); i++)
	{
		s += v[i];
	}
	return s;
}




/*
Takes:
 - map<int,vector<float> > speciesC0C1 : used by artdeco -> cost of having an adjacency or not for different species when the adjacency was not in the input
 - map<int, map<string,int> > speGeneAdjNb : used by artdeco -> degree of each gene, used to modulate the costs in speciesC0C1.
 - vector <double> adjacencyScoreVec : used by ADseq -> score between 0 and 1 of the adjacencies that denote the confidence we have in them
 - Gcost (double): cost of a Gain
 - Bcost (double): cost of a Break
 - rtree1 (ReconciledTree *): reconciled tree for the first dimension
 - rtree2 (ReconciledTree *): reconciled tree for the second dimension
 - adjacencies (vector< pair <int,int> >)
 - VERBOSE (bool)
 - boltzmann (bool) (default: false): wether to use boltzmann computation or not
 - temp (double) : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - absencePenalty double) : if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization.
 - double adjScoreLogBase : base of the log that will be used to go from adjacency confidence score to parsimony costs
*/
void AdjMatrix::AdjMatrixAux(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, vector <double> adjacencyScoreVec, double Gcost, double Bcost, ReconciledTree * rtree1, ReconciledTree * rtree2, vector< pair <int,int> > adjacencies, bool VERBOSE, bool boltzmann , double temp , double absencePenalty, double adjScoreLogBase)
{

	GainCost = Gcost;
	BreakCost = Bcost;

	matrixComputed = false;

	verbose = VERBOSE;

	Rtree1 = *rtree1;
	Rtree2 = *rtree2;


	//setting the algebra to the basic one -> most parsimonious
	useBoltzmann = boltzmann;
	Temperature = temp;
	if((Temperature == 0) && (useBoltzmann))
	{
		if(verbose)
			cout << "Asked for Boltzmann computation but with a temperature of 0 (impossible). Switching back to original algorithm instead." << endl;
		useBoltzmann = false;
	}


	if(!useBoltzmann)
	{
		scoreAggregatorfunc = &AdjMatrix::addition;
		scoreComparatorfunc = &AdjMatrix::getminimum;

		defaultScore = MYMINFINIY;
		worstScore = MYINFINIY;
		bestScore = 0;

	}
	else
		setComputationToBoltzmann();

	if(absencePenalty == -1)
	{
		worstAbsenceScore = worstScore;//default behavior...
	}
	else if(useBoltzmann)
		worstAbsenceScore = exp( - (absencePenalty) / (BOLTZMANN_K * Temperature)) ; // boltzmann case
	else
		worstAbsenceScore = absencePenalty;


	//setting reconciliation event costs
	WeightedDupCost  = bestScore;//Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score. This is defaulting at 0 and set only if requested when computing the matrix. Used to subtract a reconciliation event to the adjacency scores of co-events.
	WeightedLossCost = bestScore;//Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score. This is defaulting at 0 and set only if requested when computing the matrix. Used to subtract a reconciliation event to the adjacency scores of co-events.
	WeightedHgtCost  = bestScore;//Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score. This is defaulting at 0 and set only if requested when computing the matrix. Used to subtract a reconciliation event to the adjacency scores of co-events.	



	//we have to build the map between matrix and node ids
	BuildIdMaps();

	initMatrix(adjacencies,speciesC0C1,speGeneAdjNb, adjacencyScoreVec, adjScoreLogBase);//set adjacencies and their absence


	if((Rtree1.getNumberOfTransfer() >0) || (Rtree2.getNumberOfTransfer() >0)) // at least one transfer in one of the tree
		decoLTalgo = true;
	else
		decoLTalgo = false;

	/////WMODIF
	initRootMatrix();

	/////WMODIF

}


/*
Builds the maps between matrix and tree ids
*/
void AdjMatrix::BuildIdMaps()
{
	//first making sure thing are cleared
	TreeToMatrixId1.clear();
	TreeToMatrixId2.clear();
	MatrixToTreeId1.clear();
	MatrixToTreeId2.clear();

	vector <int> tree1Ids = Rtree1.getNodesId();

	for(unsigned i = 0; i < tree1Ids.size(); i++)
	{
		TreeToMatrixId1[tree1Ids[i]] = i;
		MatrixToTreeId1[i] = tree1Ids[i];
	}

	vector <int> tree2Ids = Rtree2.getNodesId();

	for(unsigned i = 0; i < tree2Ids.size(); i++)
	{
		TreeToMatrixId2[tree2Ids[i]] = i;
		MatrixToTreeId2[i] = tree2Ids[i];
	}


}



/*
adds an adjacency to the matrix by setting its c1 and c0 costs (respectively to 0 and MYINFINIY)

Takes:
 - adjacency (pair < int,int >) : pair of leaf ids in tre reconciled trees which form an adjacency

CAUTION: does not check wether these id exists or are leaves.
CAUTION2: won't do anything if the matrix is computed
*/
void AdjMatrix::addAdjacency(pair <int,int> adjacency)
{

	if(matrixComputed)
	{
		if(verbose)
			cout << "Trying to add an adjacency while the matrix is already computed (matrixComputed is true). Does nothing."<< endl;
		return;
	}

	int id1 = TreeToMatrixId1[adjacency.first];
	int id2 = TreeToMatrixId2[adjacency.second];

	MatrixC0.setValue(id1,id2,worstScore);
	MatrixC1.setValue(id1,id2,bestScore);
}

/*
adds an adjacency to the matrix by setting its c1 and c0 costs (respectively to 0 and MYINFINIY)

Takes:
 - adjacency (pair < int,int >) : pair of leaf ids in tre reconciled trees which form an adjacency
 - double score : confidence score associated with the adj
 - double adjScoreLogBase : base of the log that will be used to go from adjacency confidence score to parsimony costs

CAUTION: does not check wether these id exists or are leaves.
CAUTION2: won't do anything if the matrix is computed
*/
void AdjMatrix::addAdjacency(pair <int,int> adjacency, double score, double adjScoreLogBase)
{

	if(matrixComputed)
	{
		if(verbose)
			cout << "Trying to add an adjacency while the matrix is already computed (matrixComputed is true). Does nothing."<< endl;
		return;
	}

	int id1 = TreeToMatrixId1[adjacency.first];
	int id2 = TreeToMatrixId2[adjacency.second];

	//NB: normally the case where score == 1 have already been treated 
	if(score == 1) // re-check to be sure
	{
		addAdjacency(adjacency);
	}
	else if(score != 0) 
	{
		//ADseq formulas
		double c1 =  - log( score ) / log( adjScoreLogBase);
		double c0 =  - log( 1 - score ) / log( adjScoreLogBase);


		if(useBoltzmann)
		{ // boltzmann adaptation
			c1 = exp( - c1 / Temperature);
			c0 = exp( - c0 / Temperature);
		}

		MatrixC0.setValue(id1,id2,c0);
		MatrixC1.setValue(id1,id2,c1);

		if(verbose)
			cout << "adj " << adjacency.first << "-" << adjacency.second << " -> (c1/c0) "<< c1 << "/" << c0 << endl;
	}
	else
	{
		// this adj is actually an adj absence !!
		// we don't do anything and it will be treated as such
		if(verbose)
			cout << "adj " << adjacency.first << "-" << adjacency.second << " treated as absent"<< endl;
	}
}



/*
Initialize the matrix -> sets all values at -1

CAUTION: won't do anything if the matrix is computed
*/
void AdjMatrix::initMatrix()
{
	if(matrixComputed)
	{
		if(verbose)
			cout << "Trying to initialize the matrix while it is already computed (matrixComputed is true). Does nothing."<< endl;
		return;
	}

	MatrixC0 = MyMatrix();
	MatrixC1 = MyMatrix();

	SolutionMatrixC0 = MyMatrixAdjSolution();
	SolutionMatrixC1 = MyMatrixAdjSolution();

	int d1 = TreeToMatrixId1.size();
	int d2 = TreeToMatrixId2.size();

	MatrixC0.setDim(d1, d2);//setting dimension is what initializes the matrix
	MatrixC1.setDim(d1, d2);
	
	SolutionMatrixC0.setDim(d1, d2);//setting dimensions
	SolutionMatrixC1.setDim(d1, d2);



	//by default the value of all case is -1 -> we replace it by default value (which is also -1 for now...)
	
	if(defaultScore != MatrixC0.getValue(0,0)) //default value of the matrices is not the default anymore -> we have to set the correct default score everywhere
	{
		for(unsigned i=0; i< d1; i++)
			for(unsigned j=0; j< d2; j++)
			{
				MatrixC0.setValue (i, j, defaultScore);
				MatrixC1.setValue (i, j, defaultScore);
			}
	}

	return;
}



/*
Sets the absence of adjacencies between leaves of the same species that aren't already set
*/
void AdjMatrix::initAdjAbsence(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb)
{
	// setting absence of adjacencies
	vector <int> leaves1 = Rtree1.getLeavesId();
	vector <int> leaves2 = Rtree2.getLeavesId();

	for(unsigned i = 0; i < leaves1.size(); i++)
	{
		int leaf1 = leaves1[i];

		if(!Rtree1.isRealLeaf(leaf1))//case of the loss nodes
			continue;

		for(unsigned j = 0; j < leaves2.size(); j++)
		{
			int leaf2 = leaves2[j];

			if(!Rtree2.isRealLeaf(leaf2))//case of the loss nodes
				continue;

			if(!issetC1(leaf1,leaf2)) // case where the adjacency exists -> it is already set to something other than the default value
			{
				int sp1 = Rtree1.getNodeSpecies(leaf1);
				int sp2 = Rtree2.getNodeSpecies(leaf2);
				if(sp1 == sp2)//species coincides -> setting the values
				{
					if(speciesC0C1.empty())
					{
						setC0(leaf1, leaf2, bestScore);
						setC1(leaf1, leaf2, worstAbsenceScore);
					}
					else
					{
						// Get number of adjacencies of leaf1 and leaf2 to get corresponding value of C0 and C1
						string g1=Rtree1.getNode(leaf1)->getName();
						string g2=Rtree2.getNode(leaf2)->getName();
						int adjNb1=speGeneAdjNb[sp1][g1];
						int adjNb2=speGeneAdjNb[sp2][g2];
						double scoreC0;
						double scoreC1;

						if(verbose)
							cout<<"\tGene1: "<<g1<<" -> "<<adjNb1<<" | Gene2: "<<g2<<" -> "<<adjNb2<<endl;

						if(adjNb1==adjNb2==1)
						{
							scoreC0=speciesC0C1[sp1][0];
							scoreC1=speciesC0C1[sp1][3];
						}
						else if((adjNb1==0 && adjNb2==1) || (adjNb1==1 && adjNb2==0))
						{
							scoreC0=speciesC0C1[sp1][1];
							scoreC1=speciesC0C1[sp1][4];
						}
						else if(adjNb1==adjNb2==0)
						{
							scoreC0=speciesC0C1[sp1][2];
							scoreC1=speciesC0C1[sp1][5];
						}
						else
						{
							cout<<endl<<"ERROR: one of the 2 genes "<<g1<<" or "<<g2<<" of species "<<spe<<endl;
							exit(EXIT_FAILURE);
						}

						if(useBoltzmann) // applying transformation to score in order to make them probas
						{
							scoreC0 = exp( - scoreC0 / Temperature);
							scoreC1 = exp( - scoreC1 / Temperature);
						}


						if(verbose)
							cout<<"\t\t=> ScoreC0: "<<scoreC0<<" | ScoreC1: "<<scoreC1<<endl;

						setC0(leaf1, leaf2, scoreC0);
						setC1(leaf1, leaf2, scoreC1);
					}
				}
				else
				{
					setC0(leaf1, leaf2, bestScore);
					setC1(leaf1, leaf2, worstScore);	
				}
			}
			//setting default solution for each kind of leaf couple.
			vector<AdjSolution> Vsol1 = SolutionC1ExtantWithExtant(leaf1,leaf2);
			vector<AdjSolution> Vsol0 = SolutionC0ExtantWithExtant(leaf1,leaf2);
			
			setSolutionC1(leaf1,leaf2, Vsol1);
			setSolutionC0(leaf1,leaf2, Vsol0);
		}
	}
}




/*
Set the values of C1 and C0 to respectively 0 and MYINFINIY for the cases given in adjacencies; 
Sets the values of C1 and C0 to respectively MYINFINIY and 0 between leaves not given in adjacencies;
sets all other values to -1.


Takes:
 - adjacencies (vector< pair <int,int> >): vector of pairs of ids in reconciled trees 1 and 2
 - map<int,vector<float> > speciesC0C1 : used by artdeco -> cost of having an adjacency or not for different species when the adjacency was not in the input
 - map<int, map<string,int> > speGeneAdjNb : used by artdeco -> degree of each gene, used to modulate the costs in speciesC0C1.
 - vector <double> adjacencyScoreVec : used by ADseq -> score between 0 and 1 of the adjacencies that denote the confidence we have in them
 - double adjScoreLogBase : base of the log that will be used to go from adjacency confidence score to parsimony costs

CAUTION: does not check wether these id exists or are leaves.
CAUTION2: won't do anything if the matrix is computed
*/
void AdjMatrix::initMatrix(vector< pair <int,int> > adjacencies, map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, vector <double> adjacencyScoreVec, double adjScoreLogBase	)
{
	initMatrix();

	//setting adjacencies
	addAdjacencies(adjacencies, adjacencyScoreVec,adjScoreLogBase);

	// setting absence of adjacencies
	initAdjAbsence(speciesC0C1,speGeneAdjNb);

}


/*
Sets a case in the C0 matrix

Takes:
 - id1 (int): a node id in the first reconciled tree
 - id2 (int): a node id in the second reconciled tree
 - value (double): the wished value for the case

CAUTION: does not check wether these id exists or are leaves.
CAUTION2: won't do anything if the matrix is computed
*/
void AdjMatrix::setC0(int id1, int id2, double value)
{
	if(matrixComputed)
	{
		if(verbose)
			cout << "Trying to set C0 while it is already computed (matrixComputed is true). Does nothing."<< endl;
		return;
	}

	int Id1 = TreeToMatrixId1[id1];
	int Id2 = TreeToMatrixId2[id2];

	MatrixC0.setValue(Id1,Id2,value);
}

/*
Sets a case in the C1 matrix

Takes:
 - id1 (int): a node id in the first reconciled tree
 - id2 (int): a node id in the second reconciled tree
 - value (double): the wished value for the case

CAUTION: does not check wether these id exists or are leaves.
CAUTION2: won't do anything if the matrix is computed
*/
void AdjMatrix::setC1(int id1, int id2, double value)
{
	if(matrixComputed)
	{
		if(verbose)
			cout << "Trying to set C1 while it is already computed (matrixComputed is true). Does nothing."<< endl;
		return;
	}

	int Id1 = TreeToMatrixId1[id1];
	int Id2 = TreeToMatrixId2[id2];

	MatrixC1.setValue(Id1,Id2,value);
}


/*
Reset to the default value the case given as argument

Takes:
 - id1 (int): a node id in the first reconciled tree
 - id2 (int): a node id in the second reconciled tree


*/
void AdjMatrix::resetCase(int id1, int id2)
{
	setC0(id1, id2, defaultScore);
	setC1(id1, id2, defaultScore);

}

/*
Takes
	- Vsolution (vector <Adjsolution>)

Returns:
	(AdjSolution): a randomly chosen Adjsolution in Vsolution

*/
AdjSolution AdjMatrix::chooseRandomSolution(vector <AdjSolution> Vsolution)
{
	int r =  rand() % Vsolution.size(); //random index
	return Vsolution[r];
}

/*
Takes
	- Vsolution (vector <Adjsolution>)

Returns:
	(AdjSolution): a randomly chosen Adjsolution among the best solutions of Vsolution
*/
AdjSolution AdjMatrix::chooseBestSolution(vector <AdjSolution> Vsolution)
{
	vector<int> tochoosefrom;

	for(unsigned i = 0 ; i < Vsolution.size(); i++)
		tochoosefrom.push_back(i);

	vector <int> chosen = SolutionMatrixC1.getIndexOnlyBest( Vsolution, tochoosefrom, !useBoltzmann);

/*
	cout << "AdjMatrix::chooseBestSolution " << !useBoltzmann << endl;
	for(unsigned i = 0; i < Vsolution.size(); i++ )
		cout << i << ":" << Vsolution[i].score << " ";
	cout << endl;
	cout << "chosen : ";
	for(unsigned i = 0; i < chosen.size(); i++ )
		cout << chosen[i] << " ";
	cout << endl;
*/
	if(chosen.size() == 1)
		return Vsolution[chosen[0]];

	return Vsolution[chosen[rand() % chosen.size()]];
}



/*
Takes:
	- Vsolution (vector <AdjSolution>)

Returns:
	(AdjSolution): an AdjSolution randomly chosen into Vsolution using its score as weight (or exp(-score/temperature) if non-boltzmann)

*/
AdjSolution AdjMatrix::chooseRandomSolutionWeighted(vector <AdjSolution> Vsolution)
{
	vector <double> probas;
	double total = 0;
	//cout << "AdjMatrix::chooseRandomSolutionWeighted " ;
	for(unsigned i = 0; i < Vsolution.size(); i++)
	{
		double score = Vsolution[i].score;
		if(!useBoltzmann)
			score = exp(( - score ) / (BOLTZMANN_K * Temperature)) ; 

		total += score;
		probas.push_back(score);
		//cout <<" ; " << i << " -> " << score ;

	}
	//cout << endl;

	if( total < 1) // trying to avoid overflows.
	{
		//cout << "normalizing to avoid overflow (total = "<<total<<")";

		for(unsigned i = 0 ; i < probas.size() ; i++)
		{
			probas[i]/=total;
			//cout <<" ; " << i << " -> " << probas[i] ;
		}
		//cout<<endl;
		total = 1;
	}

	double r = total * ((double) rand()/ RAND_MAX); //random results-> potential overflow if total is very small.
	//cout << "AdjMatrix::chooseRandomSolutionWeighted TOTAL -> " << total << " -> " << r;

	if((r<=0)||(r>=total)) // checking potential overflows here
	{
		cerr << "Detected overflow in random weighted backtracking. Temporarily reverting to chosing best solution in order to avoid undefined behaviours."<<endl;
		return chooseBestSolution(Vsolution);
	}

	unsigned i = 0;
	while(r > probas[i])
	{
		r -= probas[i];
		i++;
	}
	//cout << " -> chosen "<< i<<endl;
	return Vsolution[i];
}



/*
Takes:
	- id1 (int): a node id in Gfamily1 tree
	- id2 (int): a node id in Gfamily2 tree
 	- AdjacencyTrees (vector< AdjTree *> *): pointer to a vector of adjacency trees vector that will be update as we build more adjacency trees
 	- stochastic (bool): true if the backtrack is to be stochastic; false if the backtrack choses the solution with the best score

Returns:
	(Node *): node for an AdjTree object representing the adjacency between id1 and id2
*/
Node * AdjMatrix::backtrackAuxC1(int id1, int id2, vector< AdjTree *> * AdjacencyTrees, bool stochastic)
{
	if(verbose)
		cout << "In backtrackC1 " << id1 << "," << id2 << endl;

	//1. create the node corresponding to this adj

	pair <int, int> AdjEvtSpe = getAdjEventAndSpecies(id1,id2); // first part of the pair is an event code ; second part is the species

	bool isLeaf;


	if(Rtree1.isExtant(AdjEvtSpe.first))
		isLeaf = true;

	//2.choose the correct solution

	vector<AdjSolution> Vsolution = getSolutionC1( id1, id2);

	if(verbose)
	{
		for(unsigned i = 0; i < Vsolution.size(); i++)
		{
			if(Vsolution[i].score == worstScore)
				cout << "worst";
			else
				cout << Vsolution[i].score ;
			//cout << "(" << Vsolution[i].coevent << ")" << " ; " ;
		}
		//cout << endl;
	}

	AdjSolution chosenSolution;




	if(Vsolution.size() == 0 )
		throw Exception("AdjMatrix::backtrackAuxC1 : empty solution set found.");
	else if(Vsolution.size() == 1 )
	{	
		chosenSolution = Vsolution[0];
	}
	else
	{
		if(stochastic)
			chosenSolution = chooseRandomSolutionWeighted(Vsolution);
		else
			chosenSolution = chooseBestSolution(Vsolution);
	}

	if(verbose)
		cout << "chosen solution: " << chosenSolution.score << "(" << chosenSolution.coevent << "): components:" <<  chosenSolution.components.size() << " B: " << chosenSolution.NbBreak  << " G: " << chosenSolution.NbGain << endl;

//	//cout << "plop " <<  AdjacencyTrees->size() << endl;

	//3. creating the node
	Node * n = AdjacencyTrees->at(AdjacencyTrees->size() - 1)->createNode(pair <int,int> (id1,id2), AdjEvtSpe.second, chosenSolution.coevent, AdjEvtSpe.first); // NB: there is at least 1 tree in Adjacency trees

	//cout << "plop2"<< endl;

	int GainToDo = chosenSolution.NbGain; // the solution are arranged so that the first C1s are the one bearing a gain.


	if(chosenSolution.components.size() > 0)
	{
		for (vector<AdjScore>::iterator it = chosenSolution.components.begin() ; it != chosenSolution.components.end(); ++it)//iterates over the AdjScore of the chosen solution
		{
			if(it->C1) // this one is C1 -> add to the current node
			{
				if(GainToDo > 0)
				{
					//there is a new gain/tree 
					AdjTree * ATree = new AdjTree(true); //new adjacency tree with a gain at its root
					AdjacencyTrees->push_back(ATree);
					ATree->setRootNode(backtrackAuxC1(it->id1, it->id2, AdjacencyTrees, stochastic));

					//one less gain
					GainToDo--;
				}
				else // not a gain -> the tree is continued
				{
					n->addSon(backtrackAuxC1(it->id1, it->id2, AdjacencyTrees, stochastic));
				}
			}
			else // C0 -> continue backtrack. if another c1 is encountered it will be part of another AdjTree anyway
				backtrackAuxC0(it->id1, it->id2, AdjacencyTrees, stochastic);
	
	
		}
	}

	//cout << "plop3"<< endl;

	//adding supplementary AdjBreakNodes
	for(unsigned i = 0; i < chosenSolution.NbBreak; i++)
	{
		n->addSon(AdjacencyTrees->at(0)->createBreakNode());
	}
	//cout << "plop4"<< endl;
	return n;
}


/*
Takes:
	- id1 (int): a node id in Gfamily1 tree
	- id2 (int): a node id in Gfamily2 tree
 	- AdjacencyTrees (vector< AdjTree *> *): pointer to a vector of adjacency trees vector that will be update as we build more adjacency trees
 	- stochastic (bool): true if the backtrack is to be stochastic; false if the backtrack choses the solution with the best score

*/
void AdjMatrix::backtrackAuxC0(int id1, int id2, vector< AdjTree *> * AdjacencyTrees, bool stochastic)
{
	if(verbose)
		cout << "In backtrackC0 " << id1 << "," << id2 << endl;

	//choose the correct solution

	vector<AdjSolution> Vsolution = getSolutionC0( id1, id2);

	AdjSolution chosenSolution;


	if(verbose)
	{
		for(unsigned i = 0; i < Vsolution.size(); i++)
		{
			if(Vsolution[i].score == worstScore)
				cout << "worst";
			else
				cout << Vsolution[i].score ;
			//cout << "(" << Vsolution[i].coevent << ")" << " ; " ;
		}
		cout << endl;
	}




	if(Vsolution.size() == 0 )
		throw Exception("AdjMatrix::backtrackAuxC0 : empty solution set found.");
	else if(Vsolution.size() == 1 )
		chosenSolution = Vsolution[0];
	else
	{
		if(stochastic)
			chosenSolution = chooseRandomSolutionWeighted( Vsolution);
		else
			chosenSolution = chooseBestSolution( Vsolution);
	}


	if(verbose)
		cout << "chosen solution: " << chosenSolution.score << "(" << chosenSolution.coevent << "): components:" <<  chosenSolution.components.size() << " B: " << chosenSolution.NbBreak  << " G: " << chosenSolution.NbGain << endl;


	for (vector<AdjScore>::iterator it = chosenSolution.components.begin() ; it != chosenSolution.components.end(); ++it) //iterates over the AdjScore of the chosen solution
	{
		if(it->C1) // this one is C1 -> create a new AdjTree
		{
			AdjTree * ATree = new AdjTree(true); //new adjacency tree with a gain at its root
			AdjacencyTrees->push_back(ATree);
			ATree->setRootNode(backtrackAuxC1(it->id1, it->id2, AdjacencyTrees, stochastic));
			
		}
		else // C0 -> continue backtrack. if another c1 is encountered it will be part of another AdjTre anyway
			backtrackAuxC0(it->id1, it->id2, AdjacencyTrees, stochastic);


	}

}

/*
Takes:
	- NodeId1 (int): id of a node in Rtree1
	- NodeId2 (int): id of a node in Rtree2

Returns:
	(pair <int,int>): the event code and the species for the adjacency
*/
pair <int,int> AdjMatrix::getAdjEventAndSpecies(int NodeId1, int NodeId2)
{

	int e1 = Rtree1.getNodeEvent(NodeId1);
	int e2 = Rtree2.getNodeEvent(NodeId2);
	int s1 = Rtree1.getNodeSpecies(NodeId1);
	int s2 = Rtree2.getNodeSpecies(NodeId2);

	int ea,sa;

	//determining adjacency event
	if(e1 == e2)
		ea = e1;
	else //use event hierarchy (same as in the formulae)
	{
		if(Rtree1.isLoss(e1))
			ea = e1;
		else if(Rtree2.isLoss(e2))
			ea = e2;
		else if(Rtree1.isDup(e1))
			ea = e1;
		else if(Rtree2.isDup(e2))
			ea = e2;
		else if(Rtree1.isSout(e1))
			ea = e1;
		else if(Rtree2.isSout(e2))
			ea = e2;
		else if(Rtree1.isBout(e1))
			ea = e1;
		else if(Rtree2.isBout(e2))
			ea = e2;
		else if(Rtree1.isRec(e1))
			ea = e1;
		else if(Rtree2.isRec(e2))
			ea = e2;
		else if(Rtree1.isNull(e1))
			ea = e1;
		else if(Rtree2.isNull(e2))
			ea = e2;
		else if(Rtree1.isSpeciation(e1))
			ea = e1;
		else
			ea = e2;
	}

	//determining adjacency species
	sa = min(s1,s2); //s1 should equal s2 except in the case of tranfers, in which case the node with species -1 (dead) is attributed

	return pair <int,int> (ea,sa);
}

///// Boltzmann specific function /////
void AdjMatrix::setComputationToBoltzmann()
{
	scoreAggregatorfunc = &AdjMatrix::multiplication;
	scoreComparatorfunc = &AdjMatrix::getsum;

	defaultScore = -1;
	worstScore = 0;
	bestScore = 1;

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////public methods////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
Constructor

Takes:
 - map<int,vector<float> > speciesC0C1 : used by artdeco -> cost of having an adjacency or not for different species when the adjacency was not in the input
 - map<int, map<string,int> > speGeneAdjNb : used by artdeco -> degree of each gene, used to modulate the costs in speciesC0C1.
 - vector <double> adjacencyScoreVec : used by ADseq -> score between 0 and 1 of the adjacencies that denote the confidence we have in them
 - Gcost (double): cost of a Gain
 - Bcost (double): cost of a Break
 - rtree1 (ReconciledTree *): reconciled tree for the first dimension
 - rtree2 (ReconciledTree *): reconciled tree for the second dimension
 - adjacencies (vector< pair <string,string> >)
 - VERBOSE (bool)
 - boltzmann (bool) (default: false): wether to use boltzmann computation or not
 - temp (double) (default: 1) : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - absencePenalty (double) (default: -1): if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization. 
 - double adjScoreLogBase [ default : 10000 ]: base of the log that will be used to go from adjacency confidence score to parsimony costs
*/
AdjMatrix::AdjMatrix(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, 
					vector <double> adjacencyScoreVec, 
					double Gcost, double Bcost, 
					ReconciledTree * rtree1, ReconciledTree * rtree2, 
					vector< pair <string,string> > adjacencies, bool VERBOSE, 
					bool boltzmann , double temp , double absencePenalty , double adjScoreLogBase) : Rtree1( *rtree1), Rtree2( *rtree2)//: Rtree1(rtree1->getRootNode(), rtree1->getTimeSliceStatus()), Rtree2(rtree2->getRootNode(), rtree2->getTimeSliceStatus())
{

	//cout << "AdjMatrix::AdjMatrix " << "speciesC0C1: " <<  speciesC0C1.size() << " speGeneAdjNb: " << speGeneAdjNb.size() << endl;

	int nbadj = adjacencies.size();
	
	//cout << "nbadj " << nbadj << endl;
	//cout << "nbadjScores " << adjacencyScoreVec.size() << endl;
	//cout << rtree1->getNumberOfNodes() << " " << rtree2->getNumberOfNodes() << endl;
	//cout << rtree1->getRootId() << " " << rtree2->getRootId() << endl;

	//cout << rtree1->NewickString() << endl;
	//cout << rtree2->NewickString() << endl;

	vector <pair <int, int> > adjs;


	for(unsigned i = 0; i < nbadj ; i++)
	{
		pair <string,string> ps = adjacencies[i];

		//cout << ps.first << "-"<< ps.second <<endl;
	}

	for(unsigned i = 0; i < nbadj ; i++)
	{
		pair <int, int> pi;
		pair <string,string> ps = adjacencies[i];

		//cout << ps.first << "-"<< ps.second <<endl;

		pi.first = rtree1->getIdWithName(ps.first);
		pi.second = rtree2->getIdWithName(ps.second);



		adjs.push_back(pi);

		//cout << ps.first << "-"<< ps.second << " -> "<< pi.first << "-"<< pi.second << endl;
	}


	AdjMatrixAux(speciesC0C1, speGeneAdjNb, adjacencyScoreVec,Gcost, Bcost, rtree1, rtree2, adjs, VERBOSE, boltzmann, temp, absencePenalty, adjScoreLogBase);
}


/*
Constructor

Takes:
 - map<int,vector<float> > speciesC0C1 : used by artdeco -> cost of having an adjacency or not for different species when the adjacency was not in the input
 - map<int, map<string,int> > speGeneAdjNb : used by artdeco -> degree of each gene, used to modulate the costs in speciesC0C1.
 - vector <double> adjacencyScoreVec : used by ADseq -> score between 0 and 1 of the adjacencies that denote the confidence we have in them
 - Gcost (double): cost of a Gain
 - Bcost (double): cost of a Break
 - rtree1 (ReconciledTree *): reconciled tree for the first dimension
 - rtree2 (ReconciledTree *): reconciled tree for the second dimension
 - adjacencies (vector< pair <int,int> >)
 - VERBOSE (bool)
 - boltzmann (bool) (default: false): wether to use boltzmann computation or not
 - temp (double) (default: 1) : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - absencePenalty (double) (default: -1): if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization.
 - double adjScoreLogBase [ default : 10000 ]: base of the log that will be used to go from adjacency confidence score to parsimony costs
*/
AdjMatrix::AdjMatrix(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, 
					vector <double> adjacencyScoreVec, 
					double Gcost, double Bcost, 
					ReconciledTree * rtree1, ReconciledTree * rtree2, vector< pair <int,int> > adjacencies, 
					bool VERBOSE, bool boltzmann , double temp, double absencePenalty , 
					double adjScoreLogBase)//: Rtree1(rtree1->getRootNode(), rtree1->getTimeSliceStatus() ), Rtree2(rtree2->getRootNode(), rtree2->getTimeSliceStatus() )
{
	//cout << Rtree1.getNumberOfNodes() << "<>"<< rtree1->getNumberOfNodes() <<  " - " << Rtree2.getNumberOfNodes() << "<>"<< rtree2->getNumberOfNodes()<< endl;

	AdjMatrixAux(speciesC0C1, speGeneAdjNb, adjacencyScoreVec, Gcost,Bcost, rtree1, rtree2,  adjacencies, VERBOSE, boltzmann,temp, absencePenalty, adjScoreLogBase);
}




/*
Sets all the values of the matrix to -1 (including leaves)
*/
void AdjMatrix::resetMatrix()
{
	matrixComputed = false;
	initMatrix();

}

/*
Reset to the default value the case gievn as argument along with all its ancestors

Takes:
 - id1 (int): a node id in the first reconciled tree
 - id2 (int): a node id in the second reconciled tree


*/
void AdjMatrix::partiallyResetMatrix(int id1, int id2)
{

	vector <int> path1 = Rtree1.idPathFromAncestorToSon(Rtree1.getRootId(), id1);
	vector <int> path2 = Rtree2.idPathFromAncestorToSon(Rtree2.getRootId(), id2);

	for(unsigned i = 0; i < path1.size(); i++)
		for(unsigned j = 0; j < path2.size(); j++)
			resetCase(path1[i],path2[j]);

	matrixComputed = false;
}

/*
Takes:
 - id1 (int): id in reconciled tree 1
 - id2 (int): id in reconciled tree 2
 - invert (bool): wether id1 is really in tree 1 (true) or is in fact in tree 2 (false) (and vice versa for id2)
Returns:
	(double): value of the case in C1

*/
double AdjMatrix::getC1(int id1, int id2, bool invert)
{
	int Id1,Id2;

	if(!invert)
	{
		Id1 = TreeToMatrixId1[id1];
		Id2 = TreeToMatrixId2[id2];
	}
	else
	{
		Id1 = TreeToMatrixId1[id2];
		Id2 = TreeToMatrixId2[id1];
	}

	return MatrixC1.getValue(Id1,Id2);
}

/*
Takes:
 - id1 (int): id in reconciled tree 1
 - id2 (int): id in reconciled tree 2
 - invert (bool): wether id1 is really in tree 1 (true) or is in fact in tree 2 (false) (and vice versa for id2)
Returns:
	(double): value of the case in C0

*/
double AdjMatrix::getC0(int id1, int id2, bool invert)
{
	int Id1,Id2;

	if(!invert)
	{
		Id1 = TreeToMatrixId1[id1];
		Id2 = TreeToMatrixId2[id2];
	}
	else
	{
		Id1 = TreeToMatrixId1[id2];
		Id2 = TreeToMatrixId2[id1];
	}

	return MatrixC0.getValue(Id1,Id2);
}


/*
Takes:
 - id1 (int): id in reconciled tree 1
 - id2 (int): id in reconciled tree 2

Returns:
	(vector <AdjSolution> ): the solution vector of C1( id1,id2 )
*/
vector <AdjSolution>  AdjMatrix::getSolutionC1(int id1, int id2)
{
	int Id1 = TreeToMatrixId1[id1];
	int Id2 = TreeToMatrixId2[id2];
	return SolutionMatrixC1.getValue(Id1,Id2);
}


/*
Takes:
 - id1 (int): id in reconciled tree 1
 - id2 (int): id in reconciled tree 2

Returns:
	(vector <AdjSolution> ): the solution vector of C0( id1,id2 )
*/
vector <AdjSolution>  AdjMatrix::getSolutionC0(int id1, int id2)
{
	int Id1 = TreeToMatrixId1[id1];
	int Id2 = TreeToMatrixId2[id2];

	return SolutionMatrixC0.getValue(Id1,Id2);
}

void AdjMatrix::setSolutionC1(int id1, int id2, vector <AdjSolution> &Vsolution)
{
	int Id1 = TreeToMatrixId1[id1];
	int Id2 = TreeToMatrixId2[id2];

	SolutionMatrixC1.setValues(Id1,Id2,Vsolution);
}
void AdjMatrix::setSolutionC0(int id1, int id2, vector <AdjSolution> &Vsolution)
{
	int Id1 = TreeToMatrixId1[id1];
	int Id2 = TreeToMatrixId2[id2];

	SolutionMatrixC0.setValues(Id1,Id2,Vsolution);
}

/*
Takes:
 - id1 (int): id in reconciled tree 1
 - id2 (int): id in reconciled tree 2

Returns:
	(bool): true if the case is set to something else than the default value in C1
*/
bool AdjMatrix::issetC1(int id1, int id2)
{
	return(getC1(id1,id2) != defaultScore);
}


/*
Takes:
 - id1 (int): id in reconciled tree 1
 - id2 (int): id in reconciled tree 2

Returns:
	(bool): true if the case is set to something else than the default value in C0
*/
bool AdjMatrix::issetC0(int id1, int id2)
{
	return(getC0(id1,id2) != defaultScore);
}

void AdjMatrix::setdecoLTalgo(bool decolt)
{
	decoLTalgo = decolt;
}

double AdjMatrix::getGainCost()
{
	return GainCost;
}


double AdjMatrix::getBreakCost()
{
	return BreakCost;
}


double AdjMatrix::getBoltzmannGainBreakCost(int nbGain, int nbBreak)
{
	return exp( - (nbGain * getGainCost() + nbBreak * getBreakCost()) / (BOLTZMANN_K * Temperature)) ; 
}


void AdjMatrix::setRtree1(ReconciledTree * rtree)
{
	Rtree1 = *rtree;
}

void AdjMatrix::setRtree2(ReconciledTree * rtree)
{
	Rtree2 = *rtree;
}

/*
Set the values of C1 and C0 at the given cases of the matrix

Takes:
 - adjacencies (vector< pair <int,int> >): vector of pairs of ids in reconciled trees 1 and 2
 - vector <double> adjacencyScoreVec : used by ADseq -> score between 0 and 1 of the adjacencies that denote the confidence we have in them
 - double adjScoreLogBase : base of the log that will be used to go from adjacency confidence score to parsimony costs

CAUTION: does not check wether these id exists or are leaves.
CAUTION2: won't do anything if the matrix is computed

*/
void AdjMatrix::addAdjacencies(vector < pair <int,int> > adjacencies, vector <double> adjacencyScoreVec, double adjScoreLogBase)
{

	for(unsigned i = 0; i < adjacencies.size(); i++)
	{
		if(adjacencyScoreVec[i] == 1) // default for 1 -> absolute confidence
			addAdjacency(adjacencies[i]);
		else
			addAdjacency(adjacencies[i], adjacencyScoreVec[i], adjScoreLogBase);
	}

	return;

}

void AdjMatrix::computeMatrix()
{
	if(matrixComputed)
	{
		if(verbose)
			cout << "Matrix already computed. Does nothing." << endl;
		return;
	}

	vector <int> Ids1 = Rtree1.getNodesId();
	vector <int> Ids2 = Rtree2.getNodesId();

	pair <double, double> scores;

	
	for(unsigned i = 0; i < Ids1.size() ; i++)
	{
		int id1 = Ids1[i];

		for(unsigned j = 0; j < Ids2.size() ; j++)
		{
			
			int id2 = Ids2[j];

//			if(verbose)
//			cout << id1 << " , " << id2 << endl;


			if(issetC1(id1,id2)) //the test could also be done on C0 
			{
				continue;
			}
			vector <AdjSolution> Vsolution1;
			vector <AdjSolution> Vsolution0;
			
			computeSolution( id1,id2, Vsolution1, Vsolution0);

			setSolutionC1(id1,id2, Vsolution1);
			setSolutionC0(id1,id2, Vsolution0);

			scores.first  = compareScore(getSolutionC1(id1,id2));
			scores.second = compareScore(getSolutionC0(id1,id2));
			//scores = computeScore(id1,id2); // only works on score and not on AdjSolution

			if(verbose)
			{
				cout << id1 << " , " << id2 << " -> C1: " << scores.first << " C0: " << scores.second << endl;//" default: " << defaultScore << endl; 
			}//plop

			if((scores.first <= defaultScore ) || (scores.second <= defaultScore )) //defaultScore is -1 or -infinity
				throw Exception("AdjMatrix::computeMatrix : used non-set score!");

			//if( useBoltzmann ) //NonSense.
			//{
			//	if( ( scores.first > bestScore ) || ( scores.second > bestScore ) )
			//	{
			//		cout << id1 << " , " << id2 << " -> C1: " << scores.first << " C0: " << scores.second << endl;//" default: " << defaultScore << endl; 
			//		throw Exception("AdjMatrix::computeMatrix : at least one score is > 1 ; Overflow?");
			//	}
			//}


			setC1(id1,id2, scores.first);
			setC0(id1,id2, scores.second);


		}
	}
	if(verbose)
		cout << "Matrix computed." << endl;
	matrixComputed = true;
}


/*
Version of matrix computing that substract reconciliation event scores to co-events in the adjacency matrix

Takes:
 - WDupCost (double): Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - WLossCost (double):Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - WHgtCOst (double): Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
*/
void AdjMatrix::computeMatrix(double WDupCost, double WLossCost, double WHgtCost)
{
	WeightedDupCost  =  - WDupCost;
	WeightedLossCost =  - WLossCost;
	WeightedHgtCost  =  - WHgtCost;

	if(useBoltzmann)
	{
		WeightedDupCost = exp( - (WeightedDupCost) / (BOLTZMANN_K * Temperature)) ; 
		WeightedLossCost = exp( - (WeightedLossCost) / (BOLTZMANN_K * Temperature)) ; 
		WeightedHgtCost = exp( - (WeightedHgtCost) / (BOLTZMANN_K * Temperature)) ; 
	}
	computeMatrix();
}



void AdjMatrix::printC0()
{
	vector <int> Ids1 = Rtree1.getNodesId();
	vector <int> Ids2 = Rtree2.getNodesId();

	for(unsigned j = 0; j < Ids2.size() ; j++)
	{
		int id2 = Ids2[j];
		cout << "\t" << id2;
	}
	cout << endl;

	
	for(unsigned i = 0; i < Ids1.size() ; i++)
	{
		int id1 = Ids1[i];

		cout << id1;

		for(unsigned j = 0; j < Ids2.size() ; j++)
		{
			int id2 = Ids2[j];

			cout << "\t";

			double s = getC0(id1, id2);
			if((s == worstScore) && (!useBoltzmann))
				cout << "inf";
			else
				cout << s;


		}
		cout << endl;
	}
}

void AdjMatrix::printC1()
{
	vector <int> Ids1 = Rtree1.getNodesId();
	vector <int> Ids2 = Rtree2.getNodesId();

	for(unsigned j = 0; j < Ids2.size() ; j++)
	{
		int id2 = Ids2[j];
		cout << "\t" << id2;
	}
	cout << endl;

	
	for(unsigned i = 0; i < Ids1.size() ; i++)
	{
		int id1 = Ids1[i];

		cout << id1;

		for(unsigned j = 0; j < Ids2.size() ; j++)
		{
			int id2 = Ids2[j];

			cout << "\t";

			double s = getC1(id1, id2);
			if((s == worstScore) && (!useBoltzmann))
				cout << "inf";
			else
				cout << s;

		}
		cout << endl;
	}
}

void AdjMatrix::printMe()
{
	cout << "C1:" << endl;
	printC1();

//	cout << "SolutionMatrixC1:" << endl;
//	SolutionMatrixC1.printMe();

	cout << "C0:" << endl;
	printC0();

//	cout << "SolutionMatrixC0:" << endl;
//	SolutionMatrixC0.printMe();

}

/*
Takes:
	- AdjacencyTrees (vector< AdjTree *> *): pointer to a vector of adjacency trees vector that will be update as we build more adjacency trees
 	- stochastic (bool): true if the backtrack is to be stochastic; false if the backtrack choses the solution with the best score
 	- alwaysGainAtTop (bool) [default: true]: there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 	- c1proba (double) [default = 0.5]: probability to choose c1 over c0 IF (and only if) they have the same score
*/
void AdjMatrix::backtrack( vector< AdjTree *> * AdjacencyTrees, bool stochastic, bool alwaysGainAtTop, double c1proba )
{
	if(rootList.size() == 0)
		rootList = getMatrixRoots();
	if(verbose)
		cout << "MultiRootAdjMatrix::backtrack . nbroots : "<< rootList.size() << endl;

	if(rootList.size() == 0) // no root found
		rootList.push_back(pair <int,int> ( Rtree1.getRootId() , Rtree2.getRootId() ) ); // putting root on top of the trees


	for(unsigned i = 0 ; i < rootList.size(); i++)
	{
		backtrackAux(AdjacencyTrees,stochastic,alwaysGainAtTop,c1proba, rootList[i].first, rootList[i].second );
	}
	/*
	int Root1 = Rtree1.getRootId();
	int Root2 = Rtree2.getRootId();

	double c1 = getC1(Root1,Root2);
	double c0 = getC0(Root1,Root2);

	if(alwaysGainAtTop) // we have to add a gain to c1
	{
		if(useBoltzmann) // for Boltzmann DeCo
		{
			c1 *= getBoltzmannGainBreakCost(1,0); 
		}
		else
		{
			c1 += getGainCost();
		}
	}


	//cout << "c1: " << c1 << ", c0: " << c0 << endl;

	bool c1best;

	if(c1 == worstScore)
		c1best = false;
	else if(c0 == worstScore)
		c1best = true;
	else if (c1 == c0) // both score are equal --> randomly determine which to choose
	{
		double r = ((double) rand()/ RAND_MAX); //random result

		if(r <  c1proba)
			c1best = true;
		else
			c1best = false;
	}
	else if(stochastic)//big score has a bigger proba to be chosen. Designed to be used with Boltzmann versio only.
	{
		double r = (c1 + c0) * ((double) rand()/ RAND_MAX); //random result

		if(r < c1)
			c1best = true;
		else
			c1best = false;
	}
	else // best score is chosen ( best means minimum in classical DeCo, and maximal in BoltzmannDeCo)
	{
		c1best = false;
		if(c0 > c1)
			c1best = true;

		if(useBoltzmann)//reverse -> we want the greater score to win
			c1best = !c1best;
	}


	if(c1best) // creation of a new AdjTree
	{
		AdjTree * ATree = new AdjTree(alwaysGainAtTop); //new adjacency tree with a gain at its root
		AdjacencyTrees->push_back(ATree);
		ATree->setRootNode(backtrackAuxC1(Root1, Root2, AdjacencyTrees, stochastic));
		
	}
	else // C0 -> continue backtrack. if another c1 is encountered it will be part of another AdjTree anyway
		backtrackAuxC0(Root1, Root2, AdjacencyTrees, stochastic);

	for(vector< AdjTree *>::iterator it = AdjacencyTrees->begin(); it != AdjacencyTrees->end(); ++it)
	{
		(*it)->resetNodesId();
	}*/
}


bool AdjMatrix::getRootMatrix(int id1, int id2)
{
	int MatrixId1 = TreeToMatrixId1[id1];
	int MatrixId2 = TreeToMatrixId2[id2];

	double rootStatus = RootMatrix.getValue(MatrixId1,MatrixId2);

	if(rootStatus == 1)
		return true;
	return false;
}

void AdjMatrix::setRootMatrixPossible(int id1, int id2)
{
	int MatrixId1 = TreeToMatrixId1[id1];
	int MatrixId2 = TreeToMatrixId2[id2];
	RootMatrix.setValue(MatrixId1,MatrixId2 ,1);	
}

void AdjMatrix::setRootMatrixNotPossible(int id1, int id2)
{
	int MatrixId1 = TreeToMatrixId1[id1];
	int MatrixId2 = TreeToMatrixId2[id2];
	RootMatrix.setValue(MatrixId1,MatrixId2 ,0);
}


AdjMatrix* AdjMatrix::getClone()
{

	vector < pair <int,int> > tmpAdjs;
	vector <double> tmpAdjScore;

	AdjMatrix * newAmat = new AdjMatrix(speciesC0C1, speGeneAdjNb, tmpAdjScore, GainCost, BreakCost, Rtree1.cloneSubtree(Rtree1.getRootId()), Rtree2.cloneSubtree(Rtree2.getRootId()),
										tmpAdjs, verbose, useBoltzmann, Temperature  , worstAbsenceScore);



	newAmat->setdecoLTalgo(decoLTalgo);


	// setting the C1 and C0 matrix
	vector <int> Ids1 = Rtree1.getNodesId();
	vector <int> Ids2 = Rtree2.getNodesId();


	for(unsigned i = 0; i < Ids1.size() ; i++)
	{
		int id1 = Ids1[i];

		for(unsigned j = 0; j < Ids2.size() ; j++)
		{
			int id2 = Ids2[j];

			if(issetC0(id1, id2))
			{
				double s0 = getC0(id1, id2);
				double s1 = getC1(id1, id2);
				newAmat->setC1(id1,id2, s1);
				newAmat->setC0(id1,id2, s0);
			}
		}
	}

	// setting the solution matrix

	for(unsigned i = 0; i < Ids1.size() ; i++)
	{
		int id1 = Ids1[i];

		for(unsigned j = 0; j < Ids2.size() ; j++)
		{
			int id2 = Ids2[j];

			if(issetC0(id1, id2))
			{
				vector <AdjSolution> V0 = getSolutionC1( id1, id2);
				vector <AdjSolution> V1 = getSolutionC0( id1, id2);
				newAmat->setSolutionC1( id1, id2, V1);
				newAmat->setSolutionC0( id1, id2, V0);
			}
		}
	}


	if(matrixComputed)
		newAmat->setComputed();

	if(rootComputed)
		newAmat->getMatrixRoots();


	return newAmat;
}

/*
Takes:
	- double threshold (default = 200)

Returns:
	(int) the number of score whose absolute log 10 are above the given threshold 
			(in both c1 and c0) 
			(ignoring the cases with worst score (infinite in parcimony))
*/
int AdjMatrix::getNumberScoreWithAbsLog10Above(double threshold )
{
	int nbAbove = 0;

	vector <int> Ids1 = Rtree1.getNodesId();
	vector <int> Ids2 = Rtree2.getNodesId();
	
	for(unsigned i = 0; i < Ids1.size() ; i++)
	{
		int id1 = Ids1[i];
		for(unsigned j = 0; j < Ids2.size() ; j++)
		{
			int id2 = Ids2[j];

			double c1 = getC1(id1,id2);
			double c0 = getC0(id1,id2);

			if((c1 != worstScore)&&(c1 != defaultScore))
			{
				if( abs(log10(c1)) > threshold )
					nbAbove++;
			}
			if((c0 != worstScore)&&(c0 != defaultScore))
			{
				if( abs(log10(c0)) > threshold )
					nbAbove++;
			}
		}
	}
	return nbAbove;
}