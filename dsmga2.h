/*
 * dsmga2.h
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */

#ifndef _DSMGA2_H_
#define _DSMGA2_H_

#include <list>
#include "chromosome.h"
#include "statistics.h"
#include "trimatrix.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"
#include "bmrecord.h"

#define MAXLEVEL 10000


class DSMGA2 {
public:
    DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff);

    ~DSMGA2 ();

    void selection (int);
    /** tournament selection without replacement*/
    void tournamentSelection(int);
    void OrigSelection();

    void initialBuildLevels (bool output = true);
    void oneRun (bool output = true);
    int doIt (bool output = true);

    void buildGraph ();
    void buildOrigGraph ();
    void mixing ();
    bool OrigRM(Chromosome&);

    // 0: <
    // 1: ==
    // 2: >
    int restrictedMixing(Chromosome&);
    int restrictedMixing(Chromosome& ch, list<int>& mask);
    int restrictedMixing(Chromosome& ch, int pos, bool init = false);

    void backMixing(Chromosome& source, list<int>& mask, Chromosome& des, bool init = false);
    void backMixingE(Chromosome& source, list<int>& mask, Chromosome& des, bool init = false);

    bool shouldTerminate ();

    bool foundOptima ();

    int getGeneration () const {
        return generation;
    }

    bool isInP(const Chromosome& ) const;
    void genOrderN();
    void genOrderELL();

    void showStatistics ();

    bool isSteadyState ();

//protected:
public:

    bool ADD;
    bool NEW;

    int ell;                                  // chromosome length
    int nCurrent;                             // population size
    int nPrev;                             // population size
    bool EQ;
    bool pFreeze;
    int tempIndex;
    bool initBuilding;
    unordered_map<unsigned long, double> pHash; // to check if a chromosome is in the population

    vector<vector<BMRecord> > BMhistory;
    vector<BMRecord> BMlevel;
    vector<vector<int> > nIndex;

    list<int> *masks;
    list<int> *orig_masks;
    vector<int> selectionIndex;
    vector<int> orig_selectionIndex;
    vector<int> orderN;                             // for random order
    int *orderELL;                             // for random order
    int selectionPressure;
    int maxGen;
    int maxFe;
    int repeat;
    int generation;
    int bestIndex;

    vector<Chromosome> population;
    vector<Chromosome> orig_popu;
    FastCounting* fastCounting;
    FastCounting* orig_fc;

    TriMatrix<double> graph;
    TriMatrix<double> orig_graph;


    double previousFitnessMean;
    Statistics stFitness;

    // methods
    double computeMI(double, double, double, double) const;


    void findClique(int startNode, list<int>& result);
    void buildFastCounting(int);
    int countXOR(int, int) const;
    int countOne(int) const;

    void findOrigClique(int startNode, list<int>& result);
    void buildOrigFastCounting();
    int countOrigXOR(int, int) const;
    int countOrigOne(int) const;

    size_t findSize(Chromosome&, list<int>&) const;
    size_t findSize(Chromosome&, list<int>&, int) const;
    size_t findOrigSize(Chromosome&, list<int>&) const;
    size_t findSize(Chromosome&, list<int>&, Chromosome&) const;

    void increaseOne();
    void increaseOne(Chromosome& );


};


#endif /* _DSMGA2_H_ */
