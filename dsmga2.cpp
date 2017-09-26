/*
 * dsmga2.cpp
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */


#include <list>
#include <vector>
#include <algorithm>
#include <iterator>

#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"


using namespace std;


DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {

    NEW = false;
    pFreeze = false;

    previousFitnessMean = -INF;
    ell = n_ell;
    nPrev = 0;
    nCurrent = n_nInitial;

    Chromosome::length = ell;
    Chromosome::lengthLong = quotientLong(ell)+1;
    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);
    orig_graph.init(ell);

    bestIndex = 0;
    masks = new list<int>[ell];
    orig_masks = new list<int>[ell];
    orderELL = new int[ell];

    selectionIndex.resize(nCurrent);
    orig_selectionIndex.resize(nCurrent);

    orderN.resize(nCurrent);

    orig_fc = new FastCounting[ell];
    for (int i = 0; i < ell; i++)
        orig_fc[i].init(nCurrent);

    fastCounting = new FastCounting[ell];
    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);

    pHash.clear();
    for (int i=0; i<nCurrent; ++i) {
        Chromosome ch;
        population.push_back(ch);
        population[i].initR();
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
    }

    if (GHC) {
        for (int i=0; i < nCurrent; i++) {
            population[i].GHC();
            orig_popu.push_back(population[i]);
        }
    }

    for (int i=0; i<MAXLEVEL; ++i) {
        nIndex.push_back(vector<int>(0));
        BMhistory.push_back(vector<BMRecord>(0));
    }
}


DSMGA2::~DSMGA2 () {
    delete []masks;
    delete []orig_masks;
    delete []orderELL;
    delete []fastCounting;
    delete []orig_fc;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + 1e-6;
        return false;
    }

    return true;
}



int DSMGA2::doIt (bool output) {
    generation = 0;

    while (!shouldTerminate ()) {
        oneRun (output);
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {

    /*
    if (CACHE)
        Chromosome::cache.clear();
        */

    mixing();


    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

        if (SHOW_POPULATION) {
            if (SHORT_HAND)
                population[i].printOut();
            else
                population[i].shortPrintOut();
            cout << endl;
        }
    }

    if (output)
        showStatistics ();

    ++generation;


}


bool DSMGA2::shouldTerminate () {
    bool
    termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }


    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    /*
    if (stFitness.getMax() - 1e-10 <= stFitness.getMean() )
        termination = true;
        */

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen:%d  N:%d  Fitness:(Max/Mean/Min):%f/%f/%f ",
            generation, nCurrent, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    printf ("best chromosome:\n");
    population[bestIndex].printOut();
    printf ("\n");


    fflush(NULL);
}


void DSMGA2::buildOrigFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                orig_fc[j].setVal(i, orig_popu[orig_selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++) {
            for (int j = 0; j < ell; j++) {
                orig_fc[j].setVal(i, orig_popu[i].getVal(j));
            }
        }
    }

}

void DSMGA2::buildFastCounting(int level) {

    if (SELECTION) {
        int counter = 0;
        for (int i=0; i<nCurrent; ++i)
            if (population[selectionIndex[i]].level == level)
                ++counter;

        for (int i = 0; i < ell; i++)
            fastCounting[i].init(counter);

        counter = 0;
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                if (population[selectionIndex[i]].level == level)
                    fastCounting[j].setVal(counter, population[selectionIndex[i]].getVal(j));
                ++counter;
            }

    } else {

        for (int i = 0; i < ell; i++)
            fastCounting[i].init(nIndex[level].size());

        int counter = 0;
        // for (int index:nIndex[level]) {
        for (int i = 0; i < nIndex[level].size(); ++i) {
            int index = nIndex[level][i];
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(counter, population[index].getVal(j));

            }
            ++counter;

        }
    }

}

int DSMGA2::countOrigOne(int x) const {

    int n = 0;

    for (int i=0; i<orig_fc[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = orig_fc[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countOrigXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<orig_fc[0].lengthLong; ++i) {

        unsigned long val = 0;


        val = orig_fc[x].gene[i];

        val ^= orig_fc[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}

// Do BM
int DSMGA2::restrictedMixing(Chromosome& ch, int pos) {

    list<int> mask = masks[pos];

    size_t size;
    size = findSize(ch, mask, ch.level);
    //size = findSize(ch, mask);

    if (size > (size_t)ell/2)
        size = ell/2;

    // prune mask to exactly size
    while (mask.size() > size)
        mask.pop_back();


    Chromosome copy = ch;
    int resultRM = restrictedMixing(copy, mask);

    EQ = true;
    if (resultRM !=0) {

        // BM to the next level
        // for (auto index:nIndex[ch.level+1]) {
        for (int i = 0; i < nIndex[ch.level+1].size(); ++i) {
            int index = nIndex[ch.level+1][i];

            if (EQ)
                backMixingE(copy, mask, population[index]);
            else
                backMixing(copy, mask, population[index]);
        }

        // BM to the current level
        // for (auto index:nIndex[ch.level]) {
        for (int i = 0; i < nIndex[ch.level].size(); ++i) {
            int index = nIndex[ch.level][i];

            if (EQ)
                backMixingE(copy, mask, population[index]);
            else
                backMixing(copy, mask, population[index]);
        }

        BMhistory[ch.level+1].push_back(BMRecord(copy, mask, EQ, 0.0));
    }

    if (resultRM == 2)  {
        // Add copy to the next level
        copy.level = ch.level+1;
        int p_size = population.size();
        increaseOne(copy);
        if (p_size == population.size())
            ++resultRM;
    }

    return resultRM;
}

int DSMGA2::restrictedMixing(Chromosome& ch) {

    int r = myRand.uniformInt(0, ell-1);
    return restrictedMixing(ch, r);

}

void DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial = des;

    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    Chromosome& real = des;

    if (RTRBM) {
        int minDist = trial.getHammingDistance(des);
        for (int i=0; i<RTRW; ++i) {
            int r = myRand.uniformInt(0, nCurrent-1);
            int dist = trial.getHammingDistance(population[r]);
            if (minDist > dist) {
                minDist = dist;
                real = population[r];
            }
        }
    }

    if (trial.getFitness() > real.getFitness()) {
        if (!pFreeze) {
            pHash.erase(real.getKey());
            pHash[trial.getKey()] = trial.getFitness();
        }
        real = trial;
        return;
    }

}

void DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial = des;

    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    Chromosome& real = des;

    if (RTRBM) {
        int minDist = trial.getHammingDistance(des);
        for (int i=0; i<RTRW; ++i) {
            int r = myRand.uniformInt(0, nCurrent-1);
            int dist = trial.getHammingDistance(population[r]);
            if (minDist > dist) {
                minDist = dist;
                real = population[r];
            }
        }
    }

    if (trial.getFitness() > real.getFitness()) {
        if (!pFreeze) {
            pHash.erase(real.getKey());
            pHash[trial.getKey()] = trial.getFitness();
        }

        EQ = false;
        real = trial;
        return;
    }

    if (trial.getFitness() >= real.getFitness()) {
        if (!pFreeze) {
            pHash.erase(real.getKey());
            pHash[trial.getKey()] = trial.getFitness();
        }

        real = trial;
        return;
    }

}

// This does NOT do BM
int DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    int result = 0;
    size_t lastUB = 0;

    for (size_t ub = 1; ub <= mask.size(); ++ub) {

        size_t size = 1;
        Chromosome trial = ch;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }

        if (isInP(trial)) break;
        //continue;

        double fff = trial.getFitness();


        if (fff >= ch.getFitness()) {
/*
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();
*/
            result = (fff == ch.getFitness())? 1:2;

            ch = trial;
        }

        if (result !=0) {
            lastUB = ub;
            break;
        }
    }

    // prune mask for backmixing
    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return result;

}

size_t DSMGA2::findOrigSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (orig_popu[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}


size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, int level) const {

    DLLA candidate(nCurrent);
    // for (auto i:nIndex[level])
    //     candidate.insert(i);

    for (int i = 0; i < nIndex[level].size(); ++i)
        candidate.insert(nIndex[level][i]);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, Chromosome& ch2) const {

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::mixing() {


    increaseOne();
    int resultRM = 0;
    int level = 0;
    do {

        if (SELECTION)
            selection(level);
        buildFastCounting(level);
        buildGraph();

        for (int i=0; i<ell; ++i)
            findClique(i, masks[i]);

        genOrderELL();
        for (int i=0; i<ell; ++i) {
            int pos = orderELL[i];
            resultRM = restrictedMixing(population[nCurrent-1], pos);
            if (resultRM >= 2) break;
        }

        ++level;
    } while (resultRM == 2);


}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildOrigGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOrigOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countOrigXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;

            double linkage;
            linkage = computeMI(p00,p01,p10,p11);
            orig_graph.write(i,j,linkage);

        }
    }

    if (SHOW_GRAPH) {
        for (int i=0; i<ell; ++i) {
            for (int j=0; j<ell; ++j)
                printf("%3f ", graph(i,j));
            printf("\n");
        }
    }

    delete []one;

}
void DSMGA2::buildGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;

            double linkage;
            linkage = computeMI(p00,p01,p10,p11);
            graph.write(i,j,linkage);
        }
    }

    if (SHOW_GRAPH) {
        for (int i=0; i<ell; ++i) {
            for (int j=0; j<ell; ++j)
                printf("%3f ", graph(i,j));
            printf("\n");
        }
    }

    delete []one;

}

// from 1 to ell, pick by max edge
void DSMGA2::findOrigClique(int startNode, list<int>& result) {


    result.clear();

    DLLA rest(ell);
    genOrderELL();
    for (int i=0; i<ell; ++i) {
        if (orderELL[i]==startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }

    double *connection = new double[ell];

    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
        connection[*iter] = orig_graph(startNode, *iter);

    while (!rest.isEmpty()) {

        double max = -INF;
        int index = -1;
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }

        rest.erase(index);
        result.push_back(index);

        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
            connection[*iter] += orig_graph(index, *iter);
    }


    delete []connection;

}

// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, list<int>& result) {


    result.clear();

    DLLA rest(ell);
    genOrderELL();
    for (int i=0; i<ell; ++i) {
        if (orderELL[i]==startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }

    double *connection = new double[ell];

    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
        connection[*iter] = graph(startNode, *iter);

    while (!rest.isEmpty()) {

        double max = -INF;
        int index = -1;
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }

        rest.erase(index);
        result.push_back(index);

        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
            connection[*iter] += graph(index, *iter);
    }


    delete []connection;

}

double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > 1e-10)
        join += a00*log(a00);
    if (a01 > 1e-10)
        join += a01*log(a01);
    if (a10 > 1e-10)
        join += a10*log(a10);
    if (a11 > 1e-10)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > 1e-10)
        p += p0*log(p0);
    if (p1 > 1e-10)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > 1e-10)
        q += q0*log(q0);
    if (q1 > 1e-10)
        q += q1*log(q1);

    return -p-q+join;

}


void DSMGA2::selection (int level) {
    tournamentSelection (level);
}


// tournamentSelection without replacement
void DSMGA2::OrigSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = orig_popu[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        orig_selectionIndex[i] = winner;
    }
}

// tournamentSelection without replacement
void DSMGA2::tournamentSelection (int level) {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}


void DSMGA2::increaseOne (Chromosome& ch) {

    if (isInP(ch)) return;

    // APPLY BMHISTORY at that level first
    pFreeze = true;
    if (BMhistory.size() > ch.level) {
        if (!BMhistory[ch.level].empty()) {
            // for (auto bm: BMhistory[ch.level]) {
            for (int i = 0; i < BMhistory[ch.level].size(); ++i) {
                BMRecord bm = BMhistory[ch.level][i];
                if (bm.eq) {
                    EQ = true;
                    backMixingE(bm.pattern, bm.mask, ch);
                    if (!EQ)
                        bm.eq = false;
                }
                else
                    backMixing(bm.pattern, bm.mask, ch);
            }
        }
    }
    pFreeze = false;

    if (isInP(ch)) return;

    if (ch.level > MAXLEVEL) {
        cout << "Exceed max level" << endl;
        exit(-1);
    }

    // really ADD
    ++nCurrent;
    population.push_back(ch);
    nIndex[ch.level].push_back(nCurrent-1);
    selectionIndex.resize(nCurrent);
    pHash[ch.getKey()] = ch.getFitness();
/****************************
    cout << "Add level " << ch.level << " : ";
    ch.printOut();
    cout << endl;
****************************/
}

void DSMGA2::increaseOne () {

    Chromosome ch;

    int old = nCurrent;
    do {
        do {
            ch.initR();
            ch.GHC();
        } while (isInP(ch));

        increaseOne(ch);
    } while (old == nCurrent);

}

