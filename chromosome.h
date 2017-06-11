/***************************************************************************
 *   Copyright (C) 2011 by TEIL                                        *
 *                                                                         *
 ***************************************************************************/
#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include <unordered_map>
#include "global.h"
#include "nk-wa.h"

using namespace std;

class Chromosome {

public:

    static enum Function {
        ONEMAX=0,
        MKTRAP=1,
        FTRAP=2,
        CYCTRAP=3,
        NK=4,
        SPINGLASS=5,
        SAT=6
    } function;


    Chromosome ();
    Chromosome (const Chromosome& c) {
	gene = NULL;
	level = 0;
        init();
        *this = c;
    }

    ~Chromosome ();

    bool hasSeen() const;

    bool GHC();
    void steepestDescent();

    void init ();
    void init0 ();
    void initR ();

    bool tryFlipping (int index);

    int getHammingDistance (const Chromosome& c) const {
        int result = 0;
        for (int i=0; i<lengthLong; ++i)
            result += myBD.getHammingDistance(gene[i], c.gene[i]);
        return result;
    }

    int getVal (int index) const {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        if ( (gene[q] & (1lu << r)) == 0 )
            return 0;
        else
            return 1;
    }

    void setVal (int index, int val) {

        assert (index >= 0 && index < length);

        if (getVal(index) == val) return;

        setValF(index, val);
        key ^= zKey[index];
    }

    unsigned long getKey () const {
        return key;
    }


    void setValF (int index, int val) {

        assert (index >= 0 && index < length);
        //outputErrMsg ("Index overrange in Chromosome::operator[]");

        int q = quotientLong(index);
        int r = remainderLong(index);

        if (val == 1)
            gene[q] |= (1lu<<r);
        else
            gene[q] &= ~(1lu<<r);

        evaluated = false;
    }

    void flip (int index) {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        gene[q] ^= (1lu<<r);
        key ^= zKey[index];

        evaluated = false;
    }

    /** real evaluator */
    double evaluate ();

    bool isEvaluated () const;

    bool operator== (const Chromosome & c) const;
    Chromosome & operator= (const Chromosome & c);

    double getFitness ();
    double trap (int u, double high, double low, int trapK) const;
    double oneMax () const;
    double mkTrap (double high, double low) const;
    double cycTrap(double fHigh, double fLow) const;
    double fTrap () const;
    double spinGlass () const;
    double nkFitness() const;
    double satFitness() const;


    void printOut () const;
    void shortPrintOut () const;

    int getLength () const;
    void setLength ();

    double getMaxFitness () const;


public:
    int level;
    static int length;
    static int lengthLong;
    static int nfe;
    static int lsnfe;
    static int hitnfe;
    static bool hit;
    static unordered_map<unsigned long, double> cache;

    bool killed;

protected:

    unsigned long *gene;
    double fitness;
    bool evaluated;
    unsigned long long key;

};


#endif
