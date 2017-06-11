
#ifndef _BMRECORD_H_
#define _BMRECORD_H_

#include <list>
#include "chromosome.h"


class BMRecord {

public:
    BMRecord() {;}
    BMRecord(const Chromosome& _pattern, const list<int>& _mask, bool _eq, double _rate) {
        pattern = _pattern;
        mask = _mask;
        eq = _eq;
        rate = _rate;
    }

    bool operator< (const BMRecord& rec) const {
        //return mask.size() < rec.mask.size();
        return rate > rec.rate;
    } 

public:
    Chromosome pattern;
    list<int> mask;
    bool eq;
    double rate;
};

#endif
