
#ifndef __ZKEY__
#define __ZKEY__

#include <cstdio>
#include <cstdlib>

using namespace std;

class ZKey {

public:
    ZKey() {
        FILE *fp = fopen("zobristkey", "rb");
        fread(keys, sizeof(unsigned long long), 1000, fp);
        fclose(fp);
    }

    unsigned long long operator[](int i) const {
        return keys[i];
    }

protected:

    unsigned long long keys[1000];

};

#endif
