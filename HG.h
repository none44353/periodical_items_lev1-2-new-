#ifndef _HG_H
#define _HG_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include <math.h>
#include <time.h>
#include <set>
#include "BOBHash32.h"
#include "BOBHash64.h"
#include "BF.h"

using namespace std;

#define INF_INT 1e9
#define MAXM 5000005
#define MAXL 10000005 

template <typename TT>
void AdjustHeap(TT* a, int k, int x){
    int idx = x;
    int lchild = 2*x + 1;
    int rchild = 2*x + 2;
    if(idx < k){
        if(lchild < k && a[idx] > a[lchild]){
            idx = lchild;
        }
        if(rchild < k && a[idx] > a[rchild]){
            idx = rchild;
        }
        if(idx != x){
            TT tmp = a[idx];
            a[idx] = a[x];
            a[x] = tmp;
            AdjustHeap(a, k, idx);
        }
    }
}

class heavyguardian{
  private:
    int w, h, l, bucket;
    BOBHash32* bobhash[2];
    uint64_t a[MAXM];
    double b;
  public:
    heavyguardian(int w, int h, int l, double b) : w(w), h(h), l(l), b(b) {
        bucket = 2 * h + l;
        srand((int)time(NULL));
        int x = 123;
        for (int k = 0; k < 2; ++k) 
            bobhash[k] = new BOBHash32(x + k);
        memset(a, 0, sizeof(a));
    }

    void insert(const uint64_t& x) {
        unsigned long long int H = bobhash[0] -> run((char *)&x, 8);
        unsigned long long int pos = H % w;
        uint64_t res = INF_INT;
        int posi = 0;
        for (int k = 0; k < h; ++k) {
            if (!a[pos * bucket + 2 * k + 1]) {
                a[pos * bucket + 2 * k] = x;
                a[pos * bucket + 2 * k + 1] = 1;
                return;
            } else if (a[pos * bucket + 2 * k] == x) {
                ++a[pos * bucket + 2 * k + 1];
                return;
            } else if (a[pos * bucket + 2 * k + 1] < res){
                res = a[pos * bucket + 2 * k + 1];
                posi = k;
            }
        }
        double e = 1.0 / pow(b, res);
        if (double(rand()) / RAND_MAX < e) {
            --a[pos * bucket + 2 * posi + 1];
            if (!a[pos * bucket + 2 * posi + 1]) {
                a[pos * bucket + 2 * posi] = x;
                a[pos * bucket + 2 * posi + 1] = 1;
            }
        } else if (l){
            unsigned long long int H1 = bobhash[1] -> run((char *)&x, 8);
            unsigned long long int pos1 = H % l; 
            ++a[pos * bucket + 2 * h + pos1];           
        }

        return;
    }

    uint64_t query(const uint64_t& x) {
        unsigned long long int H = bobhash[0] -> run((char *)&x, 8);
        unsigned long long int pos = H % w;
        for (int k = 0; k < h; ++k) {
            if (a[pos * bucket + 2 * k + 1] && a[pos * bucket + 2 * k] == x) {
                return a[pos * bucket + 2 * k + 1];
            }
        }
        return 0;
        unsigned long long int H1 = bobhash[1] -> run((char *)&x, 8);
        unsigned long long int pos1 = H % l; 
    }

    pair<uint64_t,uint64_t>* topk(int k){
        pair<uint64_t,uint64_t> * ret = new pair<uint64_t, uint64_t>[k];
        int cur = 0;
        int curi = 0;
        int curj = 0;
        while (cur < k){
            pair<uint64_t, uint64_t> p;
            p.first = a[curi * bucket + 2 * curj + 1];
            p.second = a[curi * bucket + 2 * curj];
            ret[cur] = p;
            ++cur;
            ++curj;
            if (curj == h) {
                ++curi;
                if (curi == w) {
                    return ret;
                }
                curj = 0;
            }
        }
        for (int i=(k-1)/2; i>0; --i){
            AdjustHeap(ret, k, i);
        }
        while(1) {
            pair<uint64_t, uint64_t> p;
            p.first = a[curi * bucket + 2 * curj + 1];
            p.second = a[curi * bucket + 2 * curj];
            if (ret[0] < p) {
                ret[0] = p;
                AdjustHeap(ret, k, 0);
            }
            ++cur;
            ++curj;
            if (curj == h) {
                ++curi;
                if (curi == w) {
                    return ret;
                }
                curj = 0;
            }            
        }
    }
};

#endif