#ifndef _detector_H
#define _detector_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include "BOBHash32.h"
#include "BOBHash64.h"

using namespace std;

#define inf 1e9
#define T 20
#define MAXM 500005
#define MAXL 1000005 

class AlgCM //CM  t * M 
{   

  private:
    int t, M; //t次hash 共t*M个位置
    BOBHash32* bobhashA[T];
    BOBHash32* bobhashB[T];

    double a[MAXL];
    uint64_t id[MAXL]; //A部分，一个时间戳

    int counter[MAXL]; //B部分，一个sketch结构 CM / CU
    
  public:
    AlgCM(int t, int M) : t(t), M(M) {
        cerr << "CM + CU" << endl;
        srand(time(NULL));
        int x = 123;
        for (int k = 0; k < t; ++k) bobhashA[k] = new BOBHash32(x + k);
        x = 987;
        for (int k = 0; k < t; ++k) bobhashB[k] = new BOBHash32(x + k);
        memset(a, 0, sizeof(a));
    }

    bool check_near(double c, double x) {
        return fabs(x - c) <= delta;
    }

    void insert(const uint64_t& x, const double& y) {
        double last_time = y + 1; 
        int idt = 0;
        for (int k = 0; k < t; ++k) {
            unsigned long long int H = bobhashA[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            if (a[k * M + pos] < last_time) last_time = a[k * M + pos], idt = id[k * M + pos];
            a[k * M + pos] = y, id[k * M + pos] = x;
        }

        if (y - last_time > rangeR) return;
        if (idt == x && y - last_time < rangeL) return;
        //剩下的情况都算作有一个 高估 

//CU
        int mn = inf;
        for (int k = 0; k < t; ++k) {
            unsigned long long int H = bobhashB[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            if (counter[k * M + pos] < mn) mn = counter[k * M + pos];
        }
        for (int k = 0; k < t; ++k) {
            unsigned long long int H = bobhashB[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            if (counter[k * M + pos] == mn) ++counter[k * M + pos];
        }

//CM
    /*    for (int k = 0; k < t; ++k) {
            unsigned long long int H = bobhashB[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            ++counter[k * M + pos];
        }*/

        return;
    }

    int query(const uint64_t& x) {
        int res = inf;
        for (int k = 0; k < t; ++k) {
            unsigned long long int H = bobhashB[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;

            res = min(res, counter[k * M + pos]);
        }

        return res;
    }

};
#endif
