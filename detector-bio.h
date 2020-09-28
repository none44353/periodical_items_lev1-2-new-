#ifndef _detector-bio_H
#define _detector-bio_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include <queue>
#include "BOBHash32.h"
#include "BOBHash64.h"

using namespace std;

#define INF_INT 1e9
#define T 20
//#define MAXM 5000005
#define MAXL 10000005 

typedef pair <int, uint64_t> pii;

struct node { //16Byte
    node* nxt;
    int pos;
    uint64_t id;

    node(uint64_t id, int pos) : id(id), pos(pos) {nxt = NULL; };
};
/*
class heap { //(12+16)*cap + 8 + 4  28012B
  public:
    int sze, cap;
    pair <int, uint64_t> a[MAXL];

    node *p[MAXL];
    BOBHash32* hash;

    heap(int cap) : cap(cap) {
        sze = 0;
        for (int i = 0; i < cap; ++i) p[i] = NULL;
        hash = new BOBHash32(456);
    };

    ~heap() {
        delete hash;
        for (int i = 0; i < cap; ++i) 
            while (p[i] != NULL) {
                node* np = p[i] -> nxt;
                delete p[i];
                p[i] = np;
            }

        return;
    }

    inline int size() {
        return sze;
    }

    int query(const uint64_t& x) { //询问x在不在堆中, 如果在，返回具体位置
        unsigned long long int H = hash -> run((char *)&x, 8);
        unsigned long long int pos = H % cap;
        
        node* np = p[pos];
        while (np != NULL) {
            if (np -> id == x) return np -> pos;
            np = np -> nxt;
        }
        return 0;
    }

    void hash_change(const uint64_t& x, int y) { //修改hash表里的位置信息
        unsigned long long int H = hash -> run((char *)&x, 8);
        unsigned long long int pos = H % cap;

        node* np = p[pos];
        while (np != NULL) {
            if (np -> id == x) {
                np -> pos = y;
                return;
            }
            np = np -> nxt;
        }
        return;
    }

    void hash_del(const uint64_t& x) { //删除节点x在hash表中的信息(如果x在堆中的话)
        unsigned long long int H = hash -> run((char *)&x, 8);
        unsigned long long int pos = H % cap;
        
        if (p[pos] == NULL) return;

        node* lst = NULL; 
        node* np = p[pos];
        while (np != NULL) {
            if (np -> id == x) {
                if (lst != NULL) lst -> nxt = np -> nxt;
                else p[pos] = np -> nxt;

                delete np;
                return;
            }
            lst = np; np = np -> nxt;
        }
        return;
    }

    void hash_add(const uint64_t& x, int y) { //在hash表中加入节点x的信息
        unsigned long long int H = hash -> run((char *)&x, 8);
        unsigned long long int pos = H % cap;
        
        node* np = new node(x, y);
        np -> nxt = p[pos];
        p[pos] = np;

        return;
    }

    void inv(const uint64_t& x, int pos) { //x的counter+1, x原来在堆中
        a[pos].first++;

        int id = pos;
        while (true) {
            int idl = id << 1, idr = (id << 1) | 1, c = id;
            if (idl <= sze && a[idl] < a[c]) c = idl;
            if (idr <= sze && a[idr] < a[c]) c = idr;

            if (c != id) {
                swap(a[c], a[id]);
                hash_change(a[id].second, id);
                id = c;
            }
            else break;
        }
        hash_change(x, id);
        return;
    }

    void add(const uint64_t& x, int y) { //加一个x在末尾，前提是sze < cap
        a[++sze].first = y, a[sze].second = x;
        
        int id = sze;
        while (true) {
            int c = id >> 1;
            if (c > 0 && a[id] < a[c]) {
                swap(a[c], a[id]);
                hash_change(a[id].second, id);
                id = c;
            } 
            else break;
        }
        hash_add(x, id);
        return;
    }

    pair <int, uint64_t> change(const uint64_t& x, int y) { //删去原有的堆顶，放进一个新的节点
        pair <int, uint64_t> res = a[1];
        hash_del(res.second);

        a[1] = make_pair(y, x);
        int id = 1;
        while (true) {
            int idl = id << 1, idr = (id << 1) | 1, c = id;
            if (idl <= sze && a[idl] < a[c]) c = idl;
            if (idr <= sze && a[idr] < a[c]) c = idr;

            if (c != id) {
                swap(a[c], a[id]);
                hash_change(a[id].second, id);
                id = c;
            }
            else break;
        }
        hash_add(x, id);

        return res;
    }

    pair <int, uint64_t> top() { //返回堆顶元素
        return a[1];
    }
};*/

statistic int hash_index = 123;

class Sketch3 {
    private:
        int* counter;
        int nm, nt; //这个sketch的大小nm*nt，hash函数个数为nt
        BOBHash32* bobhash[T]; 
        double L, R, cL, cR;
    
    Sketch3(int nm, int nt, double L, double R) : nm(nm), nt(nt), L(L), R(R) {
        for (int k = 0; k < nt; ++k) bobhash[k] = hash_index++;
        counter = new int[nm * nt * 3];
        memset(counter, 0, nm * nt * 3 * sizeof(int));

        double c = (L + R) / 2;
        cL = c - delta, cR = c + delta;
    }

    int query(const uint64_t& x, const int& type == -1) {
        if (type >= 0) {
            int res = INF_INT;
            for (int k = 0; k < nt; ++k) {
                unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
                unsigned long long int pos = H % M;

                res = min(res, counter[(3 * type + k) * M + pos]);
            }
            return res;
        }
        
        int res_max = 0;
        for (int t = 0; t < 3; ++t) {
            int res = INF_INT;
            for (int k = 0; k < nt; ++k) {
                unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
                unsigned long long int pos = H % M;

                res = min(res, counter[(3 * t + k) * M + pos]);
            }
            
            double len = cR - cL;
            if (t == 0) len = cL - L;
            if (t == 2) len = R - cR;

            res_max = max(res_max, (int)ceil(res * delta * 2 / len));
        }
        
        return res_max;
    }

    void insert(const uint64_t& x, const double& y) { //先暂时直接用CU 不维护堆
        int type = (y < cL ? 0 : (y > cR ? 2 : 1));

        int mn = INF_INT;
        for (int k = 0; k < nt; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            if (counter[(3 * type + k) * M + pos] < mn) mn = counter[(3 * type + k) * M + pos];
        }
       
        //CU
        for (int k = 0; k < nt; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            if (counter[(3 * type + k) * M + pos] == mn) ++counter[(3 * type + k) * M + pos];
        }

        if (type == 0 || type == 2) { //如果出现次数达到一定程度，我需要往下一级的sketch3走

        }
    }
};

class Alg //CM  t * M 
{   
    //堆里存一部分，sketch存一部分

  private:
    int nT, M, nK; //t次hash 共t*M个位置
    BOBHash32* bobhash[T];

    heap* h; //(4 + 8) * t * M Byte
    double a[MAXL];//A部分，一个时间戳
    int counter[MAXL]; //B部分，一个sketch结构 CM / CU
    
  public:
    Alg(int nT, int M, int nK) : nT(nT), M(M), nK(nK) {
        h = new heap(nK);

        int x = 123;
        for (int k = 0; k < nT; ++k) bobhash[k] = new BOBHash32(x + k);

        memset(a, 0, sizeof(a));
        memset(counter, 0, sizeof(counter));
    }

    ~Alg() {
        delete h;
        for (int k = 0; k < nT; ++k) delete bobhash[k];
    }

    bool check_near(double c, double x) {
        return fabs(x - c) <= delta;
    }

  

    void insert(const uint64_t& x, const double& y) {
        double last_time = y + 1; 
        for (int k = 0; k < nT; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            last_time = min(last_time, a[k * M + pos]);
            a[k * M + pos] = y;
        }

        if (y - last_time > rangeR) return;
        if (y - last_time < rangeL) return;

        //接下来发生一次加1
        int pos = h -> query(x);
        if (pos) { h -> inv(x, pos); return;}
        if (h -> size() < nK) {h -> add(x, 1); return;}

        pair <int, uint64_t> top = h -> top();
        for (int k = 0; k < nT; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;

        }

        int mn = INF_INT;
        for (int k = 0; k < nT; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            if (counter[k * M + pos] > top.first) counter[k * M + pos] = top.first; 
                //加之前sketch中的元素不超过堆中的最小元素
            if (counter[k * M + pos] < mn) mn = counter[k * M + pos];
        }

        if (mn + 1 > top.first) {
            h -> change(x, mn + 1);
            for (int k = 0; k < nT; ++k) {
                unsigned long long int H = bobhash[k] -> run((char *)&(top.second), 8);
                unsigned long long int pos = H % M;
                
                if (counter[k * M + pos] < top.first) counter[k * M + pos] = top.first;
            }
        }
        else {
        //CU
            for (int k = 0; k < nT; ++k) {
                unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
                unsigned long long int pos = H % M;
                
                if (counter[k * M + pos] == mn) ++counter[k * M + pos];
            }

       //CM
            /*for (int k = 0; k < nT; ++k) {
                unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
                unsigned long long int pos = H % M;

                if (counter[k * M + pos] < top.first) ++counter[k * M + pos];
            }*/
        }

        return;
    }

    int query(const uint64_t& x) {
        int pos = h -> query(x);
        if (pos) return h -> a[pos].first;

        int res = INF_INT;
        for (int k = 0; k < nT; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;

            res = min(res, counter[k * M + pos]);
        }

        if (h -> size()) {
            pair <int, uint64_t> top;
            res = min(res, top.first);
        }

        return res;
    }

    void getTopK(set <uint64_t>& s) {
        s.clear();
        for (int i = 1, sze = h -> size(); i <= sze; ++i)
            s.insert(h -> a[i].second);
        return;
    }

};
#endif
