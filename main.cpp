#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <time.h>

using namespace std;

#define Range 500
#define Centre 0.03
#define delta 0.01
#define rangeL (Centre - delta)
#define rangeR (Centre + delta)


#include "ssummary.h"
#include "BF.h"
#include "detector2.h"
#include "HG.h"
//#include "detectorCM.h"

//string datapath[60] = {"./130000.dat"};

string datapath[60] = {"../../usr/share/dataset/CAIDA2018/dataset/130000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135600.dat", 
                "../../usr/share/dataset/CAIDA2018/dataset/135700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135900.dat"};

ifstream fin;

const int M = 2e7;
const int _K = 1000;

double firstTime;

pair <uint64_t, double> Read()//新CAIDA
{   
    static bool isfirstRead = true;
    static int curFinID = 0;
    static double offset = 0;
    static double lastT = 0;

	double t; uint64_t s, _s;
    if (isfirstRead) {
        isfirstRead = false;
        firstTime = -1;
        fin.open(datapath[curFinID], std :: ios :: binary);
    }

	if (fin.eof()) {
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
        if (curFinID > 60) {
            fin.close();
            exit(0);
        }
    }
    fin.read((char*)&s, sizeof(uint64_t)); //srcip(4)+dstip(4) 
    fin.read((char*)&_s, 5);//srcport destport protcol

    fin.read((char*)&t, sizeof(double));

    
    t += offset;
    if(t < lastT) {
        offset += lastT - t;
        t += (lastT - t);
    }
    lastT = t;

    if (firstTime < 0) firstTime = t;

	return make_pair(s, t - firstTime); 
}

pair <uint64_t, double> input[M + 7];

bool check_near(double c, double x) {
    return fabs(x - c) <= delta;
}

struct getGT{
    map <uint64_t, int> count;
    map <uint64_t, set <double> > table;
    
    void init() {
        count.clear();
        table.clear();
    }

    void insert(uint64_t id, double key) {
        if (key < rangeL) return;
        if (key > rangeR) return;
        if (count.find(id) == count.end()) 
            count[id] = 0, table[id].clear();
        count[id]++, table[id].insert(key);
    }

    int query(const uint64_t& id) {
        if (count.find(id) == count.end()) return 0;
        return count[id];
    }
}intervalGT;


map <uint64_t, double> timeStamp;
map <uint64_t, int> intervalAnswer;

int total_ele, total_save;
pair <int, uint64_t> ele[M + 7], save[M + 7];
set <uint64_t> topKid, allid, topKid_ours;

void GroundTruth() {
    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        double curTime = e.second;
        double lastTime = -1;
        
        if (timeStamp.find(id) != timeStamp.end()) {
            lastTime = timeStamp[id];
            intervalGT.insert(id, curTime - lastTime);
        }
        
        timeStamp[id] = curTime;
    }

    intervalAnswer.clear();
    total_ele = 0;
    for (int i = 0; i < M; ++i){
        auto e = input[i]; 
        uint64_t id = e.first;
        
        if (intervalAnswer.find(id) != intervalAnswer.end()) continue;

        allid.insert(id);

        int result = intervalGT.query(id);
        intervalAnswer[id] = result;
        ele[++total_ele] = make_pair(result, id);
      //  printf("c %.6lf per %.6lf\n", result.second, intervalGT.queryPercentage(id, result.second));
    }

    sort(ele + 1, ele + total_ele + 1);
    for (int k = total_ele, i = 1; i <= _K; ++i, --k) {
        //cerr << '#' << '#' << "rank" << i <<' ' << ele[k].first << endl;
        topKid.insert(ele[k].second);
    }

}



//int precision[2], recall[2];
int cnt;

double areK[2], aaeK[2];
double are[2], aae[2];

void OURS(int memory, int nK) {
    cerr << memory << ' ' << nK << endl;
    Alg* Detector = new Alg(3, (memory * 1024 - nK * 28) / (12 * 3), nK);  //L * 12 / 1024 = memory KB

    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        double curTime = e.second;
        Detector -> insert(id, curTime);
    }

    Detector -> getTopK(topKid_ours);

    are[0] = are[1] = 0;
    aae[0] = aae[1] = 0;
    for (auto it = allid.begin(); it != allid.end(); ++it) {
        uint64_t id = *it;
        int x = intervalAnswer[id];
        int y = Detector -> query(id);
        aae[0] += abs(x - y), aae[1]++;
       // if (y != 0) printf("###%llu %d %d\n", id, x, y);
        if (x != 0) are[0] += abs((double)x - y) / x, are[1]++; 
    }
        
    areK[0] = areK[1] = 0;
    aaeK[0] = aaeK[1] = 0;
    for (auto it = topKid.begin(); it != topKid.end(); ++it) {
        int x = intervalAnswer[*it];
        int y = Detector -> query(*it);

        aaeK[0] += abs(x - y), aaeK[1]++;
        if (x != 0) areK[0] += abs((double)x - y) / x, areK[1]++; 
    }

    int both = 0, sze = 0;
    for (auto it = topKid_ours.begin(); it != topKid_ours.end(); ++it, ++sze) 
        if (topKid.find(*it) != topKid.end()) ++both;
    //排序取出top1K
   // printf("%.6lf, %.6lf, %d,,", (double)both / _K, (double)both / sze, sze);

    delete Detector;
    printf("%.lf, %.6lf, %.lf, %.6lf, ,%.lf, %.6lf, %.lf,%.6lf\n", aae[1], aae[0]/aae[1], are[1], are[0]/are[1]
                                                     , aaeK[1], aaeK[0]/aaeK[1], areK[1], areK[0]/areK[1]);
}


void HeavyGuardian(int memory, int nK) {
    cerr <<  '@' << memory << ' ' << nK << endl;
    bloomfliter* sketch = new bloomfliter(3, memory * 1024 / (8 * 5)); //3 * M * 8 / 1024 = 3/5 * memory KB
    heavyguardian* hg = new heavyguardian(memory * 1024 * 2 / (16 * 5 * 6), 6, 0, 1.08); //W * H * 16 /1024 = 2/5 * memory KB

    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        double curTime = e.second;
        double lstTime = sketch -> query(id);
        sketch -> insert(id, curTime);
        if (curTime - lstTime > rangeR) continue;
        if (curTime - lstTime < rangeL) continue;

        hg -> insert(id);
    }

    pair<uint64_t, uint64_t>* sim_topk = hg->topk(nK);
    topKid_ours.clear();
    for (int i = 0; i < nK; ++i)
        topKid_ours.insert(sim_topk[i].second);

     are[0] = are[1] = 0;
    aae[0] = aae[1] = 0;
    for (auto it = allid.begin(); it != allid.end(); ++it) {
        uint64_t id = *it;
        int x = intervalAnswer[id];
        int y = hg -> query(id);
        aae[0] += abs(x - y), aae[1]++;
       // if (y != 0) printf("###%llu %d %d\n", id, x, y);
        if (x != 0) are[0] += abs((double)x - y) / x, are[1]++; 
    }
        
    areK[0] = areK[1] = 0;
    aaeK[0] = aaeK[1] = 0;
    for (auto it = topKid.begin(); it != topKid.end(); ++it) {
        int x = intervalAnswer[*it];
        int y = hg -> query(*it);

        aaeK[0] += abs(x - y), aaeK[1]++;
        if (x != 0) areK[0] += abs((double)x - y) / x, areK[1]++; 
    }

    int both = 0, sze = 0;
    for (auto it = topKid_ours.begin(); it != topKid_ours.end(); ++it, ++sze) 
        if (topKid.find(*it) != topKid.end()) ++both;
    
  //  cerr << (double)both / _K << endl;
 //   printf("%.6lf, %.6lf, %d,,", (double)both / _K, (double)both / sze, sze);

    delete hg;
    delete sketch;
    
    printf("%.lf, %.6lf, %.lf, %.6lf, ,%.lf, %.6lf, %.lf,%.6lf\n", aae[1], aae[0]/aae[1], are[1], are[0]/are[1]
                                                     , aaeK[1], aaeK[0]/aaeK[1], areK[1], areK[0]/areK[1]);
}

//17213472216328176400 
//1006749987763913488 
int main() {
    freopen("out.csv", "w", stdout);
    srand(0);
    for (int i = 0; i < M + 1; ++i) input[i] = Read();

    GroundTruth(); //cerr << "calc GT" << endl;
  //  for (int i = 100; i <= 2000; i += 100) 
  //       HeavyGuardian(i, _K * 2);
       // HeavyGuardian(i, _K * 2), HeavyGuardian(i, _K * 3), HeavyGuardian(i, _K * 5), putchar('\n');
    for (int i = 60; i <= 1000; i += 20) 
        OURS(i, _K);
        //OURS(i, _K * 1.5), OURS(i, _K * 2), putchar('\n'); 
        //OURS(i, _K * 2), putchar('\n');
    
    //puts("calc OURS");
    //printf("%d\n", cnt);
    //printf("#%d %d\n", total_ele, total_save);
    //printf("precision %d %d %.6lf\n", precision[1], precision[0], (double)precision[0] / precision[1]);
    //printf("recall %d %d %.6lf\n", recall[1], recall[0], (double)recall[0] / recall[1]);
   // printf("centre_precision %.lf %.6lf\n", (double)centre_precision[1], (double)(centre_precision[0] / centre_precision[1]));

	return 0;
}