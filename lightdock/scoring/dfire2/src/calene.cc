#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <map>
using namespace std;


const int matype=167, mbin=30;
map<string, int> atomtype;
double edfire[matype][matype][mbin];
void rdlib(string fn){
	FILE *fp = fopen(fn.c_str(), "r");
	char str[501]; double dat[50];
	string ss[4];
	int na = 0;
	while(fgets(str, 500, fp) != NULL){
		if(str[0] == '#') continue;
		string line = str;
		istringstream iss(line);
		iss >> ss[0] >> ss[1] >> ss[2] >> ss[3];
		string an1, an2;
		an1 = ss[0] + ' ' + ss[1];
		an2 = ss[2] + ' ' + ss[3];
		if (atomtype.count(an1) == 0) atomtype[an1] = na ++;
		if (atomtype.count(an2) == 0) atomtype[an2] = na ++;
		int id1 = atomtype[an1], id2 = atomtype[an2];	// the atom pairs
		if(id1 >= matype || id2 >= matype) {
			fprintf(stderr, "something wrong with the energy file\n");
			exit(1);
		}

//		printf("%d %d\n", id1, id2);
		for(int m=0; m<mbin; m++){
			iss >> edfire[id1][id2][m];		// 30 bins of energy
			edfire[id2][id1][m] = edfire[id1][id2][m];  // symmetry
		}
	}
	fclose(fp);
}
double calENE(string fn){
	const int ma=50000;
	FILE *fp = fopen(fn.c_str(), "r");
	char str[201], rn[4], an[4];
	int na=0, nr=-1;
	int atomid[ma], resid[ma]; double x[ma][3];
	string rinfo0 = "";
	while(fgets(str, 200, fp) != NULL){
		if(strstr(str, "ATOM") != str) continue;
		string rinfo = string(str).substr(17, 10);
		if(rinfo != rinfo0) {rinfo0 = rinfo; nr ++;}
		sscanf(str+17, "%3s", rn); rn[3]='\0';
		sscanf(str+13, "%3s", an); an[3]='\0';
//		Here you can add any line to filter unwanted atoms, e.g. uncomment next line can calculate energy only between CB atoms
//		if (strcmp(an, "CB") != 0) continue;
//
		string an1 = string(rn) + ' ' + an;
		if(atomtype.count(an1) < 1) continue;
		atomid[na] = atomtype[an1];
		sscanf(str+30, "%lf%lf%lf", x[na], x[na]+1, x[na]+2);
		resid[na] = nr;
		na ++;
		if(na > ma) {fprintf(stderr, "increase ma\n"); exit(0);}
	}
	fclose(fp);
//
	double eall = 0.;
	int count = 0;
	for(int i=0; i<na; i++)
	for(int j=i+1; j<na; j++){
		if(resid[i] == resid[j]) continue;		// ignore interaction with same residue
		double r = 0.;
		for(int m=0; m<3; m++){
			double xd = x[i][m] - x[j][m];
			r += xd*xd;
		}
		r = sqrt(r);
		int b = int(r*2);				// distance to bin
		if(b >= 30) continue;
		//printf("%5.3f %d %d %d %5.3f %d %d\n", r, atomid[i], atomid[j], b, edfire[atomid[i]][atomid[j]][b], resid[i], resid[j]);
		++count;
		eall += edfire[atomid[i]][atomid[j]][b];
	}
	printf("%s %.1f %d %d %d\n", fn.c_str(), eall/100., na, count, nr);
	return eall;
}
int main(int argc, char *argv[]){
	if(argc < 2) {fprintf(stderr, "usage: %s dfire_pair.lib PDBs\n", argv[0]); exit(1);}
	rdlib(argv[1]);		// read the potential file into matrix
	for(int i=2; i<argc; i++){
		calENE(argv[2]);		// cal the energy
	}
}
