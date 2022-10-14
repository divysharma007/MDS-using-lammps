#include<iostream>
#include<fstream>
#include<string>
#include<set>
#include<vector>
#include<math.h>
using namespace std;
struct atom{
    int atomNo;
    int sec;
    int type;
    float charge;
    float x;
    float y;
    float z;
    string hash;
    string s1;
    string s2;
};
//Distance between Atoms with Periodic Boundary Condition
int sgn(float x){
if(x > 0){
    return 1;
}
else if(x < 0){
    return -1;
}
else {
    return 0;
}
}
float dist(struct atom a1, struct atom a2){
    float xd, yd, zd;
    xd = a1.x - a2.x;
    yd = a1.y - a2.y;
    zd = a1.z - a2.z;
    float d = 0;
    if(abs(xd) > 15.57115){
        xd = xd - (sgn(xd)*31.1423);
    }
    if(abs(yd) > 15.57115){
        yd = yd - (sgn(yd)*31.1423);
    }
    if(abs(xd) > 15.57115){
        zd = zd - (sgn(zd)*31.1423);
    }
    d = (xd)*(xd)+(yd)*(yd)+(zd)*(zd);
    return d;
}
int main(){
 ifstream atomsFile;
 cout<<"eded\n";
    atomsFile.open ("atomsAfterLink.txt");
    vector<atom> atoms;
    string line;
    while (getline(atomsFile, line))
    {
        vector<string> v;
        string s="";
        for(int i=0;i<line.size();i++){
            if(line[i]==' '){
                v.push_back(s);
                s="";
            }
            else {
                s+=line[i];
            }
        }
        v.push_back(s);
        struct atom temp;
        temp.atomNo = stoi(v[0]);
        temp.sec = stoi(v[1]);
        temp.type=stoi(v[2]);
        temp.charge=stof(v[3]);
        temp.x=stof(v[4]);
        temp.y=stof(v[5]);
        temp.z=stof(v[6]);
        temp.hash = v[7];
        temp.s1 = v[8];
        temp.s2 =  v[9];
        atoms.push_back(temp);

    }

    atomsFile.close();
    ifstream findAtomsFile;
    findAtomsFile.open("findAtoms.txt");
    vector< pair<int,int> >findAtoms;
        while (getline(findAtomsFile, line))
    {
        vector<string> v;
        string s="";
        for(int i=0;i<line.size();i++){
            if(line[i]==' '){
                v.push_back(s);
                s="";
            }
            else {
                s+=line[i];
            }
    
        }
        v.push_back(s);

        int d=stoi(v[0]);
        int e=stoi(v[1]);
        pair<int,int> p;
        p.first = d; p.second = e;
                cout<<d<<" "<<e<<endl;
        findAtoms.push_back( p );

    }

    findAtomsFile.close();
    double dis=0;
    int cnt=0;

for(int l=0;l<findAtoms.size();l++){
 for(int i=0;i<atoms.size();i++){
if(findAtoms[l].first==atoms[i].atomNo && atoms[i].type==16){
        for(int j=0;j<atoms.size();j++){
if(findAtoms[l].second==atoms[j].atomNo){
    cnt++;atoms[i].type=50;
dis+=sqrt(dist(atoms[i],atoms[j]));
        }}
}
 }}



cout<<"distance: "<<dis<<" cnt: "<<cnt<<" Average: "<<(dis/(float)cnt)<<endl;
}
