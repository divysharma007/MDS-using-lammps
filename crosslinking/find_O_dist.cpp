#include<iostream>

#include<fstream>

#include<string>

#include<map>

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
    float d = 0,h=-29.8174203697059+60.28957963029227;
    if(2*abs(xd) >h){
        xd = xd - (sgn(xd)*h);
    }
    if(2*abs(yd) >h){
        yd = yd - (sgn(yd)*h);
    }
    if(2*abs(zd) > h){
        zd = zd - (sgn(zd)*h);
    }
    d = (xd)*(xd)+(yd)*(yd)+(zd)*(zd);
    return d;
}

int main() {
    struct atom C;
  
    
ifstream atomsFile;
    atomsFile.open ("atoms.txt");
    vector<atom>atom9;
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


        if(v[2]=="9"){
            atom9.push_back(temp);
        }
        if(v[0]=="2862"){
            C=temp;
        }
    }
    int cnt=0;
    atomsFile.close();
           ofstream bout;
               bout.open("C_O_distance.txt");
float distance=9999999999;
    for(int i=0; i<atom9.size(); i++){
 

    bout<<atom9[i].atomNo<<" "<<sqrt(dist(C,atom9[i]))<<endl;
    
if(sqrt(dist(C,atom9[i]))<=5){cnt++;
cout<<cnt<<" O atom within cutoff "<<atom9[i].atomNo<<" "<<sqrt(dist(C,atom9[i]))<<endl;}
distance=min(distance,sqrt(dist(C,atom9[i])));
    }
    cout<<C.atomNo<<" nearest distance"<<distance<<endl;
    bout.close();
    
}