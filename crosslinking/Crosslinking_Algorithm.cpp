#include<iostream>
#include<fstream>
#include<string>
#include<set>
#include<vector>
#include<math.h>
using namespace std;

//Dataformat Declaration
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

struct bond{
    int bondNo;
    int type;
    int a1;
    int a2;
};

struct angle{
    int angleNo;
    int type;
    int a1;
    int a2;
    int a3;
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

// Storing Input Equillibrated Packed Mixture Data
    ifstream atomsFile;
    atomsFile.open ("atoms.txt");
    vector<atom> atoms,atom7,atom9,atom10;
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
        if(v[2]=="7"){
            atom7.push_back(temp);
        }
        if(v[2]=="9"){
            atom9.push_back(temp);
        }
        if(v[2]=="10"){
            atom10.push_back(temp);
        }
    }
    int la=atoms.size();
    atomsFile.close();

    ifstream bondsFile;
    bondsFile.open ("bonds.txt");
    vector<bond> bonds,bond7,bond8,bond9,bond10,bond11;
    while (getline(bondsFile, line))
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
        struct bond temp;
        temp.bondNo = stoi(v[0]);
        temp.type=stoi(v[1]);
        temp.a1=stoi(v[2]);
        temp.a2=stoi(v[3]);
        bonds.push_back(temp);
        if(v[1]=="8")bond8.push_back(temp);
        if(v[1]=="11")bond11.push_back(temp);
        if(v[1]=="7")bond7.push_back(temp);
        if(v[1]=="9")bond9.push_back(temp);
        if(v[1]=="10")bond10.push_back(temp);
    }
    int l=bonds.size();
    bondsFile.close();

    ifstream anglesFile;
    anglesFile.open ("angles.txt");
    vector<angle> angles;
    while (getline(anglesFile, line))
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
        struct angle temp;
        temp.angleNo = stoi(v[0]);
        temp.type=stoi(v[1]);
        temp.a1=stoi(v[2]);
        temp.a2=stoi(v[3]);
        temp.a3=stoi(v[4]);
        angles.push_back(temp);
    }
    int al=angles.size();
    anglesFile.close();

    //Crosslinking Algorithm
    int c1=0,c2=0;//Number of Single and Double Linked c1

    // Searching for c1 Atoms
    for(int i=0;i<atom10.size();i++){
        // Searching for O1 Atoms
        for(int j=0;j<atom9.size();j++){

                //Single Crosslinking of c1-O1
                if(((atom10[i].type==10) && (atom9[j].type==9) && dist(atom10[i],atom9[j])<=5*5)){

                        //Forming c1-O1 Bond
                        struct bond nb;
                        nb.bondNo = l+1;
                        nb.type = 14;
                        nb.a1=atom10[i].atomNo;
                        nb.a2=atom9[j].atomNo;
                        bonds.push_back(nb);
                        l++;

                        //Modifying Atom type of c1 and O1 Atoms
                        atom9[j].type = 14;
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                if((*ac).atomNo == atom9[j].atomNo){
                                        (*ac).type = 14;
                                }
                        }
                        atom10[i].type = 15;
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                if((*ac).atomNo == atom10[i].atomNo){
                                        (*ac).type = 15;
                                }
                        }
                        c1++;//Number of Single Linked c1

                        //c2-c1-O1 Angle
                        int c2c1=-1;
                        for (int b=0;b<bond9.size();b++){
                                if(atom10[i].atomNo==bond9[b].a1){
                                    c2c1 = bond9[b].a2;
                                }
                                else if(atom10[i].atomNo==bond9[b].a2){
                                    c2c1 = bond9[b].a1;
                                }
                        }
                        struct angle na;
                        na.angleNo = al+1;
                        na.type = 26;
                        na.a1=atom9[j].atomNo;
                        na.a2=atom10[i].atomNo;
                        na.a3=c2c1;
                        angles.push_back(na);
                        al++;

                        //C2*-O1-c1 Angle
                        int C2O=-1;
                        for (int b=0;b<bond7.size();b++){
                                if(atom9[j].atomNo==bond7[b].a1){
                                    C2O = bond7[b].a2;
                                }
                                else if(atom9[j].atomNo==bond7[b].a2){
                                    C2O = bond7[b].a1;
                                }
                        }
                        struct angle na1;
                        na1.angleNo = al+1;
                        na1.type = 29;
                        na1.a1=atom10[i].atomNo;
                        na1.a2=atom9[j].atomNo;
                        na1.a3=C2O;
                        angles.push_back(na1);
                        al++;

                        //h-c1-O1 Angle
                        int h1=-1;
                        for (int b=0;b<bond10.size();b++){
                                if(atom10[i].atomNo==bond10[b].a1){
                                    h1 = bond10[b].a2;
                                }
                                else if(atom10[i].atomNo==bond10[b].a2){
                                    h1 = bond10[b].a1;
                                }
                        }
                        struct angle na2;
                        na2.angleNo = al+1;
                        na2.type = 27;
                        na2.a1=atom9[j].atomNo;
                        na2.a2=atom10[i].atomNo;
                        na2.a3=h1;
                        angles.push_back(na2);
                        al++;

                        int H=-1;
                        int O=-1;

                        //Finding Neighbouring HO1 Atom
                        for (int b=0;b<bond8.size();b++){
                            if(bond8[b].type == 8){
                                if(atom9[j].atomNo==bond8[b].a1){
                                    H=bond8[b].a2;

                                    //Modifying Atom Type of HO1
                                    for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 17;
                                            }
                                    }
                                    for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 17;
                                            }
                                    }

                                    //Deleting HO1-O1 Bond
                                    bond8[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bond8[b].a1==(*bb).a1 && bond8[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                }
                                else if(atom9[j].atomNo==bond8[b].a2){
                                    H=bond8[b].a1;

                                    //Modifying Atom Type of HO1
                                    for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 17;
                                            }
                                    }
                                    for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 17;
                                            }
                                    }

                                    //Deleting HO1-O1 Bond
                                    bond8[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bond8[b].a1==(*bb).a1 && bond8[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                }
                            }
                        }

                        //Finding Neighbouring o Atom
                        for(int p=0;p<bond11.size();p++){
                                if(bond11[p].type == 11){
                                    if(atom10[i].atomNo==bond11[p].a1){
                                        O = bond11[p].a2;

                                        //Modifying Atom Type of o
                                        for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 18;
                                            }
                                        }
                                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 18;
                                            }
                                        }

                                        //Modifying c1-o Bond Type
                                        bond11[p].type = 16;
                                        for(auto bc = bonds.begin();bc!=bonds.end();bc++){
                                            if(bond11[p].a1==(*bc).a1 && bond11[p].a2==(*bc).a2){
                                                (*bc).type = 16;
                                            }
                                        }
                                    }
                                    else if(atom10[i].atomNo==bond11[p].a2){
                                        O = bond11[p].a1;

                                        //Modifying Atom Type of o
                                        for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 18;
                                            }
                                        }
                                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 18;
                                            }
                                        }

                                        //Modifying c1-o Bond Type
                                        bond11[p].type = 16;
                                        for(auto bc = bonds.begin();bc!=bonds.end();bc++){
                                            if(bond11[p].a1==(*bc).a1 && bond11[p].a2==(*bc).a2){
                                                (*bc).type = 16;
                                            }
                                        }
                                    }
                                }
                            }

                            //Forming HO1-o Bond
                            struct bond nb1;
                            nb1.bondNo = l+1;
                            nb1.type = 15;
                            nb1.a1=H;
                            nb1.a2=O;
                            bonds.push_back(nb1);
                            l++;

                            //c1-o-HO1
                            struct angle na3;
                            na3.angleNo = al+1;
                            na3.type = 31;
                            na3.a1=atom10[i].atomNo;
                            na3.a2=O;
                            na3.a3=H;
                            angles.push_back(na3);
                            al++;

                            //O1-c1-o
                            struct angle na5;
                            na5.angleNo = al+1;
                            na5.type = 28;
                            na5.a1=atom9[j].atomNo;
                            na5.a2=atom10[i].atomNo;
                            na5.a3=O;
                            angles.push_back(na5);
                            al++;

                            //Deleting C2*-O1-HO1 Angle
                            for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atom9[j].atomNo==(*bb).a2){
                                    if((H==(*bb).a1 && C2O==(*bb).a3)||(H==(*bb).a3 && C2O==(*bb).a1)){
                                        (*bb).type = 50;
                                    }
                                }
                            }

                            //Modifying Angle Type of h-c1-o Angle
                            for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atom10[i].atomNo==(*bb).a2){
                                    if((h1==(*bb).a1 && O==(*bb).a3)||(h1==(*bb).a3 && O==(*bb).a1)){
                                        (*bb).type = 25;
                                    }
                                }
                            }

                            //Modifying Angle type of c2-c1-o Angle
                            for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atom10[i].atomNo==(*bb).a2){
                                    if((c2c1==(*bb).a1 && O==(*bb).a3)||(c2c1==(*bb).a3 && O==(*bb).a1)){
                                        (*bb).type = 24;
                                    }
                                }
                            }
                    }

                    //Double Crosslinking of c1-O1
                    else if((atom10[i].type==15) && (atom9[j].type==9) && dist(atom10[i],atom9[j])<=5*5){
                        //Forming c1-O1 Bond
                        struct bond nb;
                        nb.bondNo = l+1;
                        nb.type = 14;
                        nb.a1=atom10[i].atomNo;
                        nb.a2=atom9[j].atomNo;
                        bonds.push_back(nb);
                        l++;

                        //Modifying Atom type of c1 and O1 Atoms
                        atom9[j].type = 14;
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                    if((*ac).atomNo == atom9[j].atomNo){
                                        (*ac).type = 14;
                                    }
                        }
                        atom10[i].type = 16;
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                    if((*ac).atomNo == atom10[i].atomNo){
                                        (*ac).type = 16;
                                    }
                        }
                        c2++;//Number of Double Linked c1

                        //c2-c1-O1 Angle
                        int c2c1=-1;
                        for (int b=0;b<bond9.size();b++){
                                if(atom10[i].atomNo==bond9[b].a1){
                                    c2c1 = bond9[b].a2;
                                }
                                else if(atom10[i].atomNo==bond9[b].a2){
                                    c2c1 = bond9[b].a1;
                                }
                        }
                        struct angle na;
                        na.angleNo = al+1;
                        na.type = 26;
                        na.a1=atom9[j].atomNo;
                        na.a2=atom10[i].atomNo;
                        na.a3=c2c1;
                        angles.push_back(na);
                        al++;

                        //C2*-O1-c1 Angle
                        int C2O=-1;
                        for (int b=0;b<bond7.size();b++){
                                if(atom9[j].atomNo==bond7[b].a1){
                                    C2O = bond7[b].a2;
                                }
                                else if(atom9[j].atomNo==bond7[b].a2){
                                    C2O = bond7[b].a1;
                                }
                        }
                        struct angle na1;
                        na1.angleNo = al+1;
                        na1.type = 29;
                        na1.a1=atom10[i].atomNo;
                        na1.a2=atom9[j].atomNo;
                        na1.a3=C2O;
                        angles.push_back(na1);
                        al++;

                        //h-c1-O1 Angle
                        int h1=-1;
                        for (int b=0;b<bond10.size();b++){
                                if(atom10[i].atomNo==bond10[b].a1){
                                    h1 = bond10[b].a2;
                                }
                                else if(atom10[i].atomNo==bond10[b].a2){
                                    h1 = bond10[b].a1;
                                }
                        }
                        struct angle na2;
                        na2.angleNo = al+1;
                        na2.type = 27;
                        na2.a1=atom9[j].atomNo;
                        na2.a2=atom10[i].atomNo;
                        na2.a3=h1;
                        angles.push_back(na2);
                        al++;

                        //O1-c1-O1 Angle
                        int O1c1 = -1;
                        for(auto ac = bonds.begin();ac!=bonds.end();ac++){
                                if((*ac).type == 14 && (*ac).a2 == atom10[i].atomNo && atom9[j].atomNo != (*ac).a1){
                                        O1c1 = (*ac).a1;
                                }
                                else if((*ac).type == 14 && (*ac).a1 == atom10[i].atomNo && atom9[j].atomNo != (*ac).a2){
                                        O1c1 = (*ac).a2;
                                }
                        }
                        struct angle na6;
                        na6.angleNo = al+1;
                        na6.type = 30;
                        na6.a1=atom9[j].atomNo;
                        na6.a2=atom10[i].atomNo;
                        na6.a3=O1c1;
                        angles.push_back(na6);
                        al++;

                        int H=-1;
                        int O=-1;

                        //Finding Neighbouring HO1 Atom
                        for (int b=0;b<bond8.size();b++){
                            if(bond8[b].type == 8){
                                if(atom9[j].atomNo==bond8[b].a1){
                                    H=bond8[b].a2;

                                    //Deleting HO1 Atom
                                    for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                                if((*ac).atomNo == H){
                                                    (*ac).type = 50;
                                                }
                                    }
                                    for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                                if((*ac).atomNo == H){
                                                    (*ac).type = 50;
                                                }
                                    }

                                    //Deleting HO1-O1 Bond
                                    bond8[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bond8[b].a1==(*bb).a1 && bond8[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                }
                                else if(atom9[j].atomNo==bond8[b].a2){
                                    H=bond8[b].a1;

                                    //Deleting HO1 Atom
                                    for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                                if((*ac).atomNo == H){
                                                    (*ac).type = 50;
                                                }
                                    }
                                    for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                                if((*ac).atomNo == H){
                                                    (*ac).type = 50;
                                                }
                                    }

                                    //Deleting HO1-O1 Bond
                                    bond8[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bond8[b].a1==(*bb).a1 && bond8[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                }
                            }
                        }

                        //Finding Neighbouring o Atom
                        for(int p=0;p<bond11.size();p++){
                            if(bond11[p].type == 16 && atom10[i].atomNo==bond11[p].a1){
                                O = bond11[p].a2;

                                //Deleting o Atom
                                for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 50;
                                            }
                                        }
                                for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                        if((*ac).atomNo == O){
                                            (*ac).type = 50;
                                        }
                                }

                                //Deleting c1-o Bond
                                bond11[p].type = 50;
                                for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                    if(bond11[p].a1==(*bb).a1 && bond11[p].a2==(*bb).a2){
                                        (*bb).type = 50;
                                    }
                                }
                            }
                            else if(bond11[p].type == 16 && atom10[i].atomNo==bond11[p].a2){
                                O = bond11[p].a1;

                                //Deleting o Atom
                                for(auto ac = atom7.begin();ac!=atom7.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 50;
                                            }
                                        }
                                for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                        if((*ac).atomNo == O){
                                            (*ac).type = 50;
                                        }
                                }

                                //Deleting c1-o Bond
                                bond11[p].type = 50;
                                for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                    if(bond11[p].a1==(*bb).a1 && bond11[p].a2==(*bb).a2){
                                        (*bb).type = 50;
                                    }
                                }
                            }
                        }

                        //Finding HO1 Atom and Deleting o-HO1 Bond
                        int HO1 = -1;
                        for(auto ac = bonds.begin();ac!=bonds.end();ac++){
                                if((*ac).type == 15 && (*ac).a2 == O){
                                    HO1 = (*ac).a1;
                                    (*ac).type = 50;
                                }
                                else if((*ac).type == 15 && (*ac).a1 == O){
                                    HO1 = (*ac).a2;
                                    (*ac).type = 50;
                                }
                        }

                        //Deleting HO1 Atom
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                if((*ac).atomNo == HO1){
                                    (*ac).type = 50;
                                }
                        }

                        //Deleting HO1-o-c1 Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(O==(*bb).a2){
                                    if((HO1==(*bb).a1 && atom10[i].atomNo==(*bb).a3)||(HO1==(*bb).a3 && atom10[i].atomNo==(*bb).a1)){
                                        (*bb).type = 50;
                                    }
                                }
                        }

                        //Deleting c2-c1-o Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atom10[i].atomNo==(*bb).a2){
                                        if((O==(*bb).a1 && c2c1==(*bb).a3)||(O==(*bb).a3 && c2c1==(*bb).a1)){
                                                (*bb).type = 50;
                                        }
                                }
                        }

                        //Deleting h-c1-o Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atom10[i].atomNo==(*bb).a2){
                                        if((h1==(*bb).a1 && O==(*bb).a3)||(h1==(*bb).a3 && O==(*bb).a1)){
                                                (*bb).type = 50;
                                        }
                                }
                        }

                        //Deleting O1-c1-o Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atom10[i].atomNo==(*bb).a2){
                                        if((O1c1==(*bb).a1 && O==(*bb).a3)||(O1c1==(*bb).a3 && O==(*bb).a1)){
                                                (*bb).type = 50;
                                        }
                                }
                        }

                        //Deleting C2*-O1-HO1 Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atom9[j].atomNo==(*bb).a2){
                                        if((H==(*bb).a1 && C2O==(*bb).a3)||(H==(*bb).a3 && C2O==(*bb).a1)){
                                                (*bb).type = 50;
                                        }
                                }
                        }
                    }
        }
    }

    //Removing Deleted Data
    for(int ac = 0; ac < atoms.size();){
            if(atoms[ac].type == 50){
                atoms.erase(atoms.begin()+ac);
            }
            else{
                ac++;
            }
    }
    for(int ac = 0; ac < bonds.size();){
            if(bonds[ac].type == 50){
                bonds.erase(bonds.begin()+ac);
            }
            else{
                ac++;
            }
    }
    for(int ac = 0; ac < angles.size();){
            if(angles[ac].type == 50){
                angles.erase(angles.begin()+ac);
            }
            else{
                ac++;
            }
    }

    //Counting Types of Data
    int atomsn = atoms.size();
    set<int> atype;
    for(auto ac = atoms.begin();ac!=atoms.end();ac++){
            atype.insert((*ac).type);
    }
    int atomst = atype.size();

    int bondsn = bonds.size();
    set<int> btype;
    for(auto ac = bonds.begin();ac!=bonds.end();ac++){
            btype.insert((*ac).type);
    }
    int bondst = btype.size();

    int anglesn = angles.size();
    set<int> agtype;
    for(auto ac = angles.begin();ac!=angles.end();ac++){
            agtype.insert((*ac).type);
    }
    int anglest = agtype.size();

    //Printing Output
    cout<<"Data"<<" "<<"Type"<<" "<<"Value"<<endl;
    cout << "Atoms" << " " << atomst << " " << atomsn << endl;
    cout << "Bonds" << " " << bondst << " " << bondsn << endl;
    cout << "Angles" << " " << anglest << " " << anglesn << endl;

    //Creating Output Files
    ofstream bout;
    bout.open("bondsOutput.txt");
    for(int i=0;i<bonds.size();i++){
        bout<<bonds[i].bondNo<<" "<<bonds[i].type<<" "<<bonds[i].a1<<" "<<bonds[i].a2<<endl;
    }
    bout.close();
    bout.open("atomsOutput.txt");
    for(int i=0;i<atoms.size();i++){
        bout<<atoms[i].atomNo<<" "<<atoms[i].sec<<" "<<atoms[i].type<<" "<<atoms[i].charge<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<" "<<atoms[i].hash<<" "<<atoms[i].s1<<" "<<atoms[i].s2<<endl;
    }
    bout.close();

    bout.open("anglesOutput.txt");
    for(int i=0;i<angles.size();i++){
        bout<<angles[i].angleNo<<" "<<angles[i].type<<" "<<angles[i].a1<<" "<<angles[i].a2<<" "<<angles[i].a3<<endl;
    }
    bout.close();

    bout.open("PVA4_GA20_CROSSLINKED.data");
    bout<<"Atoms"<<endl;
    bout<<""<<endl;
    for(int i=0;i<atoms.size();i++){
        bout<<atoms[i].atomNo<<" "<<atoms[i].sec<<" "<<atoms[i].type<<" "<<atoms[i].charge<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<" "<<atoms[i].hash<<" "<<atoms[i].s1<<" "<<atoms[i].s2<<endl;
    }
    bout<<""<<endl;
    bout<<"Bonds"<<endl;
    bout<<""<<endl;
    for(int i=0;i<bonds.size();i++){
        bout<<bonds[i].bondNo<<" "<<bonds[i].type<<" "<<bonds[i].a1<<" "<<bonds[i].a2<<endl;
    }
    bout<<""<<endl;
    bout<<"Angles"<<endl;
    bout<<""<<endl;
    for(int i=0;i<angles.size();i++){
        bout<<angles[i].angleNo<<" "<<angles[i].type<<" "<<angles[i].a1<<" "<<angles[i].a2<<" "<<angles[i].a3<<endl;
    }
    bout.close();

    //Printing Number of Crosslinks
    cout<<"Single Linked : "<<c1<<"\nDouble Linked : "<<c2<<endl;


}
