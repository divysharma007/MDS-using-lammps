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
    float d = 0,h=-0.529760685729876 + 89.85923931427982;
    if(2*abs(xd) > h){
        xd = xd - (sgn(xd)*h);
    }

    
    if(2*abs(yd) > h){
        yd = yd - (sgn(yd)*h);
    }
    if(2*abs(zd) >h){
        zd = zd - (sgn(zd)*h);
    }
    d = (xd)*(xd)+(yd)*(yd)+(zd)*(zd);
    return d;
}

vector<atom> atoms, atoms8, atoms11, atoms14;
vector<bond> bonds,bonds7,bonds11,bonds14,bonds15,bonds16;
vector<angle> angles;
// void countTypes(){
//     map<int,int>mp;
//     for (int i = 0;)
// }
void getAtoms(){
ifstream atomsFile;
atomsFile.open("atoms.txt");
string line;
while (getline(atomsFile, line))
{
    vector<string> entries;
    string s = "";
    for (int i = 0; i < line.size(); i++)
    {
        if (line[i] == ' ')
        {
            entries.push_back(s);

            s = "";
        }
        else
        {
            s += line[i];
        }
    }

    entries.push_back(s);
    struct atom temp;
    temp.atomNo = stoi(entries[0]);
    temp.sec = stoi(entries[1]);
    temp.type = stoi(entries[2]);
    temp.charge = stof(entries[3]);
    temp.x = stof(entries[4]);
    temp.y = stof(entries[5]);
    temp.z = stof(entries[6]);
    temp.hash = entries[7];
    temp.s1 = entries[8];
    temp.s2 = entries[9];

    atoms.push_back(temp);
    switch(stoi(entries[2])){
        case 8:
        atoms8.push_back(temp);
        break;
        case 11:
        atoms11.push_back(temp);
        break;
        case 14:
        atoms14.push_back(temp);
        break;
    }

}
    atomsFile.close();
}
void getBonds(){
    ifstream bondsFile;
    bondsFile.open ("bonds.txt");
    string line;
    while (getline(bondsFile, line))
    {
        vector<string> entries;
        string s="";
        for(int i=0;i<line.size();i++){
            if(line[i]==' '){
                entries.push_back(s);

                s="";
            }
            else {
                s+=line[i];
            }
        }
        entries.push_back(s);
        struct bond temp;
        temp.bondNo = stoi(entries[0]);
        temp.type=stoi(entries[1]);
        temp.a1=stoi(entries[2]);
        temp.a2=stoi(entries[3]);
        bonds.push_back(temp);

        switch(stoi(entries[1])){
        case 7:
        bonds7.push_back(temp);
        break;
        case 11:
        bonds11.push_back(temp);
        break;
        case 14:
        bonds14.push_back(temp);
        break;
          case 15:
        bonds15.push_back(temp);
        break;
          case 16:
        bonds16.push_back(temp);
        break;
    }
       
    }
      bondsFile.close();


}
void getAngles(){
ifstream anglesFile;
anglesFile.open ("angles.txt");
string line;
while (getline(anglesFile, line))
    {
        vector<string> entries;
        string s="";
        for(int i=0;i<line.size();i++){
            if(line[i]==' '){
                entries.push_back(s);
                s="";
            }
            else {
                s+=line[i];

            }
        }
        entries.push_back(s);
        struct angle temp;
        temp.angleNo = stoi(entries[0]);
        temp.type=stoi(entries[1]);
        temp.a1=stoi(entries[2]);
        temp.a2=stoi(entries[3]);
        temp.a3=stoi(entries[4]);
        angles.push_back(temp);
    }
      anglesFile.close();
}
void generateOutputFiles(){
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

    bout.open("PVA4_GA20_CROSSLINKED_latest.data");
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
}
int main(){
      getAtoms();
      getBonds();
      getAngles();
      int totAtoms=atoms.size();
      int totBonds = bonds.size();
      int totAngles = angles.size();
      int cutoffDist = 5;
      // Crosslinking Algorithm

      // Number of Single and Double Linked c1
      int c1 = 0, c2 = 0;
      // Searching for c1 Atoms
      for (int i = 0; i < atoms14.size(); i++)
      {
        // Searching for O1 Atoms
            for(int j=0;j<atoms11.size();j++){
//Single Crosslinking of c1-O1
                if(((atoms14[i].type==14) && (atoms11[j].type==11) && dist(atoms14[i],atoms11[j])<=cutoffDist*cutoffDist)){
                    //Forming c1-O1 Bond type->19
                        struct bond nb;
                        nb.bondNo = totBonds+1;
                        nb.type = 19;
                        nb.a1=atoms14[i].atomNo;
                        nb.a2=atoms11[j].atomNo;
                        bonds.push_back(nb);
                        totBonds++;
                         //Modifying Atom type of c1 and O1 Atoms
                        atoms11[j].type = 18;
                         for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                if((*ac).atomNo == atoms11[j].atomNo){
                                        (*ac).type = 18;
                                }
                        }
                        atoms14[i].type = 19;
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                if((*ac).atomNo == atoms14[i].atomNo){
                                        (*ac).type = 19;
                                }
                        }
                        c1++;//Number of Single Linked c1

                        int c2c1=-1;
                        for (int b=0;b<bonds14.size();b++){
                                if(atoms14[i].atomNo==bonds14[b].a1){
                                    c2c1 = bonds14[b].a2;
                                }
                                else if(atoms14[i].atomNo==bonds14[b].a2){
                                    c2c1 = bonds14[b].a1;
                                }
                        }
                         struct angle na;
                        na.angleNo = totAngles+1;
                        na.type = 34;
                        na.a1=atoms11[j].atomNo;
                        na.a2=atoms14[i].atomNo;
                        na.a3=c2c1;
                        angles.push_back(na);
                        totAngles++;

                         //C2*-O1-c1 Angle
                        int C2O=-1;
                        for (int b=0;b<bonds7.size();b++){
                                if(atoms11[j].atomNo==bonds7[b].a1){
                                    C2O = bonds7[b].a2;
                                }
                                else if(atoms11[j].atomNo==bonds7[b].a2){
                                    C2O = bonds7[b].a1;
                                }
                        }
                            struct angle na1;
                        na1.angleNo = totAngles+1;
                        na1.type = 37;
                        na1.a1=atoms14[i].atomNo;
                        na1.a2=atoms11[j].atomNo;
                        na1.a3=C2O;
                        angles.push_back(na1);
                        totAngles++;

                        //h-c1-O1 Angle
                        int h1=-1;
                        for (int b=0;b<bonds15.size();b++){
                                if(atoms14[i].atomNo==bonds15[b].a1){
                                    h1 = bonds15[b].a2;
                                }
                                else if(atoms14[i].atomNo==bonds15[b].a2){
                                    h1 = bonds15[b].a1;
                                }
                        }
                        struct angle na2;
                        na2.angleNo = totAngles+1;
                        na2.type = 35;
                        na2.a1=atoms11[j].atomNo;
                        na2.a2=atoms14[i].atomNo;
                        na2.a3=h1;
                        angles.push_back(na2);
                        totAngles++;

                        int H=-1;
                        int O=-1;
                        //Finding Neighbouring HO1 Atom
                        for (int b=0;b<bonds11.size();b++){
                            if(bonds11[b].type == 11){
                                if(atoms11[j].atomNo==bonds11[b].a1){
                                    H=bonds11[b].a2;

                                    //Modifying Atom Type of HO1
                                    for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 21;
                                            }
                                    }
                                    for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 21;
                                            }
                                    }

                                    //Deleting HO1-O1 Bond
                                    bonds11[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bonds11[b].a1==(*bb).a1 && bonds11[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                } else if(atoms11[j].atomNo==bonds11[b].a2){
                                    H=bonds11[b].a1;

                                    //Modifying Atom Type of HO1
                                    for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 21;
                                            }
                                    }
                                    for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == H){
                                                    (*ac).type = 21;
                                            }
                                    }

                                    //Deleting HO1-O1 Bond
                                    bonds11[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bonds11[b].a1==(*bb).a1 && bonds11[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                }
                            }
                        }

                         //Finding Neighbouring o Atom
                        for(int p=0;p<bonds16.size();p++){
                                if(bonds16[p].type == 16){
                                    if(atoms14[i].atomNo==bonds16[p].a1){
                                        O = bonds16[p].a2;

                                        //Modifying Atom Type of o
                                        for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 22;
                                            }
                                        }
                                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 22;
                                            }
                                        }

                                        //Modifying c1-o Bond Type
                                        bonds16[p].type = 21;
                                        for(auto bc = bonds.begin();bc!=bonds.end();bc++){
                                            if(bonds16[p].a1==(*bc).a1 && bonds16[p].a2==(*bc).a2){
                                                (*bc).type = 21;
                                            }
                                        }
                                    } else if(atoms14[i].atomNo==bonds16[p].a2){
                                        O = bonds16[p].a1;

                                        //Modifying Atom Type of o
                                        for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 22;
                                            }
                                        }
                                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                            if((*ac).atomNo == O){
                                                    (*ac).type = 22;
                                            }
                                        }

                                        //Modifying c1-o Bond Type
                                        bonds16[p].type = 21;
                                        for(auto bc = bonds.begin();bc!=bonds.end();bc++){
                                            if(bonds16[p].a1==(*bc).a1 && bonds16[p].a2==(*bc).a2){
                                                (*bc).type = 21;
                                            }
                                        }
                                    }
                                }
                        }
                             //Forming HO1-o Bond
                            struct bond nb1;
                            nb1.bondNo = totBonds+1;
                            nb1.type = 20;
                            nb1.a1=H;
                            nb1.a2=O;
                            bonds.push_back(nb1);
                            totBonds++;

                              //c1-o-HO1
                            struct angle na3;
                            na3.angleNo = totAngles+1;
                            na3.type = 39;
                            na3.a1=atoms14[i].atomNo;
                            na3.a2=O;
                            na3.a3=H;
                            angles.push_back(na3);
                            totAngles++;

                            //O1-c1-o
                            struct angle na5;
                            na5.angleNo = totAngles+1;
                            na5.type = 36;
                            na5.a1=atoms11[j].atomNo;
                            na5.a2=atoms14[i].atomNo;
                            na5.a3=O;
                            angles.push_back(na5);
                            totAngles++;

                            //Deleting C2*-O1-HO1 Angle
                            for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atoms11[j].atomNo==(*bb).a2){
                                    if((H==(*bb).a1 && C2O==(*bb).a3)||(H==(*bb).a3 && C2O==(*bb).a1)){
                                        (*bb).type = 50;
                                    }
                                }
                            }
                            
                            //Modifying Angle Type of h-c1-o Angle
                            for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atoms14[i].atomNo==(*bb).a2){
                                    if((h1==(*bb).a1 && O==(*bb).a3)||(h1==(*bb).a3 && O==(*bb).a1)){
                                        (*bb).type = 33;
                                    }
                                }
                            }
                                //Modifying Angle type of c2-c1-o Angle
                            for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atoms14[i].atomNo==(*bb).a2){
                                    if((c2c1==(*bb).a1 && O==(*bb).a3)||(c2c1==(*bb).a3 && O==(*bb).a1)){
                                        (*bb).type = 32;
                                    }
                                }
                            }

            }

             //Double Crosslinking of c1-O1
                    else if((atoms14[i].type==19) && (atoms11[j].type==11) && dist(atoms14[i],atoms11[j])<=cutoffDist*cutoffDist){
                        //Forming c1-O1 Bond
                        struct bond nb;
                        nb.bondNo = totBonds+1;
                        nb.type = 19;
                        nb.a1=atoms14[i].atomNo;
                        nb.a2=atoms11[j].atomNo;
                        bonds.push_back(nb);
                        totBonds++;

                        //Modifying Atom type of c1 and O1 Atoms
                        atoms11[j].type = 18;
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                    if((*ac).atomNo == atoms11[j].atomNo){
                                        (*ac).type = 18;
                                    }
                        }
                        atoms14[i].type = 20;
                        for(auto ac = atoms.begin();ac!=atoms.end();ac++){
                                    if((*ac).atomNo == atoms14[i].atomNo){
                                        (*ac).type = 20;
                                    }
                        }
                        c2++;//Number of Double Linked c1

                          //c2-c1-O1 Angle
                        int c2c1=-1;
                        for (int b=0;b<bonds14.size();b++){
                                if(atoms14[i].atomNo==bonds14[b].a1){
                                    c2c1 = bonds14[b].a2;
                                }
                                else if(atoms14[i].atomNo==bonds14[b].a2){
                                    c2c1 = bonds14[b].a1;
                                }
                        }

                         struct angle na;
                        na.angleNo = totAngles+1;
                        na.type = 34;
                        na.a1=atoms11[j].atomNo;
                        na.a2=atoms14[i].atomNo;
                        na.a3=c2c1;
                        angles.push_back(na);
                        totAngles++;

                        //C2*-O1-c1 Angle
                        int C2O=-1;
                        for (int b=0;b<bonds7.size();b++){
                                if(atoms11[j].atomNo==bonds7[b].a1){
                                    C2O = bonds7[b].a2;
                                }
                                else if(atoms11[j].atomNo==bonds7[b].a2){
                                    C2O = bonds7[b].a1;
                                }
                        }

                        struct angle na1;
                        na1.angleNo = totAngles+1;
                        na1.type = 37;
                        na1.a1=atoms14[i].atomNo;
                        na1.a2=atoms11[j].atomNo;
                        na1.a3=C2O;
                        angles.push_back(na1);
                        totAngles++;

                         //h-c1-O1 Angle
                        int h1=-1;
                        for (int b=0;b<bonds15.size();b++){
                                if(atoms14[i].atomNo==bonds15[b].a1){
                                    h1 = bonds15[b].a2;
                                }
                                else if(atoms14[i].atomNo==bonds15[b].a2){
                                    h1 = bonds15[b].a1;
                                }
                        }

                        struct angle na2;
                        na2.angleNo = totAngles+1;
                        na2.type = 35;
                        na2.a1=atoms11[j].atomNo;
                        na2.a2=atoms14[i].atomNo;
                        na2.a3=h1;
                        angles.push_back(na2);
                        totAngles++;

                        //O1-c1-O1 Angle
                        int O1c1 = -1;
                        for(auto ac = bonds.begin();ac!=bonds.end();ac++){
                                if((*ac).type == 19 && (*ac).a2 == atoms14[i].atomNo && atoms11[j].atomNo != (*ac).a1){
                                        O1c1 = (*ac).a1;
                                }
                                else if((*ac).type == 19 && (*ac).a1 == atoms14[i].atomNo && atoms11[j].atomNo != (*ac).a2){
                                        O1c1 = (*ac).a2;
                                }
                        }
                         struct angle na6;
                        na6.angleNo = totAngles+1;
                        na6.type = 38;
                        na6.a1=atoms11[j].atomNo;
                        na6.a2=atoms14[i].atomNo;
                        na6.a3=O1c1;
                        angles.push_back(na6);
                        totAngles++;

                         int H=-1;
                        int O=-1;

                        //Finding Neighbouring HO1 Atom
                        for (int b=0;b<bonds11.size();b++){
                            if(bonds11[b].type == 11){
                                if(atoms11[j].atomNo==bonds11[b].a1){
                                    H=bonds11[b].a2;

                                    //Deleting HO1 Atom
                                    for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
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
                                    bonds11[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bonds11[b].a1==(*bb).a1 && bonds11[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                }else if(atoms11[j].atomNo==bonds11[b].a2){
                                    H=bonds11[b].a1;

                                    //Deleting HO1 Atom
                                    for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
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
                                    bonds11[b].type = 50;
                                    for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                        if(bonds11[b].a1==(*bb).a1 && bonds11[b].a2==(*bb).a2){
                                            (*bb).type = 50;
                                        }
                                    }
                                }
                            }
                        }

                        //Finding Neighbouring o Atom
                        for(int p=0;p<bonds16.size();p++){

                            if(bonds16[p].type == 21 && atoms14[i].atomNo==bonds16[p].a1){
                             
                                O = bonds16[p].a2;

                                // Deleting o Atom
                                for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
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
                                bonds16[p].type = 50;
                                for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                    if(bonds16[p].a1==(*bb).a1 && bonds16[p].a2==(*bb).a2){
                                        (*bb).type = 50;
                                    }
                                }
                            }
                             else if(bonds16[p].type == 21 && atoms14[i].atomNo==bonds16[p].a2){
                                O = bonds16[p].a1;

                                  //Deleting o Atom
                                for(auto ac = atoms8.begin();ac!=atoms8.end();ac++){
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
                                bonds16[p].type = 50;
                                for(auto bb = bonds.begin();bb!=bonds.end();bb++){
                                    if(bonds16[p].a1==(*bb).a1 && bonds16[p].a2==(*bb).a2){
                                        (*bb).type = 50;
                                    }
                                }
                            }
                        }
//148 30 2867 2862 2868
                         //Finding HO1 Atom and Deleting o-HO1 Bond
                        int HO1 = -1;
                        for(auto ac = bonds.begin();ac!=bonds.end();ac++){
                                if((*ac).type == 20 && (*ac).a2 == O){
                                    HO1 = (*ac).a1;
                                    (*ac).type = 50;
                                }
                                else if((*ac).type == 20 && (*ac).a1 == O){
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
                                    if((HO1==(*bb).a1 && atoms14[i].atomNo==(*bb).a3)||(HO1==(*bb).a3 && atoms14[i].atomNo==(*bb).a1)){
                                        (*bb).type = 50;
                                    }
                                }
                        }
                        
                        //Deleting c2-c1-o Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atoms14[i].atomNo==(*bb).a2){
                                        if((O==(*bb).a1 && c2c1==(*bb).a3)||(O==(*bb).a3 && c2c1==(*bb).a1)){
                                                (*bb).type = 50;
                                        }
                                }
                        }

                          //Deleting h-c1-o Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atoms14[i].atomNo==(*bb).a2){
                                        if((h1==(*bb).a1 && O==(*bb).a3)||(h1==(*bb).a3 && O==(*bb).a1)){
                                                (*bb).type = 50;
                                        }
                                }
                        }

                         //Deleting O1-c1-o Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atoms14[i].atomNo==(*bb).a2){
                                        if((O1c1==(*bb).a1 && O==(*bb).a3)||(O1c1==(*bb).a3 && O==(*bb).a1)){
                                                (*bb).type = 50;
                                        }
                                }
                        }

                                   //Deleting C2*-O1-HO1 Angle
                        for(auto bb = angles.begin();bb!=angles.end();bb++){
                                if(atoms11[j].atomNo==(*bb).a2){
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

    //Printing Number of Crosslinks
    cout<<"Single Linked : "<<c1<<"\nDouble Linked : "<<c2<<endl;
    generateOutputFiles();
}

