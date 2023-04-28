#include <bits/stdc++.h>
using namespace std;

void updateAtoms()
{
    ifstream atomsFile;
    ofstream atomsOutputFile;
    atomsFile.open("atoms.txt");
    atomsOutputFile.open("atomsWithUpdatedCharge.txt");
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

        int atomNo, sec, type;
        float charge, x, y, z;
        string hash, s1, s2;

        atomNo = stoi(entries[0]);
        sec = stoi(entries[1]);
        type = stoi(entries[2]);
        charge = stof(entries[3]);
        x = stof(entries[4]);
        y = stof(entries[5]);
        z = stof(entries[6]);
        hash = entries[7];
        s1 = entries[8];
        s2 = entries[9];

        if(type==9)
            charge = 0;

        atomsOutputFile <<atomNo<< " " <<sec<< " " <<type<< " " <<charge<< " " <<x<< " " <<y<< " " <<z<< " " <<hash<< " " <<s1<< " " <<s2<< endl;
    }
    atomsFile.close();
    atomsOutputFile.close();
}

int main(){
    cout << "Updating....\n";
    updateAtoms();
}