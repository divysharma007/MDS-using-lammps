#include<iostream>

#include<fstream>

#include<string>

#include<map>

#include<vector>

#include<math.h>

using namespace std;

void printAtoms(){
ifstream atomsFile;
  cout << "counting Atoms\n";
  vector<int>mp(19);
  atomsFile.open("atomsOutput.txt");
  string line;
  while (getline(atomsFile, line)) {
    vector < string > v;
    string s = "";
    for (int i = 0; i < line.size(); i++) {
      if (line[i] == ' ') {
        v.push_back(s);
        s = "";
      } else {
        s += line[i];
      }
    }
    v.push_back(s);
    mp[stoi(v[2])]++;

  }



  atomsFile.close();
  for(int i=1;i<=18;i++){
    cout<<i<< " "<<mp[i]<<endl;
  }
}
void printBonds(){
  ifstream atomsFile;
  cout << "counting Bonds\n";
  vector<int>mp(17);
  atomsFile.open("bondsOutput.txt");
  string line;
  while (getline(atomsFile, line)) {
    vector < string > v;
    string s = "";
    for (int i = 0; i < line.size(); i++) {
      if (line[i] == ' ') {
        v.push_back(s);
        s = "";
      } else {
        s += line[i];
      }
    }
    v.push_back(s);
    mp[stoi(v[1])]++;

  }



  atomsFile.close();
  for(int i=1;i<=16;i++){
    cout<<i<< " "<<mp[i]<<endl;
  }
}

void printAngles(){
  ifstream atomsFile;
  cout << "counting Angels\n";
  vector<int>mp(32);
  atomsFile.open("anglesOutput.txt");
  string line;
  while (getline(atomsFile, line)) {
    vector < string > v;
    string s = "";
    for (int i = 0; i < line.size(); i++) {
      if (line[i] == ' ') {
        v.push_back(s);
        s = "";
      } else {
        s += line[i];
      }
    }
    v.push_back(s);
    mp[stoi(v[1])]++;

  }



  atomsFile.close();
  for(int i=1;i<=31;i++){
    cout<<i<< " "<<mp[i]<<endl;
  }
}
int main() {
  printAngles();
}