#include<iostream>

#include<fstream>

#include<string>

#include<set>

#include<vector>

#include<math.h>

using namespace std;
struct atom {
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
int sgn(float x) {
  if (x > 0) {
    return 1;
  } else if (x < 0) {
    return -1;
  } else {
    return 0;
  }
}
float dist(struct atom a1, struct atom a2) {
  float xd, yd, zd;
  xd = a1.x - a2.x;
  yd = a1.y - a2.y;
  zd = a1.z - a2.z;
  float d = 0, h = -29.8174203697059 + 60.28957963029227;
  if (2 * abs(xd) > h) {
    xd = xd - (sgn(xd) * h);
  }
  if (2 * abs(yd) > h) {
    yd = yd - (sgn(yd) * h);
  }
  if (2 * abs(xd) > h) {
    zd = zd - (sgn(zd) * h);
  }
  d = (xd) * (xd) + (yd) * (yd) + (zd) * (zd);
  return d;
}
int main() {
  ifstream atomsFile;
  cout << "working\n";
  atomsFile.open("atomsAfterLink.txt");
  vector < atom > atoms;
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
    struct atom temp;
    temp.atomNo = stoi(v[0]);
    temp.sec = stoi(v[1]);
    temp.type = stoi(v[2]);
    temp.charge = stof(v[3]);
    temp.x = stof(v[4]);
    temp.y = stof(v[5]);
    temp.z = stof(v[6]);
    temp.hash = v[7];
    temp.s1 = v[8];
    temp.s2 = v[9];
    atoms.push_back(temp);

  }

  atomsFile.close();
  ifstream findAtomsFile;
  findAtomsFile.open("bondsOutput.txt");
  vector < pair < int, int > > findAtoms;
  while (getline(findAtomsFile, line)) {
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
    int t = stoi(v[1]);
    int d = stoi(v[2]);
    int e = stoi(v[3]);
    if (t == 14) {
      pair < int, int > p;
      p.first = d;
      p.second = e;

      findAtoms.push_back(p);
    }

  }

  findAtomsFile.close();
  double dis = 0;
  int cnt = 0;

  for (int l = 0; l < findAtoms.size(); l++) {
    for (int i = 0; i < atoms.size(); i++) {
      if (findAtoms[l].first == atoms[i].atomNo && (atoms[i].type == 16 || atoms[i].type == 15)) {
        for (int j = 0; j < atoms.size(); j++) {
          if (findAtoms[l].second == atoms[j].atomNo) {
            cnt++;
            cout << atoms[i].atomNo << " " << atoms[j].atomNo << endl;
            dis += sqrt(dist(atoms[i], atoms[j]));
          }
        }
      }
    }
  }

  cout << "distance: " << dis << " cnt: " << cnt << " Average: " << (dis / (float) cnt) << endl;
}