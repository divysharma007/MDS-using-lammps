#include<bits/stdc++.h>
using namespace::std ;

const long Na = 20 ;
const long Nt = 100000 ;
const double dt = 0.1 ;

int main (int argc, char *argv[]){
	
	ifstream trajStream ;
	ofstream avgMsdStream ;
	
	long atomID ;
	int atomType,cnt=0 ;
	double MSDsum = 0.0,xcor,ycor,zcor ;
	static double x[Na][Nt], y[Na][Nt], z[Na][Nt], avgMSD[Nt] ;

	
	avgMsdStream.open("avgMSD_link.dat", ios::out) ;
	
	trajStream.open("MSD_link_win.lammpstrj", ios::in) ;
	
	for(long i=0 ; i<Nt ; i++) {
	
		for(int j=0 ; j<9 ; j++) trajStream.ignore(10000,'\n') ;
		map<long, vector<double>> mp;
		for(long j=0 ; j<Na ; j++) {
			
			trajStream>>atomType>>atomID>>xcor>>ycor>>zcor;
			mp[atomID]={xcor,ycor,zcor};
			trajStream.ignore(10000,'\n') ;
			
		}
		cnt = 0;
		for(auto val:mp){
			x[cnt][i] = val.second[0];
			y[cnt][i] = val.second[1];
			z[cnt][i] = val.second[2];
			cnt++;
			//cout<<"mp "<<val.first<<" "<<val.second[0]<<" "<<val.second[1]<<" "<<val.second[2]<<endl;
		}

	}

	trajStream.close() ;
	
	avgMSD[0] = 0.0 ;
	
	for(long i = 1 ; i < Nt ; i++) {
	
		for (long k = 0 ; k < (Nt-i) ; k++) {
		
			for(long j = 0 ; j < Na ; j++) {
			
				MSDsum += (x[j][i+k] - x[j][k])*(x[j][i+k] - x[j][k]) ; 
				MSDsum += (y[j][i+k] - y[j][k])*(y[j][i+k] - y[j][k]) ;
				MSDsum += (z[j][i+k] - z[j][k])*(z[j][i+k] - z[j][k]) ;
		
			}	
	
		}

		avgMSD[i] = MSDsum/(((double)Na)*((double)(Nt-i))) ;
	MSDsum = 0.0 ;
		if(i< 10000 && i%100 == 0) cout<<"\nCalculating MSD for time step "<<i*dt<<" ps ..."<<" MSD is "<<avgMSD[i]<<endl ;
	else if(i >= 10000 && i%1000 == 0) cout<<"\nCalculating MSD for time step "<<i*dt<<" ps ..."<<" MSD is "<<avgMSD[i]<<endl ;

	
	}
	
	for(long i = 0 ; i < Nt ; i++) 
		avgMsdStream<<dt*i<<"\t"<<avgMSD[i]<<"\n" ;
	
	avgMsdStream.close() ;
	
	return 0 ;
	
}

