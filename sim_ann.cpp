/*
    sim_ann.cpp
    Supporting material software to the article 
    "Probing single cell fermentation flux and intercellular exchange networks via
    pH-microenvironment sensing and inverse modeling"
    
    Copyright (C) 2022 V. Onesto et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
                                                                           */



#include "iostream"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "fstream"
#include "vector"
#include "string"

using namespace std;

const int maxi=300;
const int T=36;


int N[T];
int P[T];
int who[maxi][T][2];
double M[maxi][maxi][T];
double flux[maxi][T];
double fluxmax[maxi][T];
double fluxmin[maxi][T];
double fluxerr[maxi][T];
double fluxp[maxi][T];
double PH[maxi][T];
double ERR[maxi][T];
double axis[maxi][maxi][T];
double chi2[T];
double g=1e9;
double g2=1e10;       //  Lagrange Multipliers see SI 3.2 
double g3=1e13;
double beta=100000;    
double factor=800000;       // conversion factor
int track[T][2];
double tot_flux[T];



double casual(){                  // a random number uniform in (0,1)
    double   x = double  (random())/(RAND_MAX+1.0);
    return x;
}

double CHI(int t){
                    double chi2=0; 
                    for(int i=0;i<=P[t]-1;i++){
                        double sum=0;
                        for(int j=0;j<=N[t]-1;j++) sum += fluxp[j][t]*M[i][j][t];
                        
                        double res;
                        if(sum>0) res =  PH[i][t] + log10(sum);
                        else res = 100000;

                        chi2 += res*res/(2.*ERR[i][t]*ERR[i][t]);
                        }
                        double sum2=0;
                        for(int j=0;j<=N[t]-1;j++) sum2+=fluxp[j][t]*fluxp[j][t];
                        chi2 += g*sum2;
                        
                        double sum3=0;
                        for(int j=0;j<=N[t]-1;j++) sum3 += fluxp[j][t];
                        chi2 +=  g3*( sum3/double(N[t])-tot_flux[t] )*( sum3/double(N[t])-tot_flux[t]);
                        
                        return chi2; 
}









void init(){                                  //INPUTS

                      fstream file;
                      
                      file.open("numbers.dat",ios::in);
                      for(int i=0;i<=T-1;i++) file >> N[i] >> P[i];
                      file.close();
                      
                      
                      file.open("all_mat.dat",ios::in);
                      for(int t=0;t<=T-1;t++) for(int i=0;i<=P[t]-1;i++) for(int j=0;j<=N[t]-1;j++) file >> M[i][j][t];
                      file.close();
                      
                      file.open("all_flux.dat",ios::in);
                      for(int t=0;t<=T-1;t++) for(int j=0;j<=N[t]-1;j++){
                                file >> flux[j][t];
                                flux[j][t] /=factor;
                                fluxp[j][t]=flux[j][t];
                                }
                      file.close();                      
                      
                      file.open("all_PH.dat",ios::in);
                      for(int t=0;t<=T-1;t++) for(int j=0;j<=P[t]-1;j++) file >> PH[j][t] >> ERR[j][t];
                      file.close(); 
                      
                      
                      for(int t=0;t<=T-1;t++)  chi2[t] = CHI(t);
                     
                      file.open("all_pca.dat",ios::in);
                      for(int t=0;t<=T-1;t++){
                        for(int j=0;j<=N[t]-1;j++){
                                     for(int i=0;i<=N[t]-1;i++) file >> axis[j][i][t];
                                     }
                            }
                      file.close();
                     
                     
                     
                     
                      for(int t=0;t<=T-1;t++)for(int j=0;j<=maxi-1;j++) who[j][t][0]=who[j][t][1]=-1000;
                      
                      file.open("tracking.dat",ios::in);
                      
                      
                      
                      
                      int a,b,c,d;
                      while(file >> a >> b >> c >> d){
                              who[a][c][0]=b;
                              who[b][d][1]=a;
                      
                      }
                      file.close();
                   
                     file.open("tot_flux.dat",ios::in);
                     for(int t=0;t<=T-1;t++){
                                   file >> tot_flux[t];
                                   tot_flux[t]/=factor;
                                   }
                     file.close();
                   
                   
                      int l1=72;
                      int l2=74;                       
                      track[0][0]=l1;
                      track[0][1]=l2;
                      int tt=0;
                      do{
                         track[tt+1][0] = who[l1][tt][0];
                         track[tt+1][1] = who[l2][tt][0];
                         l1 = who[l1][tt][0];
                         l2 = who[l2][tt][0];
                         tt++;
                        }while(tt<T-1);   
                      for(int t=0;t<=T-1;t++) chi2[t]=CHI(t);  
}   
   


void step(double sizep){

       int t = int(casual()*T);
       
       
     for(int n=0;n<=N[t]-1;n++){
       double x =casual();
       
       

       
       
       for(int i=0;i<=N[t]-1;i++) fluxp[i][t] = flux[i][t] + sizep*(0.5-x)*axis[n][i][t];
       
       
       double diff=0;
       
       double chip = CHI(t);
       diff = chi2[t]-chip; 
       double constr=0;
       
       
       for(int i=0;i<=N[t]-1;i++){
       
       if(who[i][t][0]>-1)  constr += (flux[i][t]-flux[who[i][t][0]][t+1] )*(flux[i][t]-flux[who[i][t][0]][t+1]) -  (fluxp[i][t]-flux[who[i][t][0]][t+1])*(fluxp[i][t]-flux[who[i][t][0]][t+1]);
        if(who[i][t][1]>-1)  constr += (flux[i][t]-flux[who[i][t][1]][t-1])*(flux[i][t]-flux[who[i][t][1]][t-1]) -  (fluxp[i][t]-fluxp[who[i][t][1]][t-1])*(fluxp[i][t]-fluxp[who[i][t][1]][t-1]);   
        }
        
        
        
       
       diff += g2*constr;        
       
       
       
       double boltz = exp(beta*diff);
       x = casual();  
       if(casual()<boltz){
                   for(int i=0;i<=N[t]-1;i++) flux[i][t]=fluxp[i][t];
                   chi2[t]=chip;
                   }
       else   for(int i=0;i<=N[t]-1;i++) fluxp[i][t]=flux[i][t]; 

       }
      }
 
void Nstep(double size){

             for(int i=0;i<=3*T-1;i++) step(size);
} 
 
   
int main (){
     srand(time(0));
     init();
     int TT=100;                                  
     for(int i=0;i<=TT-1;i++){
                               
                             
                               
                Nstep(0.0001); 
                
            
                } 
                  
                
       
       
       double dx=0.1;                 
       for(int t=0;t<=T-1;t++){         
          for(int j=0;j<=N[t]-1;j++){

          fluxp[j][t]=flux[j][t];
          double diff=0;
          int okerr=0;
          double deqqa;
          do{
          fluxp[j][t]+=dx/factor;
          diff = chi2[t]-CHI(t);
          double constr=0;
       
       
       for(int i=0;i<=N[t]-1;i++){
       
       if(who[i][t][0]>-1)  constr += (flux[i][t]-flux[who[i][t][0]][t+1] )*(flux[i][t]-flux[who[i][t][0]][t+1]) -  (fluxp[i][t]-fluxp[who[i][t][0]][t+1])*(fluxp[i][t]-fluxp[who[i][t][0]][t+1]);
        if(who[i][t][1]>-1)  constr += (flux[i][t]-flux[who[i][t][1]][t-1])*(flux[i][t]-flux[who[i][t][1]][t-1]) -  (fluxp[i][t]-fluxp[who[i][t][1]][t-1])*(fluxp[i][t]-fluxp[who[i][t][1]][t-1]);   
        }
          diff += g2*constr;
          if(diff<=-0.5 && okerr==0){
                       deqqa=fluxp[j][t];
                       okerr=1;
                       }
          }while(diff>-2);
                   fluxmax[j][t] = fluxp[j][t];
          fluxp[j][t]=flux[j][t];
          diff=0;
          double della;
          okerr=0;
          do{
          fluxp[j][t]-=dx/factor;
          diff = chi2[t]-CHI(t);
          double constr=0;
       
       
       for(int i=0;i<=N[t]-1;i++){
       
       if(who[i][t][0]>-1)  constr += (flux[i][t]-flux[who[i][t][0]][t+1] )*(flux[i][t]-flux[who[i][t][0]][t+1]) -  (fluxp[i][t]-fluxp[who[i][t][0]][t+1])*(fluxp[i][t]-fluxp[who[i][t][0]][t+1]);
        if(who[i][t][1]>-1)  constr += (flux[i][t]-flux[who[i][t][1]][t-1])*(flux[i][t]-flux[who[i][t][1]][t-1]) -  (fluxp[i][t]-fluxp[who[i][t][1]][t-1])*(fluxp[i][t]-fluxp[who[i][t][1]][t-1]);   
        }
          diff += g2*constr;
           if(diff<=-0.5 && okerr==0){
                       della=fluxp[j][t];
                       okerr=1;
                       }
         }while(diff>-2);        
              fluxerr[j][t] = 0.5*(deqqa-della);
              fluxmin[j][t]=fluxp[j][t];
              fluxp[j][t]=flux[j][t];
             }  

            }         

         
           for(int t=0;t<=T-1;t++){         
                for(int j=0;j<=N[t]-1;j++) cout << factor*flux[j][t] << "  ";
                cout << endl;
                }
                
            for(int t=0;t<=T-1;t++){         
                for(int j=0;j<=N[t]-1;j++) cout << factor*fluxmax[j][t] << "  ";
                cout << endl;
                }     
             for(int t=0;t<=T-1;t++){         
                for(int j=0;j<=N[t]-1;j++) cout << factor*fluxmin[j][t] << "  ";
                cout << endl;
                }     
             for(int t=0;t<=T-1;t++){         
                for(int j=0;j<=N[t]-1;j++) cout << factor*fluxerr[j][t] << "  ";
                cout << endl;
                }    
                  



     return 0 ;

}
