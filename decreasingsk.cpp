//_____________________________________LIBRERIE______________________________________//
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <iomanip>
using namespace std;

#include "bnldev_mod.c"
#include "poidev_mod.c"
#include "gammln_mod.c"
#include "gammln.c"
#include "mt19937-64.c"

//__________________________________________________________________________________________//
//______________________________________*** MAIN ***________________________________________//
//__________________________________________________________________________________________//

int main() {

  init_genrand64 ( time(NULL) );
  ofstream wmeantotFile;   wmeantotFile.open("hoc_wmeantot.dat", ios::trunc);
  ofstream kmeantotFile;   kmeantotFile.open("hoc_kmeantot.dat", ios::trunc); 
  ofstream advtotFile;     advtotFile.open("hoc_advtot.dat", ios::trunc);
  ofstream advtot2File;    advtot2File.open("hoc_advtot2.dat", ios::trunc); 
  ofstream vktotFile;      vktotFile.open("hoc_vktot.dat", ios::trunc);
  ofstream vstotFile;      vstotFile.open("hoc_vstot.dat", ios::trunc);
  ofstream histoFile;      histoFile.open("hoc_histo.dat", ios::trunc);
  ofstream histowFile;     histowFile.open("hoc_histow.dat", ios::trunc);
  ofstream numclassFile;   numclassFile.open("hoc_numclass.dat", ios::trunc);

  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(8); 

  //Parametri dal fitting dei dati di Enoch
  /*  double N0 = 3E+10, N = N0;
  double s0 = 0.42, s = s0;
  double slog = 0.297;
  double sgeo = 0.25;
  double alpha = 0.19;
  double qq = 0.60;*/  

 //____Dichiarazione delle variabili globali (iterazioni ensemble)____// 
  int VarPar = 1;
  int iter = 25;
  double N0 =3.3E+7, N = N0;
  double U0 = 1.51E-05, U = U0;
  int t0 = 0;
  int tmax = 50001;
  int step = floor (tmax / 20);
  int cell = tmax / step;
  int cc = 0;
  int counter0 = 20000; int counter = counter0;

  double s0 = 0.16, s = s0;
  //double slog = 0.083;
  //double sgeo = 0.16;
  double alpha = 0.26;
  //double qq = 0.02521;  
  
  //____Dichiarazione delle variabili di popolazione (iterazioni temporali)____//
  int k = 1;
  int kmax = 2*tmax;
  int MaxMut=0;
  double wmean = 0.;
  double wmean2 = 0.;
  double kmean = 0.;
  double kmean2 = 0.;
  double adv = 0.;
  int numclass = 0;

   //____Dichiarazione delle variabili di classe (iterazioni sulle frequenze)____//
  double nr, Nres = N;
  double pr=0., qr=0., sumpr=0.;
  
  double *Vwmeantot = new double[tmax];
  double *Vwmean2tot = new double [tmax];
  double *Vkmeantot = new double[tmax+1];
  double *Vkmean2tot = new double[tmax];
  double *Vwmeanvar = new double[tmax];
  double *Vkmeanvar = new double[tmax];

  double *Vadvtot = new double[tmax];
  double *Vvel_k = new double[cell]; 
  double *Vvel_s = new double[cell]; 
  double *Vnumclasstot = new double[tmax];
 
  double *Vfitness = new double[kmax]; //vettore delle wr fitness delle classi r {w0,w1,w2,w3...}
  double *Vadvantage = new double[kmax]; //vettore delle s(k) delle classi r {s0(k),s1(k),s2(k)...}



  //_____________________________ ||||| FOR SUI PARAMETRI ||||| _____________________________//
  for (int m = 1; m<=VarPar; m++) {

    cout << "Variazione di parametro #" << m << "/" << VarPar << endl;

    //___________Inizializzazione della fitness_____________//
    for (int i=0; i<kmax; i++) {
      //Vadvantage[i] = s*i;                                   //fitness a vantaggio fisso
      Vadvantage[i] = s * pow((double)i, alpha);               //fitness power law
      //Vadvantage[i] = slog * log(i + 1);                          //fitness logaritmica
      //Vadvantage[i] = sgeo * (1 - pow(qq,(double)i)) / (1-qq);           //fitness geometrica
      
      Vfitness[i] = exp(Vadvantage[i]); 
    }
    
       //___________________________ **** FOR SUGLI ENSEMBLE ****_____________________________//
    for (int a=1; a<=iter; a++) {
      cout << "Iterazione #" << a << "/" << iter << "\r" << flush;

	 //___________Inizializzazione dei vettori_____________//

  vector<int> Vmutation (kmax,0); //semplice vettore del # di mutazioni {0,1,2,3...}	 
  for (int i=0; i<kmax; i++) {
    Vmutation[i]=i;
	   //Vprob[i]=0.;
	   //Vqrob[i]=0.;
	   //Vfreq[i]=0.;
	   //Vprob_new[i]=0.;
	   //Vqrob_new[i]=0.;
	   //Vfreq_new[i]=0.;
  }

  vector <double> Vprob (kmax,0); //vettore delle pr probabilita'  delle classi r {p0,p1,p2,p3...}
  vector <double> Vqrob (kmax,0); //vettore delle qr prob ridotte delle classi r {q0,q1,q2,q3...}
  vector <double> Vfreq (kmax,0); // vettore delle fr frequenze delle classi r {f0,f1,f2,f3...}
  vector <double> Vprob_new (kmax,0); //vettore delle fr delle classi r {f0,f1,f2,f3...}
  vector <double> Vqrob_new (kmax,0); //vettore delle qr delle classi r {q0,q1,q2,q3...}
  vector <double> Vfreq_new (kmax,0); //vettore delle fr delle classi r {f0,f1,f2,f3...}

  Vfreq[0] = 1.;
  k = 1;
  counter = counter0; 

     //___________________________ °°° FOR TEMPORALE °°° ____________________________//
  for(int t=t0; t<tmax; t++) {
   k++;

   wmean = 0.;
   wmean2 = 0.;
   kmean = 0.;
   kmean2 = 0.;
   adv = 0.;
   for (int i=0; i<k; i++) { 
    if(Vfitness[i]<5000){
     wmean += (Vfreq[i]*Vfitness[i]);
     wmean2 += (Vfreq[i]*Vfitness[i]*Vfitness[i]);
   };
   kmean += (Vfreq[i]*Vmutation[i]);
   kmean2 += (Vfreq[i]*Vmutation[i]*Vmutation[i]);
   adv += (Vfreq[i]*Vadvantage[i]);
 };

	//_______________ ^^ FOR SULLE CLASSI DI FREQUENZA ^^ _________________//
 for (int r=0; r<k; r++) {

   if(Vfitness[r]>10) {break;};

   if (Vfreq[r] != 0.) {numclass++;};

   if (r!=0) {
     pr  = ((Vfreq[r-1] * U * Vfitness[r-1]) + (Vfreq[r] * (1-U) * Vfitness[r]))/wmean;
   }
   else 
     pr = (Vfreq[r] * (1-U) * Vfitness[r])/wmean;

   sumpr += pr;
   if (pr!=0.)
     { qr= pr/sumpr; 
       MaxMut=r+1;
	                //  if(kmin==0){ kmin=r;};
     }
     else {
       qr = 0.;
     };
	  //if (pr==0.) continue;

     Vprob_new[r] = pr;
     Vqrob_new[r] = qr;
     Vprob[r] = Vprob_new[r];
     Vqrob[r] = Vqrob_new[r];

     if (pr!=0.){
	            // cout << "t= " << t << "; k=" << Vmutation[r] << "; fr=" << Vfreq[r] <<  "; pr= " << Vprob_new[r]  << "; sumpr=" << sumpr << "; qr=" << Vqrob_new[r]  << "; wr=" << Vfitness[r] << endl;
     };

     if (r>0   &&   Vprob[r-1] != 0   &&   Vprob[r] == 0){ r=k+1;};	 
   };

   qr = 0.;
	  	//cerr << " qui 3" << endl;

	  //MaxMut = k;
	  //	  kminimo=kmin;
	  //cout << "sono kmin " << kmin << " sono kminimo " << kminimo << endl;

	  //___________ + FOR SULLA BINOMIALE + ____________//
   for (int i = MaxMut; i>=0; i--) {
     qr = Vqrob_new[i];
     if((qr*Nres) < 3){	
       if((qr*Nres)< 1){
        if((qr*Nres)< 1E-4){ nr=0;
        }else{ 
          if((qr*Nres)< 1E-2)
          {
            nr= floor(double(poidev_mod(Nres*qr*1000))/1000.);
          }else{
            nr= floor(double(poidev_mod(Nres*qr*10))/10.);
          };
        };
      }else{
        nr= poidev_mod(Nres*qr);
      };
    }else{
     nr= bnldev_mod(qr, Nres);
   };
	    /* if((qr*Nres)< 1E-2)
	       nr=0.; 
	       else
	       nr = bnldev_mod(qr, Nres);*/

         double freq =nr/N;	  
       Vfreq_new[i] = freq;
       Nres -= nr;
       nr = -1.;
       Vfreq[i] = Vfreq_new[i];

       if (t == counter) {
	      	   //cerr << " qui 4" << endl;

         histoFile << i << " " << Vfreq[i] << endl;
         histowFile << Vfitness[i] << " " << Vfreq[i] << endl;
       };

	  }; //__ + Chiusura del for sulla binomiale + __ 


	  //Vfreq[r] = pr; //Scommentare questa e commentare il for sulla binomiale per l'algoritmo deterministico
	  pr = 0.;
	  qr = 0.;
	  Nres = N;


	  //if (r>0   &&   Vprob[r-1] != 0   &&   Vprob[r] == 0) break;

	  //}; // __ ^^ Chiusura del for sulle classi di frequenza ^^ __

	  if (t == counter){
      histoFile << endl;
      counter += counter0;
    };

    Vwmeantot[t] += wmean;
    Vwmean2tot[t] += wmean2;
    Vkmeantot[t] += kmean;
    Vkmean2tot[t] += kmean2;
    Vadvtot[t] += adv;
	   Vnumclasstot[t] += /*(double)*/ numclass;
    sumpr = 0.;
    numclass = 0; 

	//			if (t % 500 == 0) {
	//	  cout << t << " " << kmean << endl;
	//	}
	//cout << endl;

      }; // __ °°° Chiusura del for temporale °°° __



    };  // __ **** Chiusura della singola iterazione **** __

    cc = 0; 

    for (int n=0; n<tmax; n++) {

      Vwmeantot[n] /= iter;
      Vwmean2tot[n] /= iter;
      Vkmeantot[n] /= iter;
      Vkmean2tot[n] /= iter;
      Vadvtot[n] /= iter;
      Vnumclasstot[n] /= iter;
      Vwmeanvar[n] = sqrt(abs(Vwmean2tot[n]-(Vwmeantot[n]*Vwmeantot[n])));
      Vkmeanvar[n] = sqrt(abs(Vkmean2tot[n]-(Vkmeantot[n]*Vkmeantot[n])));

      if ((n % 500==0 && n<=5000) || (n % 1000 == 0 && n<=20000) || (n % 2000 == 0))  {
	 //cerr << n << " " << Vwmeantot[n] << " " << Vwmeanvar[n] << endl;
        wmeantotFile << n << " " << Vwmeantot[n] << " " << Vwmeanvar[n] << endl;
        kmeantotFile << n << " " << Vkmeantot[n] << " " << Vkmeanvar[n] << endl;
        advtotFile << n << " " << Vadvtot[n] << endl;
        advtot2File << Vkmeantot[n] << " " << Vadvtot[n] << endl;
        numclassFile << n << " " << Vnumclasstot[n] << endl;

        if (n>0) {
        //cerr << "step " << step << endl;
        //cerr << "n " << n << endl;
        //cerr << (n-step) <<endl;
	    //Vvel_k[cc] = (Vkmeantot[n] - Vkmeantot[n-step])/step;

        //vktotFile << n << " " << Vvel_k[cc] << endl;
          vktotFile << n << " " << ((Vkmeantot[n] - Vkmeantot[n-step])/double(step)) << endl;
	    //Vvel_s[cc] = (Vadvtot[n] - Vadvtot[n-step])/double(step);
	    //vstotFile << n << " " << Vvel_s[cc] << endl;
          vstotFile << n << " " << (Vadvtot[n] - Vadvtot[n-step])/double(step) << endl;

          if (cc == cell){
           vktotFile << endl;
           vstotFile << endl;
         };

       };
	 //cout << "qui no 2" << endl;
       cc++;
     };

     if (n == (tmax-1)) {
       wmeantotFile << endl;
       kmeantotFile << endl;
       advtotFile << endl;
       numclassFile << endl;
     }
   }

   vktotFile << endl;
   vstotFile << endl;

  //cout << U << " " << vtheo << endl;
   U *= 2;

  }; // __ ||||| Chiusura del for sui parametri ||||| __

  wmeantotFile.close();
  kmeantotFile.close();
  histoFile.close();
  histowFile.close();
  advtotFile.close();
  advtot2File.close();
  vktotFile.close();
  vstotFile.close();
  numclassFile.close();
  
  delete Vwmeantot;
  delete Vwmean2tot;
  delete Vkmeantot;
  delete Vkmean2tot;
  delete Vkmeanvar;
  
  delete Vadvtot;
  delete Vvel_k;
  delete Vvel_s;
  delete Vnumclasstot;

  /*
  delete Vmutation;
  delete Vprob;
  delete Vqrob;
  delete Vfreq;
    */
  delete Vfitness;
  delete Vadvantage;

  //delete Vprob_new;
  //delete Vqrob_new;
  //delete Vfreq_new;

  return 0;
  
}//Chiusura del main






// v teorica, con la gaussiana
  //    for (int i=0; i < cell; i++){ 
  //  vgauss = 2*s*s*log(N)/(log(U)*log(U));
  //vktotFile << i*step << " " << vgauss << endl; 
  //}
