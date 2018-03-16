/*                                                                                                                            
   :module   :: Models for the CMS HCAL showers development                                                                                                                        
   :synopsis ::                                                                                                                                                           
   :input    :: The energy of the incident particles, and a data file containing the vaues of the Sampling Fraction                                                               
                                                                                                                                                      
.. moduleauthor:: Federico Ambrogi <fambrogi@cern.ch>                                                 
.. moduleauthor:: Aditee Prabhakar Rane <aditee.prabhakar.rane@cern.ch>                                                   

***************** WORKS WITH ROOT > 6.00 !!!                                                                                                                 
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
//#include <Riostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>


//include <TString.h>

//#include "Riostream.h"
//#include <TFile.h>
//#include <TSystem.h>
#include <TRandom.h>
//#include <time.h>
//#include <TMath.h>
//#include <math.h> 

using namespace std;


//****************************************************************************************************************
// Module Shower Start - adapted from Aditee:
// - I put the values for EcalEntrance, EcalExit, HcalEntrance, HcalExit, EcalintLen, HcalintLen, Normalization as fixed (theta=0)
// - I don't know where the normalization factor comes from (it is fixed in Aditee example)
// - The only input is the Energy of the incident Pion [GeV]
//****************************************************************************************************************                                                             


double generateShowerStart(double Energy) {
  double Normalization = 0.3466;
  double Ecal_Entrance=0;
  double Ecal_Exit=23;
  double Hcal_Entrance=50;

  double Hcal_Exit=141;
  double Ecal_intLen=22.527;
  double Hcal_intLen=20.258;

  srand((unsigned)time(NULL));
  double EcalLambda = (Ecal_Exit-Ecal_Entrance)/Ecal_intLen;
  //cout << EcalLambda << endl;

  TF1 * mE = new TF1("mE", "[0]+[1]*x", 0, 100);
  mE->FixParameter(0,-5.534);
  mE->FixParameter(1,82.77);
  double MeanStop=mE->Eval(Energy);

  TF1 * sE = new TF1("sE", "[0]+[1]*log(x)", 0, 100);
  sE->FixParameter(0,3.502);
  sE->FixParameter(1,1.13);
  double StdDev=sE->Eval(Energy);

  if(MeanStop>Ecal_Entrance && MeanStop<Ecal_Exit) { MeanStop +=Ecal_Entrance; }
  else { MeanStop +=Ecal_Entrance+(Hcal_Entrance-Ecal_Exit); }

  if(MeanStop>Hcal_Exit) {
    MeanStop=Hcal_Exit;
    StdDev=2;             }

  double sample_gauss=gRandom->Gaus(MeanStop,StdDev);
  double sample_exp=gRandom->Exp(1); //exp(-x/1) 
    
  if(sample_exp<EcalLambda)
    sample_exp =Ecal_Entrance + sample_exp*Ecal_intLen;
  else
    sample_exp =Hcal_Entrance + (sample_exp-EcalLambda)*Hcal_intLen;
  //cout << "FF FIX !!!" << endl;
  if(sample_exp<sample_gauss) {
    return sample_exp; }
    else { return 0; } 


}


//void ciao(double Energy) {
//    double a = generateShowerStart(Energy);
//    cout << "Shower start is: " << a ;  }
    

//This function takes as argument the showerStart 'x' position, and returns the startingPair (integer number between 1 and 16)                                                     //Thesimplified geometry is described in the file "geometry_barrel.txt": Positions is the array of the central value of the locations of the layers, and Half_Thick the half      
//width of the layers
  

int Find_startingPair(double showerStart){
  double Positions[]  = {180.6 , 185.495, 190.915, 196.335, 201.755, 207.175, 212.595, 218.015, 223.435, 229.155, 235.175, 241.195, 247.215, 253.235, 259.255, 266.2 };
  double Half_Thick[] = {2.0   , 2.525  , 2.525  , 2.525  , 2.525  , 2.525  , 2.525  , 2.525  , 2.525  , 2.825  , 2.825  , 2.825  , 2.825  , 2.825  , 2.825  , 3.75  };
  if ( showerStart > (180.6-2.0) && showerStart < (266.2-3.75) ) {                     // define the DOWN limit of the HCAL                                                        
    for (int i = 0; i <=15; i++) {
      // cout << i << " FF Pos " << Positions[i] << " Thick " <<  Half_Thick[i] << endl ;
      double Abs_down = Positions[i] - Half_Thick[i] ;
      double Abs_up   = Positions[i+1] - Half_Thick[i+1] ;
      if (showerStart > Abs_down && showerStart <= Abs_up) {
	cout << " The showerStart is in the absorber in pair number " << i+1 << endl;
	return i+1;
      }
    }
  }
  else if ( showerStart > (266.2-3.75) && showerStart <= (266.2+3.75) ) {
    cout << " The showerStart is in the last pair of the HCAL n. 16 " << endl;
    return 16 ;
  }
  else { cout << " The showerStart is outside the HCAL !!!" << endl; // define the UP limit of the HCAL                                                                            
    //exit (EXIT_FAILURE);
    return 0;
  }
}


// This function returns a list of integer corresponding to the numbers of pairs including
// following the starting one, and for which we want to evaluate the Sampling Fraction
vector<int> Determine_Pairs(int Start) {
  vector<int> Pairs;
  for (int i = Start; i<= 16; i++) {
    //cout << "Pair: " << i << endl;
    Pairs.push_back(i);
  }
  return Pairs;
}



// Aux function to read txt files
vector<string> parseLine(string line) {
	vector<string> parsed;
	stringstream s(line);
	string tmp;
	while(s >> tmp) {
		parsed.push_back(tmp);
	}
	return parsed; }


// This function extracts the SF data from the files:
// the input is a strig "looked_pair" of the form 10_HAD_START_PAIR
// where START is the number of the starting pair and PAIR is the PAIR taken into consideration 
vector<float> extractSigmaMean (string looked_pair) {    
    vector<float> MeanSigma;

    //ifstream file ("data.txt");  // data file containing al the sigmas and means                                                                                           
    ifstream file ("SF_HAD.txt");  // data file containing al the sigmas and means
    string line;
    string name, sigma, mean;
    int position = 0;
    bool found = false;
    
    cout << "STRING TO LOOK FOR: " << looked_pair << " " << endl;   // string to be constructed from the starting point and layers e.g 10_HAD_4_7
    
    if(file.is_open()) {
        while(getline(file, line)) {                // getline extracts characters and stores them into 'line' until it reaches the character '\n'
            vector<string> parsedLine = parseLine(line);

            name  = parsedLine.at(0);         // the first entry is the wanted pair, so I compare the read entry with the input of the function
            //sigma = parsedLine.at(1);
            //mean  = parsedLine.at(2);
            if(name.compare(looked_pair) == 0) {    // string::compare returns 0 if the string are equal
              mean = parsedLine.at(2);
              sigma = parsedLine.at(1);
	      cout << "The line is: " << parsedLine.at(1) << "_" << parsedLine.at(2) << endl;
              float sigma_f = stod(sigma);           // reading the values of the parameters from the array
	      float mean_f  = stod(mean);
	      //cout << mean_f << " " << endl;
                MeanSigma.push_back(mean_f); MeanSigma.push_back(sigma_f); // preserve the order! mean then sigma
                found = true;
                break;
            }            
            position++;
        } 
	//  if(found) {
        //    cout << "DATA FOR: " << looked_pair << '\n';
        //    cout << "SIGMA:\t" << sigma << '\n';
        //    cout << "MEAN:\t" << mean << '\n';
        //    cout << "Found at line: " << position << '\n';
        //} else {
        //    cout << "Could not find: (" << looked_pair << ")\n";
        //    cout << position << " lines checked.\n";
	// }        
        file.close();
    } else {
        cout << "CANNOT OPEN FILE.\n";
    }
    return MeanSigma;
}




// This function returns a number, distributed as a Landau function with the choosen Mean value and Sigma
double Landau_SF(double Mean, double Sigma)
{
  //srand((unsigned)time(NULL));
  TRandom *r0 = new TRandom();
  double SF;
  SF = gRandom->Landau(Mean,Sigma);
  SF = r0->Landau(Mean,Sigma);
  return SF;
}
 


// This temporary function picks up a random energy from the samples available.
// This will disapear since the energy wil be read from the simulated particle 
// (Will implement interpolation between energies)
int Pick_Energy() {
  // List of energies available from the samples
  //float arrayNum[23] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2,3,4,5,10,20,30,40,50,100,200,400,600};
  //int RandIndex = rand() % 23; //generates a random number between 0 and 3

  int arrayNum[14] = {1,2,3,5,10,20,30,40,50,100,150};                                                                           
  int RandIndex = rand() % 11;
  int energy = arrayNum[RandIndex];
  cout << " I selected the energy " << energy << " " << endl;
  return energy;
}
 

// Running a test 
void Test_Workflow() {  

    cout << "Testing the workflow" << endl;

    // *** Random Energy (will be the MC input)
    int E = Pick_Energy();

    string En = std::to_string(E);

    // *** Generate shower start and extract the starting pair 
    double Start_Aditee = generateShowerStart(E);

    Start_Aditee = Start_Aditee + 130; // ********************   until I fix the units
    int S = Find_startingPair(Start_Aditee);

    string St = std::to_string(S);
    if (St == "0") {
      cout << "The shower did not start in the HCAL; skipping event" << endl;
      return; }

    // *** Select which particle (will be a loop on [HAD, MIP, PI0_1, PI0_2]
    string Part = "HAD";

    // *** Loop over the possible pairs after the starting one                                                                                                        
    vector<int> Pairs = Determine_Pairs(S);

    for (int i=0; i < Pairs.size(); i++) {
      //cout << "i" << i << " pair" << Pairs[i] << endl ; 
      string P = std::to_string(i+S);
      
      // String to look for enad extract Mean and Sigma:
      string En_Par_St_P = En + "_" + Part + "_" + St + "_" + P ;
    
      cout << " I am looking for the Landau Sigma and Mean of: " << En_Par_St_P << endl; 
      //vector<float> MeanSigma = extractSigmaMean ("50_MIP_NUM_27");
      vector<float> MeanSigma = extractSigmaMean (En_Par_St_P);

      float Mean  = MeanSigma[0] ;
      float Sigma = MeanSigma[1] ;
      double SF = Landau_SF(Mean,Sigma);

      cout << "   The MEAN is " << Mean << " and the SIGMA is " << Sigma << "; the random Landau SF is" << SF << endl;
    

    }
     


}

