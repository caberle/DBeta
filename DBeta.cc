// DBeta Generator.
//
// Note: Currently neglecting the recoil of the nucleus in the energy 
// conservation given the angle.
//
// An event in HEPEVT ASCII format looks like this:
//
//    NHEP 
//    ISTHEP IDHEP JDA1 JDA2 PX PY PZ PMASS T0 X0 Y0 Z0 POLX POLY POLZ 
//    ISTHEP IDHEP JDA1 JDA2 PX PY PZ PMASS T0 X0 Y0 Z0 POLX POLY POLZ 
//
// Author: Lindley Winslow
// Date: May 24, 2012
// 

#include <stdlib.h>
#include <iostream>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Random/Randomize.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"

using namespace std;
using namespace CLHEP;

//input
TString isotope_name="128te"; //We have 116Cd, 128Te, 130Te and 136Xe. For 160Gd only 0NuBB input data is currently available. 
Int_t mass_number=128;        
Int_t atomic_number=52;       //Cd:48, Te:52, Xe:54, Gd:64
//end of input 

TH2F* his2D;
Hep3Vector IsotropicDist();

Double_t getEndpoint(Int_t A, Int_t Z){

  Double_t endpoint=1.0;

  //Values from J. Kotila and F. Iachello paper: PHYSICAL REVIEW C 85, 034316 (2012)
  if(A==160 && Z==64) endpoint=1.72969; //160Gd, we didn't get the differential probability file for this isotope which is needed for 2NuBB. 
  else if(A==116 && Z==48) endpoint=2.81350; //116Cd, value from COBRA: 2.809
  else if(A==128 && Z==52) endpoint=0.86587; //128Te
  else if(A==130 && Z==52) endpoint=2.52697; //130Te
  else if(A==136 && Z==54) endpoint=2.45783; //136Xe
  else cout << "The specified isotope is not supported. Have a look at the function getEndpoint(A,Z) for a list of supported isotopes." << endl; 
  return endpoint;

}

void fillEnergy(Bool_t isSum, Bool_t is0NuBB, TString isotope, TGraph* gE){

  char* datadir = getenv("DBETA_DATA");
  
  //Open up a file - Sum First.
  ifstream infile;

  //Get the right data file;
  TString fileName=datadir;
  fileName.Append(Form("/%s",isotope.Data()));
  if(isSum){
    fileName.Append("sumssfs");
  }else{
    fileName.Append("sessfs");
  }
  if(is0NuBB)  fileName.Append("_0nuBB");
  fileName.Append(".dat");

  //Go Ahead and Open.
  infile.open(fileName.Data());
  if(!infile){
    cout << "[ERROR] Could not find energy file "<< fileName.Data() <<endl;
    return;
  }


  //These should be reset.
  Int_t npoints=0;
  Double_t xp=-1;
  Double_t sum=0;
  Double_t max=0;
  //
  Int_t index;
  Double_t x,dx;
  Double_t y;
  string line;
  while (getline (infile, line)){
    istringstream linestream(line);
    linestream >> index >> x >> y;
    if(xp>0) {
      dx=x-xp;
      sum+=y*dx;
      if(y>max) max=y;
    }
    if(x>0) {
      gE->SetPoint(npoints, x, y);
      npoints++;
    }
    xp=x;
  }
  
  //Always Normalize.
  for(Int_t i=0; i<gE->GetN(); i++){
    gE->GetPoint(i, x, y);
    gE->SetPoint(i, x, y/max);
  }
  infile.close();

} 

void fillEnergy2D(Bool_t is0NuBB, TString isotope, TGraph2D* gE){

  char* datadir = getenv("DBETA_DATA");

  //Open up a file - Sum First.
  ifstream infile;

  //Get the right data file; for example 116cdtwodimspectsfs.dat
  TString fileName=datadir;
  fileName.Append("/2D_differential_probabilities");
  fileName.Append(Form("/%s",isotope.Data()));
  fileName.Append("twodimspectsfs");
  if(is0NuBB)  fileName.Append("_0nuBB");
  fileName.Append(".dat");

  //Go Ahead and Open.
  infile.open(fileName.Data());
  if(!infile){
    cout << "[ERROR] Could not find the 2D differential probability file. " << fileName.Data() <<endl;
    return;
  }
  //cout << "Open file " << fileName.Data() << endl;

  //These should be reset.
  Int_t npoints=0;
  //Double_t xp=-1;
  //Double_t sum=0;
  Double_t max=0;
  //
  Int_t index1,index2;
  Double_t x,dx;
  Double_t y;
  Double_t prob;
  string line;
  Double_t last_x,last_y;
  Double_t stepsize_x=0.0;
  Double_t stepsize_y=0.0;
  Int_t stepsize_x_counter=0;
  Int_t stepsize_y_counter=0;
  
  while (getline (infile, line)){    //get the input data from file and fill it into TGraph2D, also get max. probability, and the stepsize in both x and y directions  
    istringstream linestream(line);
    //cout << "npoints= " << npoints << endl;
    linestream >> index1 >> index2 >> x >> y >> prob;
    //cout << "index1= " << index1 << ", index2= " << index2 << ", x= " << x << ", y= " << y << ", prob= " << prob << endl;
    if(prob>max) max=prob;

    if(x>0 && y>0) {
      gE->SetPoint(npoints, x, y, prob);
      npoints++;
    }
    else {
    cout << "The input file format is not what we expect. At least one of the variables used for the electron energies is negative." << endl;}

    //get the stepsize in both x and y direction: search file for changing x and y and use the difference when they change the second time. 
    if(last_x!=x && stepsize_x_counter<2){
      stepsize_x=x-last_x;
      cout <<"stepsize_x= " << stepsize_x << endl;
      stepsize_x_counter+=1;}
    if(last_y!=y && stepsize_y_counter<2){
      stepsize_y=y-last_y;
      cout <<"stepsize_y= " << stepsize_y << endl;
      stepsize_y_counter+=1;}
    last_x=x;
    last_y=y;
  }


  //With the stepsizes and the endpoint we can determine the binning of the 2D histogram 
  Double_t lower_limit_x=0.0; //this is how the input data files from J. Kotila are organized
  Double_t lower_limit_y=0.0;
  Int_t nbins_x = getEndpoint(mass_number,atomic_number)/stepsize_x +10;  //add 10 to be on the safe side
  Int_t nbins_y = getEndpoint(mass_number,atomic_number)/stepsize_y +10;  //add 10 to be on the safe side
  Double_t upper_limit_x=stepsize_x*nbins_x;
  Double_t upper_limit_y=stepsize_y*nbins_y;
  //cout << "upper_limit_x= " << upper_limit_x << endl;
  //cout << "upper_limit_y= " << upper_limit_y << endl;
  his2D = new TH2F("his2D","2D differential probability distribution;E1 [MeV];E2 [MeV]; normalized diff. prob.",nbins_x,lower_limit_x,upper_limit_x,nbins_y,lower_limit_y,upper_limit_y);
  //cout << "gE->GetN()= " << gE->GetN() << ", max = " << max << ", gE->GetZ()[0]= " << gE->GetZ()[0] << endl;

  //normalize the TGraph2D to max. probability.
  for(Int_t i=0; i<gE->GetN(); i++){
    //gE->GetPoint(i, x, y);
    Double_t x_temp = gE->GetX()[i];
    Double_t y_temp = gE->GetY()[i];
    Double_t prob_temp = gE->GetZ()[i];
    gE->SetPoint(i, x_temp, y_temp, prob_temp/max);
  }


  //Fill the TH2F here using the TGraph2D
  for(Int_t i=0; i<gE->GetN(); i++){
    Double_t x_temp = gE->GetX()[i];
    Double_t y_temp = gE->GetY()[i];
    Double_t prob_temp = gE->GetZ()[i];
    //cout << "x_temp = " << x_temp <<", y_temp= " << y_temp <<", prob_temp= " << prob_temp <<endl;
    his2D->Fill(x_temp,y_temp,prob_temp);
    //his2D->Fill(1.0,1.0,1.0);
  } 

  //cout << "gE->GetN()= " << gE->GetN() << ", gE->GetZ()[0]= " << gE->GetZ()[0] << endl;

  //Run Interpolate command to calculate the delaunay triangles, from previous experience we know that otherwise the results can be unstable for the first Interpolate calls. For the actual interpolation we use the TH2F object now, so this step is not critical.  
  Int_t n_interpolate=10;
  Double_t x_test,y_test;
  for(Int_t jj=0;jj<n_interpolate;jj++){
    x_test = HepUniformRand()*getEndpoint(mass_number,atomic_number);
    y_test = HepUniformRand()*getEndpoint(mass_number,atomic_number); 
    if(x_test + y_test < getEndpoint(mass_number,atomic_number)) {
      //for(Int_t hh=0;hh<10;hh++){
      gE->Interpolate(x_test,y_test);
      //cout << "Interpolation # = " << jj << ", x_test= " << x_test << ", y_test= " << y_test << ", gE->Interpolate(x_test,y_test)= " << gE->Interpolate(x_test,y_test) << endl; 
    //}
    }
  }

  infile.close();

}



void fillAngle(Bool_t is0NuBB, TString isotope, TGraph* gE){

  char* datadir = getenv("DBETA_DATA");
  
  //Open up a file - Sum First.
  ifstream infile;
  
  //Get the right data file;
  TString fileName=datadir;
  fileName.Append(Form("/%s",isotope.Data()));
  fileName.Append("correlationsfs");
  if(is0NuBB)  fileName.Append("_0nuBB");
  fileName.Append(".dat");

  //Open the file.
  infile.open(fileName.Data());
  if(!infile){
    cout << "[ERROR] Could not find angle file: " << fileName.Data() << endl;
    return;
  }

  //These should be reset.
  Int_t npoints=0;
  Double_t xp=-1;
  Double_t sum=0;
  Double_t max=0;
  //
  Int_t index;
  Double_t x,dx;
  Double_t y;
  string line;
  while (getline (infile, line)){
    istringstream linestream(line);
    linestream >> index >> x >> y;
    if(xp>0) {
      dx=x-xp;
      sum+=y*dx;
      if(y>max) max=y;
    }
    if(x>0) {
      gE->SetPoint(npoints, x, y);
      npoints++;
    }
    xp=x;
  }
  infile.close();

} 

int main(int argc, char **argv)
{

  Int_t nevent=10;
  Bool_t is0NuBB=false;
  //unsigned Int_t user_sd=0;  

  //Get the number of events.
  if (argc > 1) nevent= atoi(argv[1]);
  //Get the decay type: 0 for 2NuBB or 1 for 0NuBB
  if (argc > 2) is0NuBB= atoi(argv[2]);
  cout << "# Welcome to DBeta " << endl;
  cout << "# Running N_events= " << nevent;
  if(is0NuBB){
    cout << " with 0NuBB. " << endl;
  }else{
    cout << " with 2NuBB. " << endl;
  }

  //Give it a random start...
  time_t seconds;
  seconds = time (NULL);
  HepRandom::getTheEngine()->setSeed(seconds,3); 

  //output files
  //HEPEVT file 
  std::ofstream hepevt_file;
  char* output_dir = getenv("DBETA_OUTPUT");
  TString fileName=output_dir;
  TString outfile = "hepevt.EVT";
  TString fileName2=fileName;
  fileName2.Append(Form("/%s",outfile.Data()));
  hepevt_file.open(fileName2,std::ios::out);

  //open ROOT file to store graphs and histograms
  TString fileName3=fileName;
  outfile = "DBeta_Debug.root";
  fileName3.Append(Form("/%s",outfile.Data()));
  TFile fout(fileName3, "recreate");


  //
  //These are for Debugging.
  //
  TGraph* gSum = new TGraph();
  gSum->SetName("gSum");
  //gSum->GetXaxis()->SetTitle("E1 + E2 [MeV]");
  TGraph* gSingle = new TGraph();
  gSingle->SetName("gSingle");
  //gSingle->GetXaxis()->SetTitle("single electron energy [MeV]");
  TGraph* gAngle = new TGraph();
  gAngle->SetName("gAngle");
  //gAngle->GetXaxis()->SetTitle("cos(#Theta)");
  TGraph2D* g2D = new TGraph2D();
  g2D->SetName("g2D");
  //g2D->GetXaxis()->SetTitle("E1 [MeV]");
  //g2D->GetYaxis()->SetTitle("E2 [MeV]");
  
  //TH2F* his2D;
  TNtuple *tn_dir = new TNtuple("tn_dir","","e1:e2:cos_Theta:momentum_dir_x:momentum_dir_y:momentum_dir_z:momentum_dir2_x:momentum_dir2_y:momentum_dir2_z");

  //
  //Sum is not used but we fill it anyways for debugging.
  //
  
  //get endpoint of decay
  Double_t endpoint = getEndpoint(mass_number,atomic_number); //Cd116

  //here we fill TGraphs. They contain the input data which was obtained from the authors of Phys. Rev. C 85, 034316 (2012). 
  if(is0NuBB){
    gSum->SetPoint(0, endpoint, 1.0);
  }else{
    fillEnergy(true, is0NuBB, isotope_name, gSum);
  }
  fillEnergy(false, is0NuBB, isotope_name, gSingle);
  fillAngle(is0NuBB, isotope_name, gAngle);
  if(!is0NuBB) {
    fillEnergy2D(is0NuBB, isotope_name, g2D);}


  //Histograms for Debugging.
  TH1D* hSum = new TH1D("hSum", ";E1 + E2 [MeV];entries", 100, 0, endpoint);
  TH1D* hSingle = new TH1D("hSingle",";E1 [MeV];entries", 100, 0, endpoint);
  TH1D* hSingle2 = new TH1D("hSingle2",";E2 [MeV];entries", 100, 0, endpoint);
  TH1D* hAngle = new TH1D("hAngle","Angle between the two electrons;cos(#Theta);entries", 100, -1.0, 1.0);
  TH1D* hAngleRot = new TH1D("hAngleRot","Angle between the two electrons after final rotation;cos(#Theta);entries", 100, -1.0, 1.0);

  //These are needed for Hepevt.
  Int_t PDGcode = 11; //electron.
  Double_t melec = electron_mass_c2; //MeV already
 

   
  //Here is the real loop.
  //new method with new input data. For 2NuBB we now have the differential probability as a function of both E1 and E2 (provided by J. Kotila). 
  Double_t e1, e2, esum;
  Double_t r1, r2, r3; // These are the energies.
  Double_t r1s, rsum_check; //these are the norm.
  Double_t a1, a2, phi; //relative angles.
  Double_t rphi, rtheta, rpsi; //final rotation.
  Double_t diff;
  Bool_t doneE, doneE1, doneEs, doneA, doneA1;
  for (Int_t ievent=0; ievent<nevent; ievent++) {    

    cout << "Event number = " << ievent+1 << endl;
    //
    //First the energy of the particles.
    //
    doneE=false;
    while(!doneE){
      //doneE1=false;
      //doneEs=false;
      //while(!doneE1){
      if(!is0NuBB){  //for 2NuBB we use the 2D input data: draw E1 and E2 together, reject with box method, input data is stored in TGraph2D or TH2F
        r1 = HepUniformRand()*endpoint;
	r2 = HepUniformRand()*endpoint;
        r3 = HepUniformRand(); // differential probabilities in TGraph2D have been normalized to 1 at the maximum 
	//if(g2D->Interpolate(r1,r2) >= r3){ //Using the TGraph to interpolate is very slow
        if(his2D->Interpolate(r1,r2) >= r3){ //Using the TH2F to interpolate is faster
          e1=r1;
          e2=r2;
          esum=e1+e2;
	  doneE=true;
	}
      } 
      //} 
      // 0nuBB decay: use single electron spectrum to draw e1
      if(is0NuBB){
        doneE1=false;
        //doneEs=false;
        //first electron.
        while(!doneE1){
        r1 = HepUniformRand()*endpoint;
        r2 = HepUniformRand();
        if(gSingle->Eval(r1) >= r2){
          doneE1=true;
        }
      }
        
      //Is this overly simple, new physics?
      r1s = endpoint - r1;
      
      // Kinematics. (gives distorted sum spectrum). 
      //diff = endpoint - r1s - r1;
      //if(diff >= 0) {
      	doneE=true;
      	e1=r1;
      	esum=r1s+r1;
       	e2=r1s;
      //}
      }
       
    }


    //Debugging.
    hSum->Fill(esum);
    hSingle->Fill(e1);
    hSingle2->Fill(e2);   

 
    //
    //Angle.
    //
    //Start with on going up.
    Hep3Vector momentum_dir(0, 0, 1);
    //Now get the correlation.
    Double_t corrA = gAngle->Eval(e1);
    Double_t funcA = 0;
    doneA=false;
    while(!doneA){
      a1 = 2.0*HepUniformRand()-1.0; //cos_theta from -1 to 1.
      funcA = 1.0 + (corrA*a1);
      //a2 = HepUniformRand(); //shouldn't this be a2 = HepUniformRand()*2.0;  Bug ?
      a2 = HepUniformRand()*2.0;
      if(funcA >= a2){
	doneA=true;
      }
    }

    //Next: Draw random orientation in (x,y) plane
    phi = HepUniformRand()*2*TMath::Pi();

    //new implementation to get truly random directions for the electrons
    //The first electron has z direction and we got cos(Theta)=a1 above --> z=a1 for the second electron direction 
    Hep3Vector momentum_dir2(sqrt(1-a1*a1)*TMath::Cos(phi),sqrt(1-a1*a1)*TMath::Sin(phi),a1); //normalized to length = 1.0

    //cout << "----------------------------------------" << endl;
    //cout << "momentum_dir= " << momentum_dir << endl;
    //cout << "momentum_dir2= " << momentum_dir2 << endl;
    //cout << "phi= " << phi << endl;

    //now we find a random direction for the first electron
    Double_t phi2 = HepUniformRand()*2.0*TMath::Pi();
    Double_t theta2 = TMath::ACos(2.0*HepUniformRand() - 1.0);
    Double_t dirx_rand = TMath::Sin(theta2)*TMath::Cos(phi2);
    Double_t diry_rand = TMath::Sin(theta2)*TMath::Sin(phi2);
    Double_t dirz_rand = TMath::Cos(theta2);  
    momentum_dir.setX(dirx_rand);
    momentum_dir.setY(diry_rand);
    momentum_dir.setZ(dirz_rand);

    //cout << "momentum_dir= " << momentum_dir << endl;

    //find the vector (cross product) and angle (scalar product) to rotate (0,0,1) into the new random vector momentum_dir
    Hep3Vector z_dir(0.,0.,1.);
    Hep3Vector rot_axis = z_dir.cross(momentum_dir);    
    Double_t rot_axis_norm = sqrt(rot_axis.x()*rot_axis.x() + rot_axis.y()*rot_axis.y() + rot_axis.z()*rot_axis.z());
    rot_axis = rot_axis * 1./rot_axis_norm ;
    Double_t rot_angle = TMath::ACos(z_dir.dot(momentum_dir));

    //cout << "rot_axis = " << rot_axis << ", rot_angle= " << rot_angle << endl;

    //debug: 
    //Hep3Vector debug_vec = z_dir.rotate(rot_angle,rot_axis);
    //cout << "debug_vec= " << debug_vec << endl;

    //Finally, apply the same rotation to momentum_dir2
    momentum_dir2.rotate(rot_angle,rot_axis);
    
    //cout << "momentum_dir2= " << momentum_dir2 << endl;

    
    //old implementation: did not give uniformly distributed direction for the two electron directions. 
    /*
    //Here is the second vector.
    Hep3Vector momentum_dir2(momentum_dir.x(),
			     momentum_dir.y(),
			     momentum_dir.z());
    //cout << "----------------------------------------" << endl;
    //cout << "momentum_dir2= " << momentum_dir2 << endl;
    //cout << "phi= " << phi << ", TMath::ACos(a1)= " << TMath::ACos(a1) << endl; 
    momentum_dir2.rotate(phi, TMath::ACos(a1), 0.);  //with the angle definition we use, phi doesn't do anything. 
    //momentum_dir2.rotate(0., TMath::ACos(a1), 0.); //debug 
    //cout << "momentum_dir2= " << momentum_dir2 << endl;
    //cout << "momentum_dir2 components= (" << momentum_dir2.x() << "," << momentum_dir2.y() << "," << momentum_dir2.z() << ")." <<endl;

    //Now final rotation.i

    //cout << "----------------------------------------" << endl;
    //cout << "Before final rotation: momentum_dir = " << momentum_dir << endl;
    //cout << "Before final rotation: momentum_dir2 = " << momentum_dir2 << endl;
    rphi =  HepUniformRand()*2*TMath::Pi();
    rtheta =  HepUniformRand()*2*TMath::Pi();
    rpsi =  HepUniformRand()*2*TMath::Pi();
    //cout << "Rotation angles are: Phi= " << rphi << ", Theta= " << rtheta << ", Psi = " << rpsi <<"." << endl;
    momentum_dir.rotate(rphi, rtheta, rpsi);
    momentum_dir2.rotate(rphi, rtheta, rpsi);
    //cout << "After final rotation: momentum_dir = " << momentum_dir << endl;
    //cout << "After final rotation: momentum_dir2 = " << momentum_dir2 << endl;
    */

    //Debug: fill NTuple with rotated vectors.  
    tn_dir->Fill(e1,e2,a1,momentum_dir.x(),momentum_dir.y(),momentum_dir.z(),momentum_dir2.x(),momentum_dir2.y(),momentum_dir2.z());
   

    //Debug chosen correllation.
    hAngle->Fill(a1);
    //And Check.
    hAngleRot->Fill(momentum_dir.dot(momentum_dir2));

    //Get Ready to Print the Event.
    //cout << e1 << " " << e2 << endl;
    e1 += melec;
    e2 += melec;
    //cout << e1 << " " << e2 << endl;
    Double_t p1=sqrt(e1*e1 - melec*melec);
    Double_t p2=sqrt(e2*e2 - melec*melec);
    //Now convert                                                                                     
    p1*=1e-3; //Gev from MeV
    p2*=1e-3; //GeV from MeV
    melec*=1e-3; //GeV from MeV 
    
    // First, a line for the number of particles in this event.
    hepevt_file << 2 << endl;                     // NHEP
    // Next, line(s) for each particle.  Don't forget spaces between fields.
    hepevt_file << 1 << ' '                       // ISTHEP -- 1 for real particles
	 << PDGcode << ' '                        // IDHEP  -- the PDG code
	 << 0 << ' ' << 0 << ' '                  // daughter indices

	 << momentum_dir.x()*p1 << ' '
	 << momentum_dir.y()*p1 << ' '
	 << momentum_dir.z()*p1 << ' '
         << melec << endl;  
    hepevt_file << 1 << ' '                       // ISTHEP -- 1 for real particles
	 << PDGcode << ' '                        // IDHEP  -- the PDG code
	 << 0 << ' ' << 0 << ' '                  // daughter indices
	 << momentum_dir2.x()*p2 << ' '
	 << momentum_dir2.y()*p2 << ' '
	 << momentum_dir2.z()*p2 << ' '
         << melec << endl;    

    

  }

  //set axis titles for TGraphs
  gSum->GetXaxis()->SetTitle("E1 + E2 [MeV]");
  gSingle->GetXaxis()->SetTitle("single electron energy [MeV]");
  gAngle->GetXaxis()->SetTitle("E1 [MeV]");
  gAngle->GetYaxis()->SetTitle("angular correlation #alpha");
  if(!is0NuBB) {
    g2D->GetXaxis()->SetTitle("E1 [MeV]");
    g2D->GetYaxis()->SetTitle("E2 [MeV]");
    }

  //Finalize I/O. 
  //TFile fout("DBeta_Debug.root", "recreate");
  hSum->Write();
  hSingle->Write();
  hSingle2->Write();
  hAngle->Write();
  gSum->Write();
  gSingle->Write();
  gAngle->Write();
  hAngleRot->Write();
  if(!is0NuBB) {
    g2D->Write(); 
    his2D->Write(); }

  tn_dir->Write();

  fout.Close();
  hepevt_file.close();
 
  return 0;

}


Hep3Vector IsotropicDist()
{
  Double_t u,v,w,r2;
  do {
    u= HepUniformRand()*2.0-1.0;
    v= HepUniformRand()*2.0-1.0;
    w= HepUniformRand()*2.0-1.0;
    r2= u*u+v*v+w*w;
  } while (r2 > 1.0 || r2 < 0.0625);
  r2= sqrt(r2);
  u /= r2;
  v /= r2;
  w /= r2;

  return Hep3Vector(u,v,w);
}
