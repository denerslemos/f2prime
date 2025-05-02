#include "call_libraries.h" // call libraries from ROOT and C++
#include "intro.h"	    // call UIC logo and initialization
#include "ntrkoff.h"        // get Ntrk offline

std::map<unsigned long long, int> runLumiEvtToEntryMap;
unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event);

/*
Main skim f2prime pPb data

Written by Dener Lemos (dener.lemos@cern.ch)

--> Arguments
input_file: text file with a list of root input files: Forest or Skims from jets
input_V0file: text files with a list of V0 files compatible with jets
ouputfile: just a counting number to run on Condor
ntrkoff: multiplicity selection
*/
void f2prime(TString input_file, TString input_V0file, TString ouputfile, int ntrkoff){

	typedef PtEtaPhiMVector LorentzVector;

    const double pimass = 0.13957;
    const double pro_mass = 0.93827;
    const double electronMass = 0.000511;
    
	float ptmin = 0.3;
	float etamin = 2.4;

	TString outputFileName;
	outputFileName = Form("%s",ouputfile.Data());

	clock_t sec_start, sec_end;
	sec_start = clock(); // start timing measurement

	TDatime* date = new TDatime();

	printwelcome(true); // welcome message

	print_start(); // start timing print

	// for track efficiency correction
	// Track or particle efficiency file
	TFile fileeff("DATA_SAMPLES/auxfiles/Hijing_8TeV_dataBS.root","READ");
    	TH2 *eff_factor = nullptr; 
    	fileeff.GetObject("rTotalEff3D_0", eff_factor);  // eff / (1 - fake rate)
	TFile fileefftight("DATA_SAMPLES/auxfiles/Hijing_8TeV_MB_eff_v3_tight.root","READ");
    	TH2 *eff_factortight = nullptr; 
    	fileefftight.GetObject("rTotalEff3D_0", eff_factortight);  // eff / (1 - fake rate)
	TFile fileeffloose("DATA_SAMPLES/auxfiles/Hijing_8TeV_MB_eff_v3_loose.root","READ");
    	TH2 *eff_factorloose = nullptr; 
    	fileeffloose.GetObject("rTotalEff3D_0", eff_factorloose);  // eff / (1 - fake rate)

	TFile *outputFile = new TFile(outputFileName, "RECREATE");

	// Read the input jet file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of Jet input files not founded!" << endl; return;}{cout << "List of Jet input files founded! --> " << input_file.Data() << endl;}
	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(Form("root://osg-se.sprace.org.br/%s",file_chain.c_str()));}
	inputfile.close();
	TChain *heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *skimTree = new TChain("skimanalysis/HltTree");
	TChain *trackTree = new TChain("ppTrack/trackTree");
	TChain *checkFlatteningTree = new TChain("checkflattening/tree");
	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		heavyIonTree->Add(*listIterator);
		skimTree->Add(*listIterator);
		trackTree->Add(*listIterator);
		checkFlatteningTree->Add(*listIterator);
	}

	// Read the input V0 file(s)
	fstream inputfileV0;
	inputfileV0.open(Form("%s",input_V0file.Data()), ios::in);
	if(!inputfileV0.is_open()){cout << "List of V0 input files not founded!" << endl; return;}{cout << "List of V0 input files founded! --> " << input_V0file.Data() << endl;}
	// Make a chain and a vector of file names
	std::vector<TString> file_name_vectorV0;
	string file_chainV0;
	while(getline(inputfileV0, file_chainV0)){file_name_vectorV0.push_back(Form("%s",file_chainV0.c_str()));}
	inputfileV0.close();
	TChain *MainV0Tree = new TChain("K0SAnalysis/my_tree");
	TChain *K0sTree = new TChain("K0SAnalysis/my_treeK0s");

	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vectorV0.begin(); listIterator != file_name_vectorV0.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		MainV0Tree->Add(*listIterator);
		K0sTree->Add(*listIterator);
	}

	// V0 Branches
	TBranch *V0_runBranch;	 // Branch for run
	TBranch *V0_evtBranch;	 // Branch for event
	TBranch *V0_lumiBranch;	 // Branch for lumi
	Int_t V0_run;			 // Run number
	Int_t V0_evt;			 // Event number
	Int_t V0_lumi;			 // Luminosity block

	//---------------------------  K0s  -------------------------
	//daughter 1 (pi+)
	TBranch *K0s_dxy1Branch; 
	TBranch *K0s_dz1Branch; 
	TBranch *K0s_chi21Branch; 
	TBranch *K0s_d1pxBranch; 
	TBranch *K0s_d1pyBranch; 
	TBranch *K0s_d1pzBranch; 
	TBranch *K0s_d1MBranch; 
	TBranch *K0s_d1NhitBranch;
	TBranch *K0s_d1pterrBranch; 

	//daughter 2 (pi-)
	TBranch *K0s_dxy2Branch; 
	TBranch *K0s_dz2Branch; 
	TBranch *K0s_chi22Branch; 
	TBranch *K0s_d2pxBranch; 
	TBranch *K0s_d2pyBranch; 
	TBranch *K0s_d2pzBranch; 
	TBranch *K0s_d2MBranch; 
	TBranch *K0s_d2NhitBranch;
	TBranch *K0s_d2pterrBranch; 

	//Mother (K0s)
	TBranch *K0s_3DaglBranch; 
	TBranch *K0s_3DdlBranch; 
	TBranch *K0s_ptBranch;
	TBranch *K0s_etaBranch; 
	TBranch *K0s_phiBranch; 
	TBranch *K0s_massBranch; 
	TBranch *K0s_dcaBranch; 
	TBranch *K0s_vtxBranch; 	

	// Leaves for K0s particles
	//daughter 1 (pi+)
	vector<double> *K0s_dxy1; 
	vector<double> *K0s_dz1; 
	vector<double> *K0s_chi21; 
	vector<double> *K0s_d1px; 
	vector<double> *K0s_d1py; 
	vector<double> *K0s_d1pz; 
	vector<double> *K0s_d1M; 
	vector<double> *K0s_d1Nhit;
	vector<double> *K0s_d1pterr; 

	//daughter 2 (pi-)
	vector<double> *K0s_dxy2; 
	vector<double> *K0s_dz2; 
	vector<double> *K0s_chi22; 
	vector<double> *K0s_d2px; 
	vector<double> *K0s_d2py; 
	vector<double> *K0s_d2pz; 
	vector<double> *K0s_d2M; 
	vector<double> *K0s_d2Nhit;
	vector<double> *K0s_d2pterr; 

	//Mother (K0s)
	vector<double> *K0s_3Dagl; 
	vector<double> *K0s_3Ddl; 
	vector<double> *K0s_pt; 
	vector<double> *K0s_eta; 
	vector<double> *K0s_phi; 
	vector<double> *K0s_mass; 
	vector<double> *K0s_dca; 
	vector<double> *K0s_vtx; 

	//Main tree (event info)
	MainV0Tree->SetBranchStatus("*",0);
	MainV0Tree->SetBranchStatus("runNumber",1);
	MainV0Tree->SetBranchAddress("runNumber",&V0_run,&V0_runBranch);
	MainV0Tree->SetBranchStatus("evNumber",1);
	MainV0Tree->SetBranchAddress("evNumber",&V0_evt,&V0_evtBranch);
	MainV0Tree->SetBranchStatus("LumiSection",1);
	MainV0Tree->SetBranchAddress("LumiSection",&V0_lumi,&V0_lumiBranch);

	//K0s tree
	K0sTree->SetBranchStatus("*",0);

	K0sTree->SetBranchStatus("K0s_dxy1",1);
	K0sTree->SetBranchAddress("K0s_dxy1",&K0s_dxy1,&K0s_dxy1Branch);
	K0sTree->SetBranchStatus("K0s_dz1",1);
	K0sTree->SetBranchAddress("K0s_dz1",&K0s_dz1,&K0s_dz1Branch);
	K0sTree->SetBranchStatus("K0s_chi21",1);
	K0sTree->SetBranchAddress("K0s_chi21",&K0s_chi21,&K0s_chi21Branch);
	K0sTree->SetBranchStatus("K0s_d1px",1);
	K0sTree->SetBranchAddress("K0s_d1px",&K0s_d1px,&K0s_d1pxBranch);
	K0sTree->SetBranchStatus("K0s_d1py",1);
	K0sTree->SetBranchAddress("K0s_d1py",&K0s_d1py,&K0s_d1pyBranch);
	K0sTree->SetBranchStatus("K0s_d1pz",1);
	K0sTree->SetBranchAddress("K0s_d1pz",&K0s_d1pz,&K0s_d1pzBranch);
	K0sTree->SetBranchStatus("K0s_d1M",1);
	K0sTree->SetBranchAddress("K0s_d1M",&K0s_d1M,&K0s_d1MBranch);
	K0sTree->SetBranchStatus("K0s_d1Nhit",1);
	K0sTree->SetBranchAddress("K0s_d1Nhit",&K0s_d1Nhit,&K0s_d1NhitBranch);
	K0sTree->SetBranchStatus("K0s_d1pterr",1);
	K0sTree->SetBranchAddress("K0s_d1pterr",&K0s_d1pterr,&K0s_d1pterrBranch);	

	K0sTree->SetBranchStatus("K0s_dxy2",1);
	K0sTree->SetBranchAddress("K0s_dxy2",&K0s_dxy2,&K0s_dxy2Branch);
	K0sTree->SetBranchStatus("K0s_dz2",1);
	K0sTree->SetBranchAddress("K0s_dz2",&K0s_dz2,&K0s_dz2Branch);
	K0sTree->SetBranchStatus("K0s_chi22",1);
	K0sTree->SetBranchAddress("K0s_chi22",&K0s_chi22,&K0s_chi22Branch);
	K0sTree->SetBranchStatus("K0s_d2px",1);
	K0sTree->SetBranchAddress("K0s_d2px",&K0s_d2px,&K0s_d2pxBranch);
	K0sTree->SetBranchStatus("K0s_d2py",1);
	K0sTree->SetBranchAddress("K0s_d2py",&K0s_d2py,&K0s_d2pyBranch);
	K0sTree->SetBranchStatus("K0s_d2pz",1);
	K0sTree->SetBranchAddress("K0s_d2pz",&K0s_d2pz,&K0s_d2pzBranch);
	K0sTree->SetBranchStatus("K0s_d2M",1);
	K0sTree->SetBranchAddress("K0s_d2M",&K0s_d2M,&K0s_d2MBranch);
	K0sTree->SetBranchStatus("K0s_d2Nhit",1);
	K0sTree->SetBranchAddress("K0s_d2Nhit",&K0s_d2Nhit,&K0s_d2NhitBranch);
	K0sTree->SetBranchStatus("K0s_d2pterr",1);
	K0sTree->SetBranchAddress("K0s_d2pterr",&K0s_d2pterr,&K0s_d2pterrBranch);	

	K0sTree->SetBranchStatus("K0s_3Dagl",1);
	K0sTree->SetBranchAddress("K0s_3Dagl",&K0s_3Dagl,&K0s_3DaglBranch);	
	K0sTree->SetBranchStatus("K0s_3Ddl",1);
	K0sTree->SetBranchAddress("K0s_3Ddl",&K0s_3Ddl,&K0s_3DdlBranch);	
	K0sTree->SetBranchStatus("K0s_pt",1);
	K0sTree->SetBranchAddress("K0s_pt",&K0s_pt,&K0s_ptBranch);	
	K0sTree->SetBranchStatus("K0s_eta",1);
	K0sTree->SetBranchAddress("K0s_eta",&K0s_eta,&K0s_etaBranch);	
	K0sTree->SetBranchStatus("K0s_phi",1);
	K0sTree->SetBranchAddress("K0s_phi",&K0s_phi,&K0s_phiBranch);	
	K0sTree->SetBranchStatus("K0s_mass",1);
	K0sTree->SetBranchAddress("K0s_mass",&K0s_mass,&K0s_massBranch);	
	K0sTree->SetBranchStatus("K0s_dca",1);
	K0sTree->SetBranchAddress("K0s_dca",&K0s_dca,&K0s_dcaBranch);	
	K0sTree->SetBranchStatus("K0s_vtx",1);
	K0sTree->SetBranchAddress("K0s_vtx",&K0s_vtx,&K0s_vtxBranch);	

	// Bellow here is for jets
	// Branches for heavy ion tree
	TBranch *runBranch;						 	// Branch for run
	TBranch *eventBranch;					 	// Branch for event
	TBranch *lumiBranch;						// Branch for lumi
	TBranch *hiVzBranch;						// Branch for vertex z-position
	TBranch *hiVxBranch;						// Branch for vertex x-position
	TBranch *hiVyBranch;						// Branch for vertex y-position
	TBranch *hiHFplusBranch;					// Branch for HF+ energy deposity
	TBranch *hiHFminusBranch;					// Branch for HF- energy deposity
	TBranch *hiHFplusEta4Branch;				// Branch for HF+ energy deposity for |eta| > 4
	TBranch *hiHFminusEta4Branch;				// Branch for HF- energy deposity for |eta| > 4
	
	// Leaves for heavy ion tree
	UInt_t run;					 // Run number
	ULong64_t event;			 // Event number
	UInt_t lumi;				 // Luminosity block
	Float_t vertexZ;			 // Vertex z-position
	Float_t vertexX;			 // Vertex x-position
	Float_t vertexY;			 // Vertex y-position
	Float_t hiHFplus;			 // transverse energy sum of HF+ tower;
	Float_t hiHFminus;			 // transverse energy sum of HF- tower;
	Float_t hiHFplusEta4;		 // transverse energy sum of HF+ tower for |eta| > 4;
	Float_t hiHFminusEta4;		 // transverse energy sum of HF- tower for |eta| > 4;
	Int_t Ntroff;		 		 // Multiplicity --> Ntrkoffline
	Float_t Nch;		 		 	 // Multiplicity --> Corrected
	Float_t NchLoose;		 		 // Multiplicity --> Corrected Loose
	Float_t NchTight;		 		 // Multiplicity --> Corrected Loose


	// Branches for EP tree
	TBranch *eventPlaneAngleBranch;						// Branch for event plane angles
	TBranch *eventPLaneQBranch;							// Branch for Q-vector magnitude in an event plane
	TBranch *eventPlaneQxBranch;						// Branch for Q-vector x-component in an event plane
	TBranch *eventPlaneQyBranch;						// Branch for Q-vector y-component in an event plane
	TBranch *eventPlaneMultiplicityBranch;	 			// Branch for event plane multiplicity
	
	const int numberofEPleaves = 182;				    		// Event plane leaves	
	// Name of EP in pPb
	TString EPNames = "HFm1/D:HFp1/D:trackmid1/D:trackm1/D:trackp1/D:trackm122/D:trackm118/D:trackm114/D:trackm110/D:trackm106/D:trackm102/D:trackp102/D:trackp106/D:trackp110/D:trackp114/D:trackp118/D:trackp122/D:trackmid1mc/D:trackm1mc/D:trackp1mc/D:trackm122mc/D:trackm118mc/D:trackm114mc/D:trackm110mc/D:trackm106mc/D:trackm102mc/D:trackp102mc/D:trackp106mc/D:trackp110mc/D:trackp114mc/D:trackp118mc/D:trackp122mc/D:HFm1a/D:HFm1b/D:HFm1c/D:HFm1d/D:HFm1e/D:HFm1f/D:HFp1a/D:HFp1b/D:HFp1c/D:HFp1d/D:HFp1e/D:HFp1f/D:HFm2/D:HFp2/D:trackmid2/D:trackm2/D:trackp2/D:trackm222/D:trackm218/D:trackm214/D:trackm210/D:trackm206/D:trackm202/D:trackp202/D:trackp206/D:trackp210/D:trackp214/D:trackp218/D:trackp222/D:HFm2a/D:HFm2b/D:HFm2c/D:HFm2d/D:HFm2e/D:HFm2f/D:HFp2a/D:HFp2b/D:HFp2c/D:HFp2d/D:HFp2e/D:HFp2f/D:HFm3/D:HFp3/D:trackmid3/D:trackm3/D:trackp3/D:trackm322/D:trackm318/D:trackm314/D:trackm310/D:trackm306/D:trackm302/D:trackp302/D:trackp306/D:trackp310/D:trackp314/D:trackp318/D:trackp322/D:HFm3a/D:HFm3b/D:HFm3c/D:HFm3d/D:HFm3e/D:HFm3f/D:HFp3a/D:HFp3b/D:HFp3c/D:HFp3d/D:HFp3e/D:HFp3f/D:HFm4/D:HFp4/D:trackmid4/D:trackm4/D:trackp4/D:trackm422/D:trackm418/D:trackm414/D:trackm410/D:trackm406/D:trackm402/D:trackp402/D:trackp406/D:trackp410/D:trackp414/D:trackp418/D:trackp422/D:HFm4a/D:HFm4b/D:HFm4c/D:HFm43d/D:HFm4e/D:HFm4f/D:HFp4a/D:HFp4b/D:HFp4c/D:HFp4d/D:HFp4e/D:HFp4f/D:HFm5/D:HFp5/D:trackmid5/D:trackm5/D:trackp5/D:trackm522/D:trackm518/D:trackm514/D:trackm510/D:trackm506/D:trackm502/D:trackp502/D:trackp506/D:trackp510/D:trackp514/D:trackp518/D:trackp522/D:HFm6/D:HFp6/D:trackmid6/D:trackm6/D:trackp6/D:trackm622/D:trackm618/D:trackm614/D:trackm610/D:trackm606/D:trackm602/D:trackp602/D:trackp606/D:trackp610/D:trackp614/D:trackp618/D:trackp622/D:HFm7/D:HFp7/D:trackmid7/D:trackm7/D:trackp7/D:trackm722/D:trackm718/D:trackm714/D:trackm710/D:trackm706/D:trackm702/D:trackp702/D:trackp706/D:trackp710/D:trackp714/D:trackp718/D:trackp722/D"; 
	Double_t eventPlaneAngle[numberofEPleaves] = {0};			// Event plane angles
	Double_t eventPlaneQ[numberofEPleaves] = {0};				// Magnitude of Q-vector in event plane
	Double_t eventPlaneQx[numberofEPleaves] = {0};				// x-component of the Q-vector
	Double_t eventPlaneQy[numberofEPleaves] = {0};				// y-component of the Q-vector
	Double_t eventPlaneMultiplicity[numberofEPleaves] = {0};	// Particle multiplicity in an event plane

	// Branches for skim tree
	TBranch *primaryVertexBranch;						// Branch for primary vertex filter bit
	TBranch *beamScrapingBranch;				// Branch for beam scraping filter bit
	TBranch *hfCoincidenceBranch;				// Branch for energy recorded one HF tower above threshold on each side
	TBranch *pVertexFilterCutdz1p0Branch;		// Branch for PU Filter default
	TBranch *pVertexFilterCutGplusBranch;		// Branch for PU Filter GPlus
	TBranch *pVertexFilterCutVtx1Branch;		// Branch for PU Filter 1 vertex only

	// Leaves for the skim tree
	Int_t primaryVertexFilterBit;				// Filter bit for primary vertex
	Int_t beamScrapingFilterBit;				// Filter bit for beam scraping
	Int_t hfCoincidenceFilterBit;				// Filter bit or energy recorded one HF tower above threshold on each side
	Int_t pVertexFilterCutdz1p0Bit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutGplusBit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutVtx1Bit;					// Filter bit for PU Filter

	// Branches for track tree
	TBranch *nTracksBranch;									// Branch for number of tracks
	TBranch *trackPtBranch;									// Branch for track pT
	TBranch *trackPtErrorBranch;							// Branch for track pT error
	TBranch *trackPhiBranch;								// Branch for track phi
	TBranch *trackEtaBranch;								// Branch for track eta
	TBranch *trackHighPurityBranch;							// Branch for high purity of the track
	TBranch *trackVertexDistanceZBranch;			 		// Branch for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceZErrorBranch;				// Branch for error for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceXYBranch;					// Branch for track distance from primary vertex in xy-direction
	TBranch *trackVertexDistanceXYErrorBranch; 				// Branch for error for track distance from primary vertex in xy-direction
	TBranch *trackChargeBranch;								// Branch for track charge
	
	// Leaves for the track tree
	const Int_t nMaxTrack = 2000;
	Int_t nTracks;														// Number of tracks
	Float_t trackPtArray[nMaxTrack] = {0};								// Array for track pT
	Float_t trackPtErrorArray[nMaxTrack] = {0};							// Array for track pT errors
	Float_t trackPhiArray[nMaxTrack] = {0};								// Array for track phis
	Float_t trackEtaArray[nMaxTrack] = {0};								// Array for track etas
	Bool_t trackHighPurityArray[nMaxTrack] = {0};						// Array for the high purity of tracks
	Float_t trackVertexDistanceZArray[nMaxTrack] = {0};			 		// Array for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceZErrorArray[nMaxTrack] = {0};			// Array for error for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceXYArray[nMaxTrack] = {0};				// Array for track distance from primary vertex in xy-direction
	Float_t trackVertexDistanceXYErrorArray[nMaxTrack] = {0}; 			// Array for error for track distance from primary vertex in xy-direction
	Int_t trackChargeArray[nMaxTrack] = {0}; 										// Array for track charge


	// ========================================== //
	// Read all the branches from the input trees //
	// ========================================== //
	
	// Connect the branches of the heavy ion tree
	heavyIonTree->SetBranchStatus("*",0); // remove all branchs to read it fast
	heavyIonTree->SetBranchStatus("run",1);
	heavyIonTree->SetBranchAddress("run",&run,&runBranch);
	heavyIonTree->SetBranchStatus("evt",1);
	heavyIonTree->SetBranchAddress("evt",&event,&eventBranch);
	heavyIonTree->SetBranchStatus("lumi",1);
	heavyIonTree->SetBranchAddress("lumi",&lumi,&lumiBranch);
	heavyIonTree->SetBranchStatus("vz",1);
	heavyIonTree->SetBranchAddress("vz",&vertexZ,&hiVzBranch);
	heavyIonTree->SetBranchStatus("vx",1);
	heavyIonTree->SetBranchAddress("vx",&vertexX,&hiVxBranch);
	heavyIonTree->SetBranchStatus("vy",1);
	heavyIonTree->SetBranchAddress("vy",&vertexY,&hiVyBranch);
	heavyIonTree->SetBranchStatus("hiHFplus",1);
	heavyIonTree->SetBranchAddress("hiHFplus",&hiHFplus,&hiHFplusBranch);
	heavyIonTree->SetBranchStatus("hiHFminus",1);
	heavyIonTree->SetBranchAddress("hiHFminus",&hiHFminus,&hiHFminusBranch);
	heavyIonTree->SetBranchStatus("hiHFplusEta4",1);
	heavyIonTree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4,&hiHFplusEta4Branch);
	heavyIonTree->SetBranchStatus("hiHFminusEta4",1);
	heavyIonTree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4,&hiHFminusEta4Branch);
	
	// Event plane
	checkFlatteningTree->SetBranchStatus("epang",1);
	checkFlatteningTree->SetBranchAddress("epang",&eventPlaneAngle,&eventPlaneAngleBranch);
	checkFlatteningTree->SetBranchStatus("q",1);
	checkFlatteningTree->SetBranchAddress("q",&eventPlaneQ,&eventPLaneQBranch);
	checkFlatteningTree->SetBranchStatus("qx",1);
	checkFlatteningTree->SetBranchAddress("qx",&eventPlaneQx,&eventPlaneQxBranch);
	checkFlatteningTree->SetBranchStatus("qy",1);
	checkFlatteningTree->SetBranchAddress("qy",&eventPlaneQy,&eventPlaneQyBranch);
	checkFlatteningTree->SetBranchStatus("mult",1);
	checkFlatteningTree->SetBranchAddress("mult",&eventPlaneMultiplicity,&eventPlaneMultiplicityBranch);

	// Connect the branches to the skim tree
	skimTree->SetBranchStatus("*",0);
	skimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
	skimTree->SetBranchAddress("pPAprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
	skimTree->SetBranchStatus("pBeamScrapingFilter",1);
	skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapingFilterBit,&beamScrapingBranch);
	skimTree->SetBranchStatus("phfCoincFilter",1);
	skimTree->SetBranchAddress("phfCoincFilter", &hfCoincidenceFilterBit, &hfCoincidenceBranch);
	skimTree->SetBranchStatus("pVertexFilterCutdz1p0",1);
	skimTree->SetBranchAddress("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, &pVertexFilterCutdz1p0Branch);
	skimTree->SetBranchStatus("pVertexFilterCutGplus",1);
	skimTree->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplusBit, &pVertexFilterCutGplusBranch);
	skimTree->SetBranchStatus("pVertexFilterCutVtx1",1);
	skimTree->SetBranchAddress("pVertexFilterCutVtx1", &pVertexFilterCutVtx1Bit, &pVertexFilterCutVtx1Branch);

	// Connect the branches to the track tree
	
	trackTree->SetBranchStatus("*",0);
	trackTree->SetBranchStatus("nTrk",1);
	trackTree->SetBranchAddress("nTrk",&nTracks,&nTracksBranch);
	trackTree->SetBranchStatus("highPurity",1);
	trackTree->SetBranchAddress("highPurity",&trackHighPurityArray,&trackHighPurityBranch);
	trackTree->SetBranchStatus("trkPt",1);
	trackTree->SetBranchAddress("trkPt",&trackPtArray,&trackPtBranch);
	trackTree->SetBranchStatus("trkPtError",1);
	trackTree->SetBranchAddress("trkPtError",&trackPtErrorArray,&trackPtErrorBranch);
	trackTree->SetBranchStatus("trkPhi",1);
	trackTree->SetBranchAddress("trkPhi",&trackPhiArray,&trackPhiBranch);
	trackTree->SetBranchStatus("trkEta",1);
	trackTree->SetBranchAddress("trkEta",&trackEtaArray,&trackEtaBranch);
	trackTree->SetBranchStatus("trkDz1",1);
	trackTree->SetBranchAddress("trkDz1",&trackVertexDistanceZArray,&trackVertexDistanceZBranch);
	trackTree->SetBranchStatus("trkDzError1",1);
	trackTree->SetBranchAddress("trkDzError1",&trackVertexDistanceZErrorArray,&trackVertexDistanceZErrorBranch);
	trackTree->SetBranchStatus("trkDxy1",1);
	trackTree->SetBranchAddress("trkDxy1",&trackVertexDistanceXYArray,&trackVertexDistanceXYBranch);
	trackTree->SetBranchStatus("trkDxyError1",1);
	trackTree->SetBranchAddress("trkDxyError1",&trackVertexDistanceXYErrorArray,&trackVertexDistanceXYErrorBranch);
	trackTree->SetBranchStatus("trkCharge",1);
	trackTree->SetBranchAddress("trkCharge",&trackChargeArray,&trackChargeBranch);
	


	// ========================================== //
	//			 Define output trees
	// ========================================== //

	//--------------- V0s ---------------
  	TTree *F2PrimeTreeOutput = new TTree("F2PrimeTree","");
	std::vector<float> *f2prime_PtVector   = new std::vector<float>(); f2prime_PtVector->clear();
	std::vector<float> *f2prime_EtaVector  = new std::vector<float>(); f2prime_EtaVector->clear();
	std::vector<float> *f2prime_PhiVector  = new std::vector<float>(); f2prime_PhiVector->clear();
	std::vector<float> *f2prime_MassVector = new std::vector<float>(); f2prime_MassVector->clear();

	std::vector<int> *K0s1_ShareDauVector = new std::vector<int>(); K0s1_ShareDauVector->clear();	
	std::vector<int> *K0s1_NominalVector = new std::vector<int>(); K0s1_NominalVector->clear();
	std::vector<int> *K0s1_TightVector = new std::vector<int>(); K0s1_TightVector->clear();
	std::vector<int> *K0s1_LooseVector = new std::vector<int>(); K0s1_LooseVector->clear();
	std::vector<float> *K0s1_MisIDEEVector = new std::vector<float>(); K0s1_MisIDEEVector->clear();
	std::vector<float> *K0s1_MisIDPiProtonVector = new std::vector<float>(); K0s1_MisIDPiProtonVector->clear();
	std::vector<float> *K0s1_MisIDProtonPiVector = new std::vector<float>(); K0s1_MisIDProtonPiVector->clear();
	std::vector<float> *K0s1_APAlphaVector = new std::vector<float>(); K0s1_APAlphaVector->clear();
	std::vector<float> *K0s1_APQTVector = new std::vector<float>(); K0s1_APQTVector->clear();
	std::vector<float> *K0s1_PtVector = new std::vector<float>(); K0s1_PtVector->clear();
	std::vector<float> *K0s1_EtaVector = new std::vector<float>(); K0s1_EtaVector->clear();
	std::vector<float> *K0s1_PhiVector = new std::vector<float>(); K0s1_PhiVector->clear();
	std::vector<float> *K0s1_MassVector = new std::vector<float>(); K0s1_MassVector->clear();
	std::vector<float> *K0s1_D1_PtVector = new std::vector<float>(); K0s1_D1_PtVector->clear();
	std::vector<float> *K0s1_D1_EtaVector = new std::vector<float>(); K0s1_D1_EtaVector->clear();
	std::vector<float> *K0s1_D1_PhiVector = new std::vector<float>(); K0s1_D1_PhiVector->clear();
	std::vector<float> *K0s1_D2_PtVector = new std::vector<float>(); K0s1_D2_PtVector->clear();
	std::vector<float> *K0s1_D2_EtaVector = new std::vector<float>(); K0s1_D2_EtaVector->clear();
	std::vector<float> *K0s1_D2_PhiVector = new std::vector<float>(); K0s1_D2_PhiVector->clear();
	
	std::vector<int> *K0s2_ShareDauVector = new std::vector<int>(); K0s2_ShareDauVector->clear();	
	std::vector<int> *K0s2_NominalVector = new std::vector<int>(); K0s2_NominalVector->clear();
	std::vector<int> *K0s2_TightVector = new std::vector<int>(); K0s2_TightVector->clear();
	std::vector<int> *K0s2_LooseVector = new std::vector<int>(); K0s2_LooseVector->clear();
	std::vector<float> *K0s2_MisIDEEVector = new std::vector<float>(); K0s2_MisIDEEVector->clear();
	std::vector<float> *K0s2_MisIDPiProtonVector = new std::vector<float>(); K0s2_MisIDPiProtonVector->clear();
	std::vector<float> *K0s2_MisIDProtonPiVector = new std::vector<float>(); K0s2_MisIDProtonPiVector->clear();
	std::vector<float> *K0s2_APAlphaVector = new std::vector<float>(); K0s2_APAlphaVector->clear();
	std::vector<float> *K0s2_APQTVector = new std::vector<float>(); K0s2_APQTVector->clear();
	std::vector<float> *K0s2_PtVector = new std::vector<float>(); K0s2_PtVector->clear();
	std::vector<float> *K0s2_EtaVector = new std::vector<float>(); K0s2_EtaVector->clear();
	std::vector<float> *K0s2_PhiVector = new std::vector<float>(); K0s2_PhiVector->clear();
	std::vector<float> *K0s2_MassVector = new std::vector<float>(); K0s2_MassVector->clear();	
	std::vector<float> *K0s2_D1_PtVector = new std::vector<float>(); K0s2_D1_PtVector->clear();
	std::vector<float> *K0s2_D1_EtaVector = new std::vector<float>(); K0s2_D1_EtaVector->clear();
	std::vector<float> *K0s2_D1_PhiVector = new std::vector<float>(); K0s2_D1_PhiVector->clear();
	std::vector<float> *K0s2_D2_PtVector = new std::vector<float>(); K0s2_D2_PtVector->clear();
	std::vector<float> *K0s2_D2_EtaVector = new std::vector<float>(); K0s2_D2_EtaVector->clear();
	std::vector<float> *K0s2_D2_PhiVector = new std::vector<float>(); K0s2_D2_PhiVector->clear();

	F2PrimeTreeOutput->Branch("f2prime_pt","vector<float>", &f2prime_PtVector);
	F2PrimeTreeOutput->Branch("f2prime_eta","vector<float>", &f2prime_EtaVector);
	F2PrimeTreeOutput->Branch("f2prime_phi","vector<float>", &f2prime_PhiVector);
	F2PrimeTreeOutput->Branch("f2prime_mass","vector<float>", &f2prime_MassVector);

	F2PrimeTreeOutput->Branch("K0s1_sharedaughter","vector<int>", &K0s1_ShareDauVector);
	F2PrimeTreeOutput->Branch("K0s1_nominal","vector<int>", &K0s1_NominalVector);
	F2PrimeTreeOutput->Branch("K0s1_tight","vector<int>", &K0s1_TightVector);
	F2PrimeTreeOutput->Branch("K0s1_loose","vector<int>", &K0s1_LooseVector);
	F2PrimeTreeOutput->Branch("K0s1_misidee","vector<float>", &K0s1_MisIDEEVector);
	F2PrimeTreeOutput->Branch("K0s1_misidpipro","vector<float>", &K0s1_MisIDPiProtonVector);
	F2PrimeTreeOutput->Branch("K0s1_misidpropi","vector<float>", &K0s1_MisIDProtonPiVector);
	F2PrimeTreeOutput->Branch("K0s1_alpha","vector<float>", &K0s1_APAlphaVector);
	F2PrimeTreeOutput->Branch("K0s1_qt","vector<float>", &K0s1_APQTVector);
	F2PrimeTreeOutput->Branch("K0s1_pt","vector<float>", &K0s1_PtVector);
	F2PrimeTreeOutput->Branch("K0s1_eta","vector<float>", &K0s1_EtaVector);
	F2PrimeTreeOutput->Branch("K0s1_phi","vector<float>", &K0s1_PhiVector);
	F2PrimeTreeOutput->Branch("K0s1_mass","vector<float>", &K0s1_MassVector);
	F2PrimeTreeOutput->Branch("K0s1_d1_pt","vector<float>", &K0s1_D1_PtVector);
	F2PrimeTreeOutput->Branch("K0s1_d1_eta","vector<float>", &K0s1_D1_EtaVector);
	F2PrimeTreeOutput->Branch("K0s1_d1_phi","vector<float>", &K0s1_D1_PhiVector);
	F2PrimeTreeOutput->Branch("K0s1_d2_pt","vector<float>", &K0s1_D2_PtVector);
	F2PrimeTreeOutput->Branch("K0s1_d2_eta","vector<float>", &K0s1_D2_EtaVector);
	F2PrimeTreeOutput->Branch("K0s1_d2_phi","vector<float>", &K0s1_D2_PhiVector);
	
	F2PrimeTreeOutput->Branch("K0s2_sharedaughter","vector<int>", &K0s2_ShareDauVector);	
	F2PrimeTreeOutput->Branch("K0s2_nominal","vector<int>", &K0s2_NominalVector);
	F2PrimeTreeOutput->Branch("K0s2_tight","vector<int>", &K0s2_TightVector);
	F2PrimeTreeOutput->Branch("K0s2_loose","vector<int>", &K0s2_LooseVector);
	F2PrimeTreeOutput->Branch("K0s2_misidee","vector<float>", &K0s2_MisIDEEVector);
	F2PrimeTreeOutput->Branch("K0s2_misidpipro","vector<float>", &K0s2_MisIDPiProtonVector);
	F2PrimeTreeOutput->Branch("K0s2_misidpropi","vector<float>", &K0s2_MisIDProtonPiVector);
	F2PrimeTreeOutput->Branch("K0s2_alpha","vector<float>", &K0s2_APAlphaVector);
	F2PrimeTreeOutput->Branch("K0s2_qt","vector<float>", &K0s2_APQTVector);
	F2PrimeTreeOutput->Branch("K0s2_pt","vector<float>", &K0s2_PtVector);
	F2PrimeTreeOutput->Branch("K0s2_eta","vector<float>", &K0s2_EtaVector);
	F2PrimeTreeOutput->Branch("K0s2_phi","vector<float>", &K0s2_PhiVector);
	F2PrimeTreeOutput->Branch("K0s2_mass","vector<float>", &K0s2_MassVector);
	F2PrimeTreeOutput->Branch("K0s2_d1_pt","vector<float>", &K0s2_D1_PtVector);
	F2PrimeTreeOutput->Branch("K0s2_d1_eta","vector<float>", &K0s2_D1_EtaVector);
	F2PrimeTreeOutput->Branch("K0s2_d1_phi","vector<float>", &K0s2_D1_PhiVector);
	F2PrimeTreeOutput->Branch("K0s2_d2_pt","vector<float>", &K0s2_D2_PtVector);
	F2PrimeTreeOutput->Branch("K0s2_d2_eta","vector<float>", &K0s2_D2_EtaVector);
	F2PrimeTreeOutput->Branch("K0s2_d2_phi","vector<float>", &K0s2_D2_PhiVector);
	
	// Event plane
	TTree *checkFlatteningTreeOutput = new TTree("tree","");
	float epang_HFm2, epang_HFp2,epang_HFm3, epang_HFp3, epang_HFm4, epang_HFp4, epang_trackmid2, epang_trackmid3, epang_trackmid4, epang_trackm2, epang_trackm3, epang_trackm4, epang_trackp2, epang_trackp3, epang_trackp4;	
	float q_HFm2, q_HFp2,q_HFm3, q_HFp3, q_HFm4, q_HFp4, q_trackmid2, q_trackmid3, q_trackmid4, q_trackm2, q_trackm3, q_trackm4, q_trackp2, q_trackp3, q_trackp4;
	float qx_HFm2, qx_HFp2,qx_HFm3, qx_HFp3, qx_HFm4, qx_HFp4, qx_trackmid2, qx_trackmid3, qx_trackmid4, qx_trackm2, qx_trackm3, qx_trackm4, qx_trackp2, qx_trackp3, qx_trackp4;
	float qy_HFm2, qy_HFp2,qy_HFm3, qy_HFp3, qy_HFm4, qy_HFp4, qy_trackmid2, qy_trackmid3, qy_trackmid4, qy_trackm2, qy_trackm3, qy_trackm4, qy_trackp2, qy_trackp3, qy_trackp4;
	float mult_HFm2, mult_HFp2,mult_HFm3, mult_HFp3, mult_HFm4, mult_HFp4, mult_trackmid2, mult_trackmid3, mult_trackmid4, mult_trackm2, mult_trackm3, mult_trackm4, mult_trackp2, mult_trackp3, mult_trackp4;

	checkFlatteningTreeOutput->Branch("epang_HFm2",&epang_HFm2,"epang_HFm2/F");
	checkFlatteningTreeOutput->Branch("epang_HFp2",&epang_HFp2,"epang_HFp2/F");
	checkFlatteningTreeOutput->Branch("epang_trackmid2",&epang_trackmid2,"epang_trackmid2/F");
	checkFlatteningTreeOutput->Branch("epang_trackm2",&epang_trackm2,"epang_trackm2/F");
	checkFlatteningTreeOutput->Branch("epang_trackp2",&epang_trackp2,"epang_trackp2/F");
	checkFlatteningTreeOutput->Branch("epang_HFm3",&epang_HFm3,"epang_HFm3/F");
	checkFlatteningTreeOutput->Branch("epang_HFp3",&epang_HFp3,"epang_HFp3/F");
	checkFlatteningTreeOutput->Branch("epang_trackmid3",&epang_trackmid3,"epang_trackmid3/F");
	checkFlatteningTreeOutput->Branch("epang_trackm3",&epang_trackm3,"epang_trackm3/F");
	checkFlatteningTreeOutput->Branch("epang_trackp3",&epang_trackp3,"epang_trackp3/F");
	checkFlatteningTreeOutput->Branch("epang_HFm4",&epang_HFm4,"epang_HFm4/F");
	checkFlatteningTreeOutput->Branch("epang_HFp4",&epang_HFp4,"epang_HFp4/F");
	checkFlatteningTreeOutput->Branch("epang_trackmid4",&epang_trackmid4,"epang_trackmid4/F");
	checkFlatteningTreeOutput->Branch("epang_trackm4",&epang_trackm4,"epang_trackm4/F");
	checkFlatteningTreeOutput->Branch("epang_trackp4",&epang_trackp4,"epang_trackp4/F");

	checkFlatteningTreeOutput->Branch("q_HFm2",&q_HFm2,"q_HFm2/F");
	checkFlatteningTreeOutput->Branch("q_HFp2",&q_HFp2,"q_HFp2/F");
	checkFlatteningTreeOutput->Branch("q_trackmid2",&q_trackmid2,"q_trackmid2/F");
	checkFlatteningTreeOutput->Branch("q_trackm2",&q_trackm2,"q_trackm2/F");
	checkFlatteningTreeOutput->Branch("q_trackp2",&q_trackp2,"q_trackp2/F");
	checkFlatteningTreeOutput->Branch("q_HFm3",&q_HFm3,"q_HFm3/F");
	checkFlatteningTreeOutput->Branch("q_HFp3",&q_HFp3,"q_HFp3/F");
	checkFlatteningTreeOutput->Branch("q_trackmid3",&q_trackmid3,"q_trackmid3/F");
	checkFlatteningTreeOutput->Branch("q_trackm3",&q_trackm3,"q_trackm3/F");
	checkFlatteningTreeOutput->Branch("q_trackp3",&q_trackp3,"q_trackp3/F");
	checkFlatteningTreeOutput->Branch("q_HFm4",&q_HFm4,"q_HFm4/F");
	checkFlatteningTreeOutput->Branch("q_HFp4",&q_HFp4,"q_HFp4/F");
	checkFlatteningTreeOutput->Branch("q_trackmid4",&q_trackmid4,"q_trackmid4/F");
	checkFlatteningTreeOutput->Branch("q_trackm4",&q_trackm4,"q_trackm4/F");
	checkFlatteningTreeOutput->Branch("q_trackp4",&q_trackp4,"q_trackp4/F");

	checkFlatteningTreeOutput->Branch("qx_HFm2",&qx_HFm2,"qx_HFm2/F");
	checkFlatteningTreeOutput->Branch("qx_HFp2",&qx_HFp2,"qx_HFp2/F");
	checkFlatteningTreeOutput->Branch("qx_trackmid2",&qx_trackmid2,"qx_trackmid2/F");
	checkFlatteningTreeOutput->Branch("qx_trackm2",&qx_trackm2,"qx_trackm2/F");
	checkFlatteningTreeOutput->Branch("qx_trackp2",&qx_trackp2,"qx_trackp2/F");
	checkFlatteningTreeOutput->Branch("qx_HFm3",&qx_HFm3,"qx_HFm3/F");
	checkFlatteningTreeOutput->Branch("qx_HFp3",&qx_HFp3,"qx_HFp3/F");
	checkFlatteningTreeOutput->Branch("qx_trackmid3",&qx_trackmid3,"qx_trackmid3/F");
	checkFlatteningTreeOutput->Branch("qx_trackm3",&qx_trackm3,"qx_trackm3/F");
	checkFlatteningTreeOutput->Branch("qx_trackp3",&qx_trackp3,"qx_trackp3/F");
	checkFlatteningTreeOutput->Branch("qx_HFm4",&qx_HFm4,"qx_HFm4/F");
	checkFlatteningTreeOutput->Branch("qx_HFp4",&qx_HFp4,"qx_HFp4/F");
	checkFlatteningTreeOutput->Branch("qx_trackmid4",&qx_trackmid4,"qx_trackmid4/F");
	checkFlatteningTreeOutput->Branch("qx_trackm4",&qx_trackm4,"qx_trackm4/F");
	checkFlatteningTreeOutput->Branch("qx_trackp4",&qx_trackp4,"qx_trackp4/F");
	
	checkFlatteningTreeOutput->Branch("qy_HFm2",&qy_HFm2,"qy_HFm2/F");
	checkFlatteningTreeOutput->Branch("qy_HFp2",&qy_HFp2,"qy_HFp2/F");
	checkFlatteningTreeOutput->Branch("qy_trackmid2",&qy_trackmid2,"qy_trackmid2/F");
	checkFlatteningTreeOutput->Branch("qy_trackm2",&qy_trackm2,"qy_trackm2/F");
	checkFlatteningTreeOutput->Branch("qy_trackp2",&qy_trackp2,"qy_trackp2/F");
	checkFlatteningTreeOutput->Branch("qy_HFm3",&qy_HFm3,"qy_HFm3/F");
	checkFlatteningTreeOutput->Branch("qy_HFp3",&qy_HFp3,"qy_HFp3/F");
	checkFlatteningTreeOutput->Branch("qy_trackmid3",&qy_trackmid3,"qy_trackmid3/F");
	checkFlatteningTreeOutput->Branch("qy_trackm3",&qy_trackm3,"qy_trackm3/F");
	checkFlatteningTreeOutput->Branch("qy_trackp3",&qy_trackp3,"qy_trackp3/F");
	checkFlatteningTreeOutput->Branch("qy_HFm4",&qy_HFm4,"qy_HFm4/F");
	checkFlatteningTreeOutput->Branch("qy_HFp4",&qy_HFp4,"qy_HFp4/F");
	checkFlatteningTreeOutput->Branch("qy_trackmid4",&qy_trackmid4,"qy_trackmid4/F");
	checkFlatteningTreeOutput->Branch("qy_trackm4",&qy_trackm4,"qy_trackm4/F");
	checkFlatteningTreeOutput->Branch("qy_trackp4",&qy_trackp4,"qy_trackp4/F");

	checkFlatteningTreeOutput->Branch("mult_HFm2",&mult_HFm2,"mult_HFm2/F");
	checkFlatteningTreeOutput->Branch("mult_HFp2",&mult_HFp2,"mult_HFp2/F");
	checkFlatteningTreeOutput->Branch("mult_trackmid2",&mult_trackmid2,"mult_trackmid2/F");
	checkFlatteningTreeOutput->Branch("mult_trackm2",&mult_trackm2,"mult_trackm2/F");
	checkFlatteningTreeOutput->Branch("mult_trackp2",&mult_trackp2,"mult_trackp2/F");
	checkFlatteningTreeOutput->Branch("mult_HFm3",&mult_HFm3,"mult_HFm3/F");
	checkFlatteningTreeOutput->Branch("mult_HFp3",&mult_HFp3,"mult_HFp3/F");
	checkFlatteningTreeOutput->Branch("mult_trackmid3",&mult_trackmid3,"mult_trackmid3/F");
	checkFlatteningTreeOutput->Branch("mult_trackm3",&mult_trackm3,"mult_trackm3/F");
	checkFlatteningTreeOutput->Branch("mult_trackp3",&mult_trackp3,"mult_trackp3/F");
	checkFlatteningTreeOutput->Branch("mult_HFm4",&mult_HFm4,"mult_HFm4/F");
	checkFlatteningTreeOutput->Branch("mult_HFp4",&mult_HFp4,"mult_HFp4/F");
	checkFlatteningTreeOutput->Branch("mult_trackmid4",&mult_trackmid4,"mult_trackmid4/F");
	checkFlatteningTreeOutput->Branch("mult_trackm4",&mult_trackm4,"mult_trackm4/F");
	checkFlatteningTreeOutput->Branch("mult_trackp4",&mult_trackp4,"mult_trackp4/F");
  

	// Copy the heavy ion tree to the output
	TTree *heavyIonTreeOutput = new TTree("HiTree","");
	// Connect the branches of the heavy ion tree
	heavyIonTreeOutput->Branch("run",&run,"run/i");
	heavyIonTreeOutput->Branch("evt",&event,"evt/l");
	heavyIonTreeOutput->Branch("lumi",&lumi,"lumi/i");
	heavyIonTreeOutput->Branch("vz",&vertexZ,"vz/F");
	heavyIonTreeOutput->Branch("vx",&vertexX,"vx/F");
	heavyIonTreeOutput->Branch("vy",&vertexY,"vy/F");
	heavyIonTreeOutput->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
	heavyIonTreeOutput->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
	heavyIonTreeOutput->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
	heavyIonTreeOutput->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");	
	heavyIonTreeOutput->Branch("Ntroff",&Ntroff,"Ntroff/I");
	heavyIonTreeOutput->Branch("Nch",&Nch,"Nch/F");
	heavyIonTreeOutput->Branch("NchLoose",&NchLoose,"NchLoose/F");
	heavyIonTreeOutput->Branch("NchTight",&NchTight,"NchTight/F");
	
	// Copy the skim tree to the output
	TTree *skimTreeOutput = new TTree("HltTree","");
	skimTreeOutput->Branch("pPAprimaryVertexFilter",&primaryVertexFilterBit,"pPAprimaryVertexFilter/I");
	skimTreeOutput->Branch("pBeamScrapingFilter",&beamScrapingFilterBit,"pBeamScrapingFilter/I");
	skimTreeOutput->Branch("phfCoincFilter", &hfCoincidenceFilterBit, "phfCoincFilter/I");
	skimTreeOutput->Branch("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, "pVertexFilterCutdz1p0/I");
	skimTreeOutput->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplusBit,"pVertexFilterCutGplus/I");
	skimTreeOutput->Branch("pVertexFilterCutVtx1",&pVertexFilterCutVtx1Bit,"pVertexFilterCutVtx1/I");

	// ========================================== //
	//			Starting matching events	      //
	// ========================================== //

	Int_t ch_events = heavyIonTree->GetEntries(); // number of events
    // loop through jets and create a key for each event
    for(int i_entry = 0; i_entry < ch_events; i_entry++){
       heavyIonTree->GetEntry(i_entry);
       unsigned long long key = keyFromRunLumiEvent(run, lumi, event);
       runLumiEvtToEntryMap[key] = i_entry;
    }

	// ========================================== //
	//				Loop over all events 		  //
	// ========================================== //

	int nEvents = MainV0Tree->GetEntries();
	cout << "There are " << nEvents << " events" << endl;
 	int matchedevents = 0;	
	for(Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
		
		if( iEvent % 100000 == 0 )	std::cout << "iEvent: " << iEvent <<	" of " << nEvents << std::endl;

		// ========================================== //
		//			Start with the V0s	              //
		// ========================================== //
		MainV0Tree->GetEntry(iEvent);
		K0sTree->GetEntry(iEvent);

		//Find matching jet event
		if (V0_evt < 0) continue;
		unsigned long long key = keyFromRunLumiEvent((UInt_t)V0_run,(UInt_t)V0_lumi,(ULong64_t)V0_evt);
		//if (key == 0) cout<<"V0 event "<<V0_evt<<endl;
		//else cout<<"V0 key "<<key<<endl;
        long long i_entry = -1;
        if(runLumiEvtToEntryMap.count(key) == 0) continue; // skip reco event if there is no event match
        else i_entry = runLumiEvtToEntryMap.at(key);
		
		// ========================================== //
		//	Read the event to input trees	      //
		// ========================================== //
		
		heavyIonTree->GetEntry(i_entry);
		skimTree->GetEntry(i_entry);
		trackTree->GetEntry(i_entry);
		checkFlatteningTree->GetEntry(i_entry);

		matchedevents = matchedevents + 1;

		Ntroff = 0;
		Nch = 0.0;
		NchLoose = 0.0;
		NchTight = 0.0;
		
		Nch = get_Ntrkcorr(eff_factor, "nominal", nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);
		NchLoose = get_Ntrkcorr(eff_factorloose, "loose", nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);
		NchTight = get_Ntrkcorr(eff_factortight, "tight", nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);
		
		int multiplicity = get_Ntrkoff(nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);

		bool multsel = true;
		if(ntrkoff==0 || ntrkoff==1){if(iEvent==0){cout << "No multiplicity cut" << endl;}}
		if(ntrkoff==2){if(iEvent==0){cout << "HM 1 to 6: [185,250]" << endl;} if(multiplicity < 185 || multiplicity >= 250){multsel=false;}}
		if(ntrkoff==3){if(iEvent==0){cout << "HM 7: [250,inf]" << endl;} if(multiplicity < 250){multsel=false;}}
		if(multsel==false) continue;		

		Ntroff = multiplicity;

		//Event plane (just what we want EP from 2 to 4)
		epang_HFm2 = (float) eventPlaneAngle[44];
		epang_HFp2 = (float) eventPlaneAngle[45];
		epang_trackmid2 = (float) eventPlaneAngle[46];
		epang_trackm2 = (float) eventPlaneAngle[47];
		epang_trackp2 = (float) eventPlaneAngle[48];
		epang_HFm3 = (float) eventPlaneAngle[73];
		epang_HFp3 = (float) eventPlaneAngle[74];
		epang_trackmid3 = (float) eventPlaneAngle[75];
		epang_trackm3 = (float) eventPlaneAngle[76];
		epang_trackp3 = (float) eventPlaneAngle[77];
		epang_HFm4 = (float) eventPlaneAngle[102];
		epang_HFp4 = (float) eventPlaneAngle[103];
		epang_trackmid4 = (float) eventPlaneAngle[104];
		epang_trackm4 = (float) eventPlaneAngle[105];
		epang_trackp4 = (float) eventPlaneAngle[106];

		q_HFm2 = (float) eventPlaneQ[44];
		q_HFp2 = (float) eventPlaneQ[45];
		q_trackmid2 = (float) eventPlaneQ[46];
		q_trackm2 = (float) eventPlaneQ[47];
		q_trackp2 = (float) eventPlaneQ[48];
		q_HFm3 = (float) eventPlaneQ[73];
		q_HFp3 = (float) eventPlaneQ[74];
		q_trackmid3 = (float) eventPlaneQ[75];
		q_trackm3 = (float) eventPlaneQ[76];
		q_trackp3 = (float) eventPlaneQ[77];
		q_HFm4 = (float) eventPlaneQ[102];
		q_HFp4 = (float) eventPlaneQ[103];
		q_trackmid4 = (float) eventPlaneQ[104];
		q_trackm4 = (float) eventPlaneQ[105];
		q_trackp4 = (float) eventPlaneQ[106];

		qx_HFm2 = (float) eventPlaneQx[44];
		qx_HFp2 = (float) eventPlaneQx[45];
		qx_trackmid2 = (float) eventPlaneQx[46];
		qx_trackm2 = (float) eventPlaneQx[47];
		qx_trackp2 = (float) eventPlaneQx[48];
		qx_HFm3 = (float) eventPlaneQx[73];
		qx_HFp3 = (float) eventPlaneQx[74];
		qx_trackmid3 = (float) eventPlaneQx[75];
		qx_trackm3 = (float) eventPlaneQx[76];
		qx_trackp3 = (float) eventPlaneQx[77];
		qx_HFm4 = (float) eventPlaneQx[102];
		qx_HFp4 = (float) eventPlaneQx[103];
		qx_trackmid4 = (float) eventPlaneQx[104];
		qx_trackm4 = (float) eventPlaneQx[105];
		qx_trackp4 = (float) eventPlaneQx[106];
		
		qy_HFm2 = (float) eventPlaneQy[44];
		qy_HFp2 = (float) eventPlaneQy[45];
		qy_trackmid2 = (float) eventPlaneQy[46];
		qy_trackm2 = (float) eventPlaneQy[47];
		qy_trackp2 = (float) eventPlaneQy[48];
		qy_HFm3 = (float) eventPlaneQy[73];
		qy_HFp3 = (float) eventPlaneQy[74];
		qy_trackmid3 = (float) eventPlaneQy[75];
		qy_trackm3 = (float) eventPlaneQy[76];
		qy_trackp3 = (float) eventPlaneQy[77];
		qy_HFm4 = (float) eventPlaneQy[102];
		qy_HFp4 = (float) eventPlaneQy[103];
		qy_trackmid4 = (float) eventPlaneQy[104];
		qy_trackm4 = (float) eventPlaneQy[105];
		qy_trackp4 = (float) eventPlaneQy[106];

		mult_HFm2 = (float) eventPlaneMultiplicity[44];
		mult_HFp2 = (float) eventPlaneMultiplicity[45];
		mult_trackmid2 = (float) eventPlaneMultiplicity[46];
		mult_trackm2 = (float) eventPlaneMultiplicity[47];
		mult_trackp2 = (float) eventPlaneMultiplicity[48];
		mult_HFm3 = (float) eventPlaneMultiplicity[73];
		mult_HFp3 = (float) eventPlaneMultiplicity[74];
		mult_trackmid3 = (float) eventPlaneMultiplicity[75];
		mult_trackm3 = (float) eventPlaneMultiplicity[76];
		mult_trackp3 = (float) eventPlaneMultiplicity[77];
		mult_HFm4 = (float) eventPlaneMultiplicity[102];
		mult_HFp4 = (float) eventPlaneMultiplicity[103];
		mult_trackmid4 = (float) eventPlaneMultiplicity[104];
		mult_trackm4 = (float) eventPlaneMultiplicity[105];
		mult_trackp4 = (float) eventPlaneMultiplicity[106];

	    // --> Loop over 2 K0s
    	if( K0s_pt->size() > 1 ){
    	
    		int totalsize = K0s_pt->size() * K0s_pt->size();
  	 		for(int idx = 0; idx < totalsize; idx++){ // 2 loops in one
   		
    			int ik0s1 = idx / K0s_pt->size();
    			int ik0s2 = idx % K0s_pt->size(); 
    			if( ik0s2 <= ik0s1) continue; // j = i + 1			

				// start doing K0s calculations
				// K0s 1 	
				if(K0s_pt->at(ik0s1) < ptmin) continue;   //Minimum pT   		
				if(TMath::Abs(K0s_eta->at(ik0s1)) > etamin) continue; //eta acceptance
				if(K0s_d1Nhit->at(ik0s1) < 4)continue;
				if(K0s_d2Nhit->at(ik0s1) < 4)continue;			
				if(K0s_dca->at(ik0s1) > 0.5)continue;
				if(fabs(K0s_dxy1->at(ik0s1)) < 1.0)continue;
				if(fabs(K0s_dxy2->at(ik0s1)) < 1.0)continue;
				if(fabs(K0s_dz1->at(ik0s1)) < 1.0)continue;
				if(fabs(K0s_dz2->at(ik0s1)) < 1.0)continue;
				if(K0s_3Dagl->at(ik0s1) <= 0.997)continue;
				if(K0s_3Ddl->at(ik0s1) <= 3.5)continue;
				
				float Pxp_k0s1 = (float) K0s_d1px->at(ik0s1);
				float Pyp_k0s1 = (float) K0s_d1py->at(ik0s1);
				float Pzp_k0s1 = (float) K0s_d1pz->at(ik0s1);
				float Pmp_k0s1 = (float) K0s_d1M->at(ik0s1);

				float Pxn_k0s1 = (float) K0s_d2px->at(ik0s1);
				float Pyn_k0s1 = (float) K0s_d2py->at(ik0s1);
				float Pzn_k0s1 = (float) K0s_d2pz->at(ik0s1);
				float Pmn_k0s1 = (float) K0s_d2M->at(ik0s1);

				TLorentzVector pvector1_k0s1;
				pvector1_k0s1.SetXYZM(Pxp_k0s1,Pyp_k0s1,Pzp_k0s1,Pmp_k0s1);
				TLorentzVector pvector2_k0s1;
				pvector2_k0s1.SetXYZM(Pxn_k0s1,Pyn_k0s1,Pzn_k0s1,Pmn_k0s1);

				TVector3 dauvec1_k0s1(K0s_d1px->at(ik0s1), K0s_d1py->at(ik0s1), K0s_d1pz->at(ik0s1));
				TVector3 dauvec2_k0s1(K0s_d2px->at(ik0s1), K0s_d2py->at(ik0s1), K0s_d2pz->at(ik0s1));
				TVector3 dauvecsum_k0s1(dauvec1_k0s1+dauvec2_k0s1);

				// Armenteros-Podolanski
 				float Pp_k0s1=sqrt(Pxp_k0s1*Pxp_k0s1+Pyp_k0s1*Pyp_k0s1+Pzp_k0s1*Pzp_k0s1);
				float Pn_k0s1=sqrt(Pxn_k0s1*Pxn_k0s1+Pyn_k0s1*Pyn_k0s1+Pzn_k0s1*Pzn_k0s1);
				float PN_k0s1=sqrt((Pxp_k0s1+Pxn_k0s1)*(Pxp_k0s1+Pxn_k0s1)+(Pyp_k0s1+Pyn_k0s1)*(Pyp_k0s1+Pyn_k0s1)+(Pzp_k0s1+Pzn_k0s1)*(Pzp_k0s1+Pzn_k0s1));
 				float cosx1_k0s1=(Pxp_k0s1*(Pxp_k0s1+Pxn_k0s1)+Pyp_k0s1*(Pyp_k0s1+Pyn_k0s1)+Pzp_k0s1*(Pzp_k0s1+Pzn_k0s1))/(Pp_k0s1*PN_k0s1);
				float cosx2_k0s1=(Pxn_k0s1*(Pxp_k0s1+Pxn_k0s1)+Pyn_k0s1*(Pyp_k0s1+Pyn_k0s1)+Pzn_k0s1*(Pzp_k0s1+Pzn_k0s1))/(Pn_k0s1*PN_k0s1);
				float QT_k0s1 = Pp_k0s1*sqrt(1.0-cosx1_k0s1*cosx1_k0s1);
				float Alpha_k0s1 = (Pp_k0s1*cosx1_k0s1-Pn_k0s1*cosx2_k0s1)/(Pp_k0s1*cosx1_k0s1+Pn_k0s1*cosx2_k0s1);
 
				float energyd1e_k0s1 = sqrt(electronMass*electronMass+pvector1_k0s1.P()*pvector1_k0s1.P());
				float energyd2e_k0s1 = sqrt(electronMass*electronMass+pvector2_k0s1.P()*pvector2_k0s1.P());
				float invmass_ee_k0s1 = sqrt((energyd1e_k0s1+energyd2e_k0s1)*(energyd1e_k0s1+energyd2e_k0s1)-dauvecsum_k0s1.Mag2());

				float massd1_k0s1, massd2_k0s1, energyd1_k0s1, energyd2_k0s1, invmasspiproton_k0s1, invmassprotonpi_k0s1;

				massd1_k0s1 = pimass;
				massd2_k0s1 = pro_mass;
				energyd1_k0s1 = sqrt(massd1_k0s1*massd1_k0s1+pvector1_k0s1.P()*pvector1_k0s1.P());
				energyd2_k0s1 = sqrt(massd2_k0s1*massd2_k0s1+pvector2_k0s1.P()*pvector2_k0s1.P());
				invmasspiproton_k0s1 = sqrt((energyd1_k0s1+energyd2_k0s1)*(energyd1_k0s1+energyd2_k0s1)-dauvecsum_k0s1.Mag2());

				massd1_k0s1 = pro_mass;
				massd2_k0s1 = pimass;
				energyd1_k0s1 = sqrt(massd1_k0s1*massd1_k0s1+pvector1_k0s1.P()*pvector1_k0s1.P());
				energyd2_k0s1 = sqrt(massd2_k0s1*massd2_k0s1+pvector2_k0s1.P()*pvector2_k0s1.P());
				invmassprotonpi_k0s1 = sqrt((energyd1_k0s1+energyd2_k0s1)*(energyd1_k0s1+energyd2_k0s1)-dauvecsum_k0s1.Mag2());
				
				// K0s 2
				if(K0s_pt->at(ik0s2) < ptmin) continue;   //Minimum pT   		
				if(TMath::Abs(K0s_eta->at(ik0s2)) > etamin) continue; //eta acceptance
				if(K0s_d1Nhit->at(ik0s2) < 4)continue;
				if(K0s_d2Nhit->at(ik0s2) < 4)continue;			
				if(K0s_dca->at(ik0s2) > 0.5)continue;
				if(fabs(K0s_dxy1->at(ik0s2)) < 1.0)continue;
				if(fabs(K0s_dxy2->at(ik0s2)) < 1.0)continue;
				if(fabs(K0s_dz1->at(ik0s2)) < 1.0)continue;
				if(fabs(K0s_dz2->at(ik0s2)) < 1.0)continue;
				if(K0s_3Dagl->at(ik0s2) <= 0.997)continue;
				if(K0s_3Ddl->at(ik0s2) <= 3.5)continue;
				
				float Pxp_k0s2 = (float) K0s_d1px->at(ik0s2);
				float Pyp_k0s2 = (float) K0s_d1py->at(ik0s2);
				float Pzp_k0s2 = (float) K0s_d1pz->at(ik0s2);
				float Pmp_k0s2 = (float) K0s_d1M->at(ik0s2);

				float Pxn_k0s2 = (float) K0s_d2px->at(ik0s2);
				float Pyn_k0s2 = (float) K0s_d2py->at(ik0s2);
				float Pzn_k0s2 = (float) K0s_d2pz->at(ik0s2);
				float Pmn_k0s2 = (float) K0s_d2M->at(ik0s2);

				TLorentzVector pvector1_k0s2;
				pvector1_k0s2.SetXYZM(Pxp_k0s2,Pyp_k0s2,Pzp_k0s2,Pmp_k0s2);
				TLorentzVector pvector2_k0s2;
				pvector2_k0s2.SetXYZM(Pxn_k0s2,Pyn_k0s2,Pzn_k0s2,Pmn_k0s2);

				TVector3 dauvec1_k0s2(K0s_d1px->at(ik0s2), K0s_d1py->at(ik0s2), K0s_d1pz->at(ik0s2));
				TVector3 dauvec2_k0s2(K0s_d2px->at(ik0s2), K0s_d2py->at(ik0s2), K0s_d2pz->at(ik0s2));
				TVector3 dauvecsum_k0s2(dauvec1_k0s2+dauvec2_k0s2);

				// Armenteros-Podolanski
 				float Pp_k0s2=sqrt(Pxp_k0s2*Pxp_k0s2+Pyp_k0s2*Pyp_k0s2+Pzp_k0s2*Pzp_k0s2);
				float Pn_k0s2=sqrt(Pxn_k0s2*Pxn_k0s2+Pyn_k0s2*Pyn_k0s2+Pzn_k0s2*Pzn_k0s2);
				float PN_k0s2=sqrt((Pxp_k0s2+Pxn_k0s2)*(Pxp_k0s2+Pxn_k0s2)+(Pyp_k0s2+Pyn_k0s2)*(Pyp_k0s2+Pyn_k0s2)+(Pzp_k0s2+Pzn_k0s2)*(Pzp_k0s2+Pzn_k0s2));
 				float cosx1_k0s2=(Pxp_k0s2*(Pxp_k0s2+Pxn_k0s2)+Pyp_k0s2*(Pyp_k0s2+Pyn_k0s2)+Pzp_k0s2*(Pzp_k0s2+Pzn_k0s2))/(Pp_k0s2*PN_k0s2);
				float cosx2_k0s2=(Pxn_k0s2*(Pxp_k0s2+Pxn_k0s2)+Pyn_k0s2*(Pyp_k0s2+Pyn_k0s2)+Pzn_k0s2*(Pzp_k0s2+Pzn_k0s2))/(Pn_k0s2*PN_k0s2);
				float QT_k0s2 = Pp_k0s2*sqrt(1.0-cosx1_k0s2*cosx1_k0s2);
				float Alpha_k0s2 = (Pp_k0s2*cosx1_k0s2-Pn_k0s2*cosx2_k0s2)/(Pp_k0s2*cosx1_k0s2+Pn_k0s2*cosx2_k0s2);
 
				float energyd1e_k0s2 = sqrt(electronMass*electronMass+pvector1_k0s2.P()*pvector1_k0s2.P());
				float energyd2e_k0s2 = sqrt(electronMass*electronMass+pvector2_k0s2.P()*pvector2_k0s2.P());
				float invmass_ee_k0s2 = sqrt((energyd1e_k0s2+energyd2e_k0s2)*(energyd1e_k0s2+energyd2e_k0s2)-dauvecsum_k0s2.Mag2());

				float massd1_k0s2, massd2_k0s2, energyd1_k0s2, energyd2_k0s2, invmasspiproton_k0s2, invmassprotonpi_k0s2;

				massd1_k0s2 = pimass;
				massd2_k0s2 = pro_mass;
				energyd1_k0s2 = sqrt(massd1_k0s2*massd1_k0s2+pvector1_k0s2.P()*pvector1_k0s2.P());
				energyd2_k0s2 = sqrt(massd2_k0s2*massd2_k0s2+pvector2_k0s2.P()*pvector2_k0s2.P());
				invmasspiproton_k0s2 = sqrt((energyd1_k0s2+energyd2_k0s2)*(energyd1_k0s2+energyd2_k0s2)-dauvecsum_k0s2.Mag2());

				massd1_k0s2 = pro_mass;
				massd2_k0s2 = pimass;
				energyd1_k0s2 = sqrt(massd1_k0s2*massd1_k0s2+pvector1_k0s2.P()*pvector1_k0s2.P());
				energyd2_k0s2 = sqrt(massd2_k0s2*massd2_k0s2+pvector2_k0s2.P()*pvector2_k0s2.P());
				invmassprotonpi_k0s2 = sqrt((energyd1_k0s2+energyd2_k0s2)*(energyd1_k0s2+energyd2_k0s2)-dauvecsum_k0s2.Mag2());
	
				// Construct Lorentz vectors
    			LorentzVector neutralkaon1(K0s_pt->at(ik0s1), K0s_eta->at(ik0s1), K0s_phi->at(ik0s1), K0s_mass->at(ik0s1));
    			LorentzVector neutralkaon2(K0s_pt->at(ik0s2), K0s_eta->at(ik0s2), K0s_phi->at(ik0s2), K0s_mass->at(ik0s2));

			    // Combine them and compute invariant mass
			    LorentzVector system = neutralkaon1 + neutralkaon2;
				if( system.M() <= 1.2 ) continue;
    			if( system.M() >= 1.8 ) continue;
    			if( system.Pt() <= 0.5 ) continue;
    			if( fabs(system.Eta()) > 2.4 ) continue;    			
				
				// K0s1 part
				K0s1_LooseVector->push_back(1);
				if( K0s_3Dagl->at(ik0s1) > 0.999 && K0s_3Ddl->at(ik0s1) > 5.0 ){ K0s1_NominalVector->push_back(1); }else{ K0s1_NominalVector->push_back(0); }
				if( K0s_3Dagl->at(ik0s1) > 0.9995 && K0s_3Ddl->at(ik0s1) > 6.5 && fabs(K0s_dxy1->at(ik0s1)) >= 1.2 && fabs(K0s_dxy2->at(ik0s1)) >= 1.2 && fabs(K0s_dz1->at(ik0s1)) >= 1.2 && fabs(K0s_dz2->at(ik0s1)) >= 1.2 ){ K0s1_TightVector->push_back(1); }else{ K0s1_TightVector->push_back(0); }
				K0s1_MisIDEEVector->push_back(invmass_ee_k0s1);
				K0s1_MisIDPiProtonVector->push_back(invmasspiproton_k0s1);
				K0s1_MisIDProtonPiVector->push_back(invmassprotonpi_k0s1);
				K0s1_APAlphaVector->push_back(Alpha_k0s1);
				K0s1_APQTVector->push_back(QT_k0s1);
				K0s1_PtVector->push_back(K0s_pt->at(ik0s1));
				K0s1_EtaVector->push_back(K0s_eta->at(ik0s1));
				K0s1_PhiVector->push_back(K0s_phi->at(ik0s1));
				K0s1_MassVector->push_back(K0s_mass->at(ik0s1));
				K0s1_D1_PtVector->push_back(dauvec1_k0s1.Pt());
				K0s1_D1_EtaVector->push_back(dauvec1_k0s1.Eta());
				K0s1_D1_PhiVector->push_back(dauvec1_k0s1.Phi());
				K0s1_D2_PtVector->push_back(dauvec2_k0s1.Pt());
				K0s1_D2_EtaVector->push_back(dauvec2_k0s1.Eta());
				K0s1_D2_PhiVector->push_back(dauvec2_k0s1.Phi());

				// K0s2 part
				K0s2_LooseVector->push_back(1);
				if( K0s_3Dagl->at(ik0s2) > 0.999 && K0s_3Ddl->at(ik0s2) > 5.0 ){ K0s2_NominalVector->push_back(1); }else{ K0s2_NominalVector->push_back(0); }
				if( K0s_3Dagl->at(ik0s2) > 0.9995 && K0s_3Ddl->at(ik0s2) > 6.5 && fabs(K0s_dxy1->at(ik0s2)) >= 1.2 && fabs(K0s_dxy2->at(ik0s2)) >= 1.2 && fabs(K0s_dz1->at(ik0s2)) >= 1.2 && fabs(K0s_dz2->at(ik0s2)) >= 1.2 ){ K0s2_TightVector->push_back(1); }else{ K0s2_TightVector->push_back(0); }
				K0s2_MisIDEEVector->push_back(invmass_ee_k0s2);
				K0s2_MisIDPiProtonVector->push_back(invmasspiproton_k0s2);
				K0s2_MisIDProtonPiVector->push_back(invmassprotonpi_k0s2);
				K0s2_APAlphaVector->push_back(Alpha_k0s2);
				K0s2_APQTVector->push_back(QT_k0s2);
				K0s2_PtVector->push_back(K0s_pt->at(ik0s2));
				K0s2_EtaVector->push_back(K0s_eta->at(ik0s2));
				K0s2_PhiVector->push_back(K0s_phi->at(ik0s2));
				K0s2_MassVector->push_back(K0s_mass->at(ik0s2));
				K0s2_D1_PtVector->push_back(dauvec1_k0s2.Pt());
				K0s2_D1_EtaVector->push_back(dauvec1_k0s2.Eta());
				K0s2_D1_PhiVector->push_back(dauvec1_k0s2.Phi());
				K0s2_D2_PtVector->push_back(dauvec2_k0s2.Pt());
				K0s2_D2_EtaVector->push_back(dauvec2_k0s2.Eta());
				K0s2_D2_PhiVector->push_back(dauvec2_k0s2.Phi());
								
				if( (K0s_chi21->at(ik0s1) == K0s_chi21->at(ik0s2)) || (K0s_chi22->at(ik0s1) == K0s_chi22->at(ik0s2)) || (K0s_chi21->at(ik0s1) == K0s_chi22->at(ik0s2)) || (K0s_chi22->at(ik0s1) == K0s_chi21->at(ik0s2)) ) { K0s1_ShareDauVector->push_back(1); K0s2_ShareDauVector->push_back(1); } else { K0s1_ShareDauVector->push_back(0); K0s2_ShareDauVector->push_back(0); }

				// f2 prime part
				f2prime_PtVector->push_back(system.Pt());
				f2prime_EtaVector->push_back(system.Eta());
				f2prime_PhiVector->push_back(system.Phi());
				f2prime_MassVector->push_back(system.M());
				
      		}
      	
      	}
		if(f2prime_PtVector->size() > 0){ // only fill tree if event has a F2Prime candidate
			heavyIonTreeOutput->Fill(); // fill event information
			skimTreeOutput->Fill();		// filter information
			checkFlatteningTreeOutput->Fill(); // fill EP information	
     		F2PrimeTreeOutput->Fill(); // f2 prime information
		}
		
		// Clear the vectors before the next event! Otherwise all the K0s pile up cumulatively		
		K0s1_ShareDauVector->clear();
		K0s1_LooseVector->clear();
		K0s1_NominalVector->clear();
		K0s1_TightVector->clear();
		K0s1_MisIDEEVector->clear();
		K0s1_MisIDPiProtonVector->clear();
		K0s1_MisIDProtonPiVector->clear();
		K0s1_APAlphaVector->clear();
		K0s1_APQTVector->clear();
		K0s1_PtVector->clear();
		K0s1_EtaVector->clear();
		K0s1_PhiVector->clear();
		K0s1_MassVector->clear();

		K0s2_ShareDauVector->clear();
		K0s2_LooseVector->clear();
		K0s2_NominalVector->clear();
		K0s2_TightVector->clear();
		K0s2_MisIDEEVector->clear();
		K0s2_MisIDPiProtonVector->clear();
		K0s2_MisIDProtonPiVector->clear();
		K0s2_APAlphaVector->clear();
		K0s2_APQTVector->clear();
		K0s2_PtVector->clear();
		K0s2_EtaVector->clear();
		K0s2_PhiVector->clear();
		K0s2_MassVector->clear();
		
		f2prime_PtVector->clear();
		f2prime_EtaVector->clear();
		f2prime_PhiVector->clear();
		f2prime_MassVector->clear();

	} // End loop over events

	cout << "Number of matched events: " << matchedevents << endl;

	// Write the skimmed trees to the output file

	gDirectory->mkdir("hiEvtAnalyzer");
	gDirectory->cd("hiEvtAnalyzer");
	heavyIonTreeOutput->Write();
  
	gDirectory->cd("../");
	gDirectory->mkdir("skimanalysis");
	gDirectory->cd("skimanalysis");
	skimTreeOutput->Write();
	
	gDirectory->cd("../");
	gDirectory->mkdir("checkflattening");
	gDirectory->cd("checkflattening");
	checkFlatteningTreeOutput->Write();

	gDirectory->cd("../");
	gDirectory->mkdir("F2Prime");
	gDirectory->cd("F2Prime");
	F2PrimeTreeOutput->Write();	
	gDirectory->cd("../");
  
	outputFile->Close();

	cout << endl;
	cout << "------------------------------------- SKIMMER DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

	print_stop(); // Print time, date and hour when it stops

}

unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event){

  const unsigned long long runMult = 1;
  const unsigned long long lumiMult = 1000000;
  const unsigned long long evtMult = 10000000000;
  const unsigned long long evtLimit = 10000000000;

  unsigned long long key = 0;
  if(event >= evtLimit){
    std::cout << "RUNLUMIEVENTKEY WARNING : \'" << event << "\' is greated that event limit 10^10. returning key 0" << std::endl;
    return key;
  }

  key += runMult* static_cast<unsigned long long>(run);
  key += lumiMult* static_cast<unsigned long long>(lumi);
  key += evtMult*event;

  //std::cout << "key = " << key << std::endl;
  
  return key;
  
}

int main(int argc, char** argv){
				TString firstArgument(argv[1]);
				TString secArgument(argv[2]);				
				TString outfile(argv[3]);
				int ntrkoffline = atoi(argv[4]);
				f2prime(firstArgument, secArgument,outfile, ntrkoffline);
}
