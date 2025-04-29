bool checkBounds(double pt, double eta){
  if( TMath::Abs(eta) > 2.4 ){return false;}
  if( pt < 0 || pt > 500 ){return false;}
  return true;
}

// pT and eta are the transverse momentum and pseudorapidity of the track (considering a 2D histogram where X is eta axis and Y pT axis)
double getTrkCorrWeight(TH2 *eff_factor, double pT, double eta){

  if( !checkBounds(pT, eta) ) return 0;
  double factor = 1.0;

  // add it for each system here
  double eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(eta), eff_factor->GetYaxis()->FindBin(pT) );
  if(eff > 0.9999 || eff < 0.0001) eff = 1.0;
  factor = (1.0 / eff); //only efficiency

  return factor;

}


int get_Ntrkoff(int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr){
	int Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){ 
        if(pt[ii] <= 0.4) continue;
		if(fabs(eta[ii]) >= 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] != 1) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}


float get_Ntrkcorr(TH2 *trkeff_file, TString type, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr){
	float Ntrk_off = 0.0;
	float dxyzcut = 3.0;
	float ptrescut = 0.1;
	if(type == "nominal"){dxyzcut = 3.0; ptrescut = 0.1;}
	if(type == "tight"){dxyzcut = 2.0; ptrescut = 0.05;}
	if(type == "loose"){dxyzcut = 5.0; ptrescut = 0.1;}
	for(int ii=0; ii<size; ii++){ 
		if(pt[ii] <= 0.4) continue;
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] != 1) continue;
		if(fabs(pterr[ii]/pt[ii]) >= ptrescut) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= dxyzcut) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= dxyzcut) continue;
		Ntrk_off = Ntrk_off + getTrkCorrWeight(trkeff_file, pt[ii], eta[ii]);
	}
	return Ntrk_off;
}
