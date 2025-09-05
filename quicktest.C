int quicktest()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* file = TFile::Open("hadded_waveform_20250902.root");
  TTree* tree = (TTree*)file->Get("wft");
  TFile* afile = TFile::Open("hadded_chi2file_20250902.root");
  TTree* atree = (TTree*)afile->Get("jet_tree");

  float emadcfit[96][256];
  float ihadcfit[24][64];
  float ohadcfit[24][64];

  unsigned int emwf[96][256][12];
  unsigned int ihwf[24][64][12];
  unsigned int ohwf[24][64][12];
  
  float emtow[96][256];
  float ihtow[24][64];
  float ohtow[24][64];

  int evtnum, aevtnum, runnum, arunnum, failscut;
  
  tree->SetBranchAddress("emadcfit",emadcfit);
  tree->SetBranchAddress("ihadcfit",ihadcfit);
  tree->SetBranchAddress("ohadcfit",ohadcfit);
  tree->SetBranchAddress("emwf",emwf);
  tree->SetBranchAddress("ihwf",ihwf);
  tree->SetBranchAddress("ohwf",ohwf);
  atree->SetBranchAddress("emtow",emtow);
  atree->SetBranchAddress("ihtow",ihtow);
  atree->SetBranchAddress("ohtow",ohtow);
  tree->SetBranchAddress("evtnum",&evtnum);
  atree->SetBranchAddress("evtnum",&aevtnum);
  tree->SetBranchAddress("runnum",&runnum);
  atree->SetBranchAddress("runnum",&arunnum);
  atree->SetBranchAddress("failscut",&failscut);
  
  
  TH3D* h3_fit_wfd_cte[3];

  string texts[3] = {"#bf{#it{sPHENIX}} Internal","#sqrt{s} = 200 GeV","Events pass dijet cut"};
  TLatex* tex[3];
  float ycoord[3] = {0.96,0.91,0.86};
  for(int i=0; i<3; ++i)
    {
      tex[i] = new TLatex(i==2?0.5:0.7,ycoord[i],texts[i].c_str());
      tex[i]->SetTextFont(42);
      tex[i]->SetTextSize(0.04);
      tex[i]->SetTextColor(kBlack);
      tex[i]->SetLineWidth(1);
      tex[i]->SetNDC();
    }

  string calo[3] = {"EMCal","IHCal","OHCal"};

  for(int i=0; i<3; ++i)
    {
      h3_fit_wfd_cte[i] = new TH3D(("h3_fit_wfd"+to_string(i)).c_str(),(";"+calo[i]+" Fitted ADC Value;"+calo[i]+" Max Difference Between Samples;"+calo[i]+" Calibrated Tower E [GeV]").c_str(),200,0,20000,170,0,17000,120,0,(i==0?60:(i==1?20:120)));
      h3_fit_wfd_cte[i]->GetXaxis()->SetNdivisions(5);
      h3_fit_wfd_cte[i]->GetYaxis()->SetNdivisions(5);
      h3_fit_wfd_cte[i]->GetZaxis()->SetTitleOffset(1.5);
    }

 TCanvas* c = new TCanvas("","",1000,1000);
  c->SetLogz();
  c->SetTopMargin(0.2);
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  int jte = 0;
  for(int i=0; i<tree->GetEntries(); ++i)
    {
      tree->GetEntry(i);
      //if(evtnum != 377163) continue;
      if(i%100 == 0) cout << i << endl;
      atree->GetEntry(jte);
      int flag = 0;

      while(arunnum != runnum || evtnum != aevtnum)
	{
	  ++jte;
	  if(jte == atree->GetEntries())
	    {
	      flag = 1;
	      jte = 0;
	      break;
	    }
	  atree->GetEntry(jte);
	}
      if(flag) continue;

      if(failscut != 0 && failscut != 2) continue;
      
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<(j==0?96:24); ++k)
	    {
	      for(int l=0; l<(j==0?256:64); ++l)
		{
		  unsigned int minval = 99999;
		  unsigned int maxval = 0;
		  
		  for(int m=0; m<12; ++m)
		    {
		      unsigned int value = (j==0?emwf[k][l][m]:(j==1?ihwf[k][l][m]:ohwf[k][l][m]));
		      if(value > maxval) maxval = value;
		      if(value < minval) minval = value;
		    }
		  //cout << emadcfit[k][l] << " " << maxval - minval << endl;
		  if(maxval - minval == 0) continue;
		  if(j==0)
		    {
		      h3_fit_wfd_cte[j]->Fill(emadcfit[k][l],maxval-minval,emtow[k][l]);
		    }
		  if(j==1)
		    {
		      h3_fit_wfd_cte[j]->Fill(ihadcfit[k][l],maxval-minval,ihtow[k][l]);
		    }
		  if(j==2)
		    {
		      h3_fit_wfd_cte[j]->Fill(ohadcfit[k][l],maxval-minval,ohtow[k][l]);
		    }
		}
	    }
	}
    }

  for(int i=0; i<3; ++i)
    {
      h3_fit_wfd_cte[i]->Project3D("yx")->Draw("COLZ");
      TLine* line = new TLine(0,0,20000,20000);
      line->SetLineColor(kRed);
      line->Draw();
      for(int j=0; j<3; ++j) tex[j]->Draw();
      c->SaveAs((calo[i]+"fitwfd.png").c_str());
      h3_fit_wfd_cte[i]->Project3D("zx")->Draw("COLZ");
      for(int j=0; j<3; ++j) tex[j]->Draw();
      c->SaveAs((calo[i]+"fitcte.png").c_str());
      h3_fit_wfd_cte[i]->Project3D("zy")->Draw("COLZ");
      for(int j=0; j<3; ++j) tex[j]->Draw();
      c->SaveAs((calo[i]+"ctewfd.png").c_str());
      delete line;
    }
  return 0;
}
