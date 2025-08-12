int drawf(int runnumdraw = -1, int evtdraw = -1)
{
  TFile* evtfile = TFile::Open("../events/allwf.root","READ");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  int runnum, evtnum;
  unsigned int emwf[96][256][12];
  unsigned int ihwf[24][64][12];
  unsigned int ohwf[24][64][12];

  TTree* wft = (TTree*) evtfile->Get("wft");

  wft->SetBranchAddress("runnum",&runnum);
  wft->SetBranchAddress("evtnum",&evtnum);
  wft->SetBranchAddress("emwf",emwf);
  wft->SetBranchAddress("ihwf",ihwf);
  wft->SetBranchAddress("ohwf",ohwf);

  string calo[3] = {"EMCal","IHCal","OHCal"};
  
  TH2D* h2_wf[3];
  for(int i=0; i<3; ++i)
    {
      h2_wf[i] = new TH2D(("h2_wf"+to_string(i)).c_str(),(";Sample;"+calo[i]+" ADC Value;Counts").c_str(),12,-0.5,11.5,17000,0,17000);
    }

  TCanvas* c = new TCanvas("","",1000,1000);
  c->SetTopMargin(0.15);
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);

  string texts[3] = {"#bf{#it{sPHENIX}} Internal","#sqrt{s} = 200 GeV",""};
  TLatex* tex[3];
  float ycoord[3] = {0.96,0.91,0.86};
  for(int i=0; i<3; ++i)
    {
      tex[i] = new TLatex(0.7,ycoord[i],texts[i].c_str());
      tex[i]->SetTextFont(42);
      tex[i]->SetTextSize(0.04);
      tex[i]->SetTextColor(kBlack);
      tex[i]->SetLineWidth(1);
      tex[i]->SetNDC();
    }
  
  for(int i=0; i<wft->GetEntries(); ++i)
    {
      wft->GetEntry(i);
      if(runnumdraw >= 0)
	{
	  if(runnum != runnumdraw) continue;
	  if(evtdraw >= 0)
	    {
	      if(evtnum != evtdraw) continue;
	    }
	}
      texts[2] = "Run " + to_string(evtnum) + ", Event " + to_string(evtnum);
      for(int j=0; j<3; ++j)
	{
	  if(j==0)
	    {
	      for(int k=0; k<96; ++k)
		{
		  for(int l=0; l<256; ++l)
		    {
		      for(int m=0; m<12; ++m)
			{
			  h2_wf[j]->Fill(m,emwf[k][l][m]);
			}
		    }
		}
	    }
	  else if(j==1)
	    {
	      for(int k=0; k<24; ++k)
		{
		  for(int l=0; l<64; ++l)
		    {
		      for(int m=0; m<12; ++m)
			{
			  h2_wf[j]->Fill(m,ihwf[k][l][m]);
			}
		    }
		}
	    }
	  else if(j==2)
	    {
	      for(int k=0; k<24; ++k)
		{
		  for(int l=0; l<64; ++l)
		    {
		      for(int m=0; m<12; ++m)
			{
			  h2_wf[j]->Fill(m,ohwf[k][l][m]);
			}
		    }
		}
	    }
	}

      for(int j=0; j<3; ++j)
	{
	  h2_wf[j]->Draw("COLZ");
	  for(int k=0; k<3; ++k)
	    {
	      tex[k]->Draw();
	    }
	  gPad->SaveAs(("../images/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+calo[j]+"_wf.png").c_str());
	  h2_wf[j]->Reset();
	}
    }
  return 0;
}
