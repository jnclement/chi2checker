int max(int a, int b)
{
  return (a>b?a:b);
}

int min(int a, int b)
{
  return (a>b?b:a);
}

int drawf(int lo, int hi, int runnumdraw = -1, int evtdraw = -1)
{
  TFile* evtfile = TFile::Open("../events/allwf.root","READ");
  TFile* othfile = TFile::Open("../events/allevents_20250815.root","READ");
  
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  int runnum, evtnum, failscut, rnj, enj;
  unsigned int emwf[96][256][12];
  unsigned int ihwf[24][64][12];
  unsigned int ohwf[24][64][12];
  float jconem[24][64];
  float jconih[24][64];
  float jconoh[24][64];
  
  TTree* wft = (TTree*) evtfile->Get("wft");
  TTree* jt = (TTree*) othfile->Get("jet_tree");

  const float jetcut = 10;

  jt->SetBranchAddress("jconem",jconem);
  jt->SetBranchAddress("jconih",jconih);
  jt->SetBranchAddress("jconoh",jconoh);
  jt->SetBranchAddress("runnum",&rnj);
  jt->SetBranchAddress("evtnum",&enj);
  
  wft->SetBranchAddress("runnum",&runnum);
  wft->SetBranchAddress("evtnum",&evtnum);
  wft->SetBranchAddress("emwf",emwf);
  wft->SetBranchAddress("ihwf",ihwf);
  wft->SetBranchAddress("ohwf",ohwf);
  wft->SetBranchAddress("failscut",&failscut);
  
  string calo[3] = {"EMCal","IHCal","OHCal"};
  string ct[4] = {"fails","dijet","frac","both"};

  TH1D* singlewf = new TH1D("singlwf",";Sample;ADC Value",12,-0.5,11.5);
  
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

  string texts[5] = {"#bf{#it{sPHENIX}} Internal","#sqrt{s} = 200 GeV","","Waveforms from jets E_{T}>10 GeV",""};
  TLatex* tex[5];
  float ycoord[5] = {0.96,0.91,0.86,0.96,0.91};
  for(int i=0; i<4; ++i)
    {
      tex[i] = new TLatex((i>2)?0.1:0.7,ycoord[i],texts[i].c_str());
      tex[i]->SetTextFont(42);
      tex[i]->SetTextSize(0.04);
      tex[i]->SetTextColor(kBlack);
      tex[i]->SetLineWidth(1);
      tex[i]->SetNDC();
    }
  delete tex[2];
  int jte = lo;
  
  for(int i=lo; i<(hi>wft->GetEntries()?wft->GetEntries():hi); ++i)
    {
      wft->GetEntry(i);
      jt->GetEntry(jte);
      int flag = 0;
      while(rnj != runnum || enj != evtnum)
	{
	  ++jte;
	  if(jte == jt->GetEntries())
	    {
	      flag = 1;
	      break;
	    }
	  jt->GetEntry(jte);
	}
      if(flag) continue;
      if((failscut < 0 || failscut > 2) && i%100 != 0) continue;
      if(runnumdraw >= 0)
	{
	  if(runnum != runnumdraw) continue;
	  if(evtdraw >= 0)
	    {
	      if(evtnum != evtdraw) continue;
	    }
	}
      texts[2] = "Run " + to_string(runnum) + ", Event " + to_string(evtnum);
      tex[2] = new TLatex(0.5,ycoord[2],texts[2].c_str());
      tex[2]->SetTextFont(42);
      tex[2]->SetTextSize(0.04);
      tex[2]->SetTextColor(kBlack);
      tex[2]->SetLineWidth(1);
      tex[2]->SetNDC();
      for(int j=0; j<3; ++j)
	{
	  if(j==0)
	    {
	      for(int k=0; k<96; ++k)
		{
		  for(int l=0; l<256; ++l)
		    {
		      //cout << k << " " << l << " " << jconem[k/4][l/4] << endl;
		      //cout << k/4 << " " << l/4 << endl;
		      //continue;
		      if(jconem[k/4][l/4] < jetcut) continue;
		      int minval = 99999999;
		      int maxval = -9999999;
		      for(int m=0; m<12; ++m)
			{
			  if (emwf[k][l][m] > 0) h2_wf[j]->Fill(m,emwf[k][l][m]);
			  if(runnumdraw > -1 && evtdraw > -1) singlewf->Fill(m,emwf[k][l][m]);
			  minval = min(minval,emwf[k][l][m]);
			  maxval = max(maxval,emwf[k][l][m]);
			}
		      if(runnumdraw > -1 && evtdraw > -1 && (maxval - minval > 30))
			{
			  singlewf->GetYaxis()->SetRangeUser(0,17000);
			  singlewf->Draw("HIST");
			  texts[4] = "EMCal #eta bin " + to_string(k) + ", #phi bin " + to_string(l);
			  tex[4] = new TLatex(0.1,ycoord[4],texts[4].c_str());
			  tex[4]->SetTextFont(42);
			  tex[4]->SetTextSize(0.04);
			  tex[4]->SetLineWidth(1);
			  tex[4]->SetNDC();
			  for(int m=0; m<5; ++m)
			    {
			      tex[m]->Draw();
			    }
			  gPad->SaveAs(("../images/emwf/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+to_string(k)+"_"+to_string(l)+"_emsinglewf_"+ct[failscut+1]+".png").c_str());
			  
			}
		      singlewf->Reset();
		    }
		}
	    }
	  else if(j==1)
	    {
	      for(int k=0; k<24; ++k)
		{
		  for(int l=0; l<64; ++l)
		    {
		      if(jconih[k][l] < jetcut) continue;
		      int minval = 9999999;
		      int maxval = -999999;
		      for(int m=0; m<12; ++m)
			{
			  if(ihwf[k][l][m] > 0) h2_wf[j]->Fill(m,ihwf[k][l][m]);
			  if(runnumdraw > -1 && evtdraw > -1) singlewf->Fill(m,ihwf[k][l][m]);
			  minval = min(minval,ihwf[k][l][m]);
			  maxval = max(maxval,ihwf[k][l][m]);
			}
		      if(runnumdraw > -1 && evtdraw > -1 && maxval-minval>10)
			{
			  singlewf->Draw("HIST");
			  texts[4] = "IHCal #eta bin " + to_string(k) + ", #phi bin " + to_string(l);
			  tex[4] = new TLatex(0.1,ycoord[4],texts[4].c_str());
			  tex[4]->SetTextFont(42);
			  tex[4]->SetTextSize(0.04);
			  tex[4]->SetLineWidth(1);
			  tex[4]->SetNDC();
			  for(int m=0; m<5; ++m)
			    {
			      tex[m]->Draw();
			    }
			  gPad->SaveAs(("../images/ihwf/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+to_string(k)+"_"+to_string(l)+"_ihsinglewf_"+ct[failscut+1]+".png").c_str());
			  
			}
		      singlewf->Reset();
		    }
		}
	    }
	  else if(j==2)
	    {
	      for(int k=0; k<24; ++k)
		{
		  for(int l=0; l<64; ++l)
		    {
		      //if(jconoh[k][l] < jetcut) continue;
		      int minval = 9999999;
		      int maxval = -999999;
		      for(int m=0; m<12; ++m)
			{
			  //cout << ohwf[k][l][m] << endl;
			  if(ohwf[k][l][m] > 0) h2_wf[j]->Fill(m,ohwf[k][l][m]);
			  if(runnumdraw > -1 && evtdraw > -1) singlewf->Fill(m,ohwf[k][l][m]);
			  minval = min(minval,ohwf[k][l][m]);
			  maxval = max(maxval,ohwf[k][l][m]);
			}
		      //cout << maxval << endl;
		      if(runnumdraw > -1 && evtdraw > -1 && maxval-minval>30)
			{
			  singlewf->Draw("HIST");
			  texts[4] = "OHCal #eta bin " + to_string(k) + ", #phi bin " + to_string(l);
			  tex[4] = new TLatex(0.1,ycoord[4],texts[4].c_str());
			  tex[4]->SetTextFont(42);
			  tex[4]->SetTextSize(0.04);
			  tex[4]->SetLineWidth(1);
			  tex[4]->SetNDC();
			  for(int m=0; m<5; ++m)
			    {
			      tex[m]->Draw();
			    }
			  gPad->SaveAs(("../images/ohwf/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+to_string(k)+"_"+to_string(l)+"_ohsinglewf_"+ct[failscut+1]+".png").c_str());
			  
			}
		      singlewf->Reset();
		    }
		}
	    }
	}

      for(int j=0; j<3; ++j)
	{
	  h2_wf[j]->Draw("COLZ");
	  for(int k=0; k<4; ++k)
	    {
	      tex[k]->Draw();
	    }
	  gPad->SaveAs(("../images/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+calo[j]+"_"+ct[failscut+1]+"_wf.png").c_str());
	  h2_wf[j]->Reset();
	}
      delete tex[2];
    }
  return 0;
}
