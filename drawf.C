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


  gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;");
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
  float emtow[96][256];
  float ihtow[24][64];
  float ohtow[24][64];
  float jet_e[100];
  int jet_n;
  
  TChain* wft = new TChain("wft");
  TChain* jt = new TChain("jet_tree");

  const float jetcut = 30;

  float jetsum[3] = {0};
  
  
  string calo[3] = {"EMCal","IHCal","OHCal"};
  string ct[4] = {"fails","dijet","frac","both"};


  string tempinfilename;
  ifstream chi2list("chi2files.txt");
  for(int i=0; i<lo; ++i)
    {
      std::getline(chi2list,tempinfilename);
    }
  int counter = lo;
  while(std::getline(chi2list,tempinfilename))
    {
      jt->Add(tempinfilename.c_str());
      ++counter;
      if(counter >= hi) break;
    }

  ifstream wflist("wavefiles.txt");
  for(int i=0; i<lo; ++i)
    {
      std::getline(wflist,tempinfilename);
    }
  counter = lo;
  while(std::getline(wflist,tempinfilename))
    {
      wft->Add(tempinfilename.c_str());
      ++counter;
      if(counter >= hi) break;
    }


  jt->SetBranchAddress("jconem",jconem);
  jt->SetBranchAddress("jconih",jconih);
  jt->SetBranchAddress("jconoh",jconoh);
  jt->SetBranchAddress("runnum",&rnj);
  jt->SetBranchAddress("evtnum",&enj);
  jt->SetBranchAddress("emtow",emtow);
  jt->SetBranchAddress("ihtow",ihtow);
  jt->SetBranchAddress("ohtow",ohtow);
  jt->SetBranchAddress("jet_n",&jet_n);
  jt->SetBranchAddress("jet_et",jet_e);
  
  wft->SetBranchAddress("runnum",&runnum);
  wft->SetBranchAddress("evtnum",&evtnum);
  wft->SetBranchAddress("emwf",emwf);
  wft->SetBranchAddress("ihwf",ihwf);
  wft->SetBranchAddress("ohwf",ohwf);
  wft->SetBranchAddress("failscut",&failscut);

  
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
  int jte = 0;
  
  for(int i=0; i<wft->GetEntries(); ++i)
    {
      if(i%100==0) cout << i << endl;
      wft->GetEntry(i);
      jt->GetEntry(jte);
      int flag = 0;
      while(rnj != runnum || enj != evtnum)
	{
	  ++jte;
	  if(jte == jt->GetEntries())
	    {
	      flag = 1;
	      jte = 0;
	      break;
	    }
	  jt->GetEntry(jte);
	}
      if(flag) continue;
      if((failscut < 0 || failscut > 2)) continue;
      float maxE = 0;
      for(int j=0; j<jet_n; ++j)
	{
	  if(jet_e[j] > maxE) maxE = jet_e[j];
	}
      if(runnumdraw >= 0)
	{
	  if(runnum != runnumdraw) continue;
	  if(evtdraw >= 0)
	    {
	      if(evtnum != evtdraw) continue;
	    }
	}
      if(maxE < 60) continue;
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
			  singlewf->Fill(m,emwf[k][l][m]);//if(runnumdraw > -1 && evtdraw > -1)
			  minval = min(minval,emwf[k][l][m]);
			  maxval = max(maxval,emwf[k][l][m]);
			}
		      if(emtow[k][l] > 0)
			{
			  jetsum[j] += emtow[k][l];
			  //cout << emtow[k][l] << " " << jetsum[j] << " " << (maxval-minval)/emtow[k][l] << endl;
			}
		      //if(runnumdraw > -1 && evtdraw > -1 && (maxval-minval > 30))
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
			  TLatex* Etex = new TLatex(0.1,0.86,("Tower E="+to_string(emtow[k][l])+" GeV").c_str());
			  Etex->SetTextFont(42);
			  Etex->SetTextSize(0.04);
			  Etex->SetTextColor(kBlack);
			  Etex->SetLineWidth(1);
			  Etex->SetNDC();
			  Etex->Draw();
			  
			  gPad->SaveAs(("../images/wf/emwf/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+to_string(k)+"_"+to_string(l)+"_emsinglewf_"+ct[failscut+1]+".png").c_str());
			  
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
		      if(ihtow[k][l] > 0)
			{
			  jetsum[j] += ihtow[k][l];
			  //cout << ihtow[k][l] << " " << jetsum[j] << " " << (maxval-minval)/ihtow[k][l] << endl;
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
			  TLatex* Etex = new TLatex(0.1,0.86,("Tower E="+to_string(ihtow[k][l])+" GeV").c_str());
			  Etex->SetTextFont(42);
			  Etex->SetTextSize(0.04);
			  Etex->SetTextColor(kBlack);
			  Etex->SetLineWidth(1);
			  Etex->SetNDC();
			  Etex->Draw();
			  
			  gPad->SaveAs(("../images/wf/ihwf/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+to_string(k)+"_"+to_string(l)+"_ihsinglewf_"+ct[failscut+1]+".png").c_str());
			  
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
		      if(jconoh[k][l] < jetcut) continue;
		      
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
		      if(ohtow[k][l] > 0)
			{		       
			  jetsum[j] += ohtow[k][l];
			  //cout << ohtow[k][l] << " " << jetsum[j] << " " << (maxval-minval)/ohtow[k][l] << endl;
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
			  TLatex* Etex = new TLatex(0.1,0.86,("Tower E="+to_string(ohtow[k][l])+" GeV").c_str());
			  Etex->SetTextFont(42);
			  Etex->SetTextSize(0.04);
			  Etex->SetTextColor(kBlack);
			  Etex->SetLineWidth(1);
			  Etex->SetNDC();
			  Etex->Draw();
			  
			  gPad->SaveAs(("../images/wf/ohwf/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+to_string(k)+"_"+to_string(l)+"_ohsinglewf_"+ct[failscut+1]+".png").c_str());
			  
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
	  gPad->SaveAs(("../images/wf/"+to_string(runnum)+"_"+to_string(evtnum)+"_"+calo[j]+"_"+ct[failscut+1]+"_wf.png").c_str());
	  h2_wf[j]->Reset();
	}
      delete tex[2];
    }
  for(int i=0; i<3; ++i)
    {
      cout << jetsum[i] << endl;
    }
  return 0;
}
