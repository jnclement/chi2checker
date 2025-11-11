int draw_fake(int toad = 0)
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  TFile* f[4];
  string typestr[4] = {"dat","jet10","jet30","jet50"};

  for(int i=0; i<4; ++i)
    {
      f[i] = TFile::Open(("../hists_out_"+typestr[i]+".root").c_str(),"READ");
    }
  const int ntype = 2;
  
  TH1D* hspectra[ntype][20];
  TH1D* hleadtimeNOMBD[ntype];
  TH1D* hleadtimeYESMBD[ntype];
  TH1D* hleadtimeNOMBDwdijet[ntype];
  TH1D* hleadtimeYESMBDwdijet[ntype];
  TH1D* hleadtimeNOMBDwdijetP[ntype];
  TH1D* hleadtimeYESMBDwdijetP[ntype];
  TH1D* hratio[ntype][10];

  float scfs[3] = {3.997e6/7.2695,2.502e3/7.2695,1};
  
  for(int i=0; i<4; ++i)
    {
      int index = (i==0?0:1);
      if(i<2)
	{
	  hleadtimeNOMBD[index] = (TH1D*)(f[i]->Get(("hleadtimeNOMBD"+typestr[i]).c_str()));
	  if(i==1) hleadtimeNOMBD[index]->Scale(scfs[i-1]);
	  hleadtimeYESMBD[index] = (TH1D*)(f[i]->Get(("hleadtimeYESMBD"+typestr[i]).c_str()));
	  if(i==1) hleadtimeYESMBD[index]->Scale(scfs[i-1]);
	  hleadtimeNOMBDwdijet[index] = (TH1D*)(f[i]->Get(("hleadtimeNOMBDwdijet"+typestr[i]).c_str()));
	  if(i==1) hleadtimeNOMBDwdijet[index]->Scale(scfs[i-1]);
	  hleadtimeYESMBDwdijet[index] = (TH1D*)(f[i]->Get(("hleadtimeYESMBDwdijet"+typestr[i]).c_str()));
	  if(i==1) hleadtimeYESMBDwdijet[index]->Scale(scfs[i-1]);
	  hleadtimeNOMBDwdijetP[index] = (TH1D*)(f[i]->Get(("hleadtimeNOMBDwdijetP"+typestr[i]).c_str()));
	  if(i==1) hleadtimeNOMBDwdijetP[index]->Scale(scfs[i-1]);
	  hleadtimeYESMBDwdijetP[index] = (TH1D*)(f[i]->Get(("hleadtimeYESMBDwdijetP"+typestr[i]).c_str()));
	  if(i==1) hleadtimeYESMBDwdijetP[index]->Scale(scfs[i-1]);
	  for(int j=0; j<20; ++j)
	    {
	      hspectra[index][j] = (TH1D*)(f[i]->Get(("hspectra"+to_string(j)+typestr[i]).c_str()));
	      if(i==1) hspectra[index][j]->Scale(scfs[i-1]);
	      if(j<10)
		{
		  hratio[index][j] = (TH1D*)(f[i]->Get(("hratio"+to_string(j)+typestr[i]).c_str()));
		  if(i==1) hratio[index][j]->Scale(scfs[i-1]);
		}
	    }
	}
      else
	{
	  hleadtimeNOMBD[index]->Add((TH1D*)(f[i]->Get(("hleadtimeNOMBD"+typestr[i]).c_str())),scfs[i-1]);
	  hleadtimeYESMBD[index]->Add((TH1D*)(f[i]->Get(("hleadtimeYESMBD"+typestr[i]).c_str())),scfs[i-1]);
	  hleadtimeNOMBDwdijet[index]->Add((TH1D*)(f[i]->Get(("hleadtimeNOMBDwdijet"+typestr[i]).c_str())),scfs[i-1]);
	  hleadtimeYESMBDwdijet[index]->Add((TH1D*)(f[i]->Get(("hleadtimeYESMBDwdijet"+typestr[i]).c_str())),scfs[i-1]);
	  hleadtimeNOMBDwdijetP[index]->Add((TH1D*)(f[i]->Get(("hleadtimeNOMBDwdijetP"+typestr[i]).c_str())),scfs[i-1]);
	  hleadtimeYESMBDwdijetP[index]->Add((TH1D*)(f[i]->Get(("hleadtimeYESMBDwdijetP"+typestr[i]).c_str())),scfs[i-1]);
	  for(int j=0; j<20; ++j)
	    {
	      hspectra[index][j]->Add((TH1D*)(f[i]->Get(("hspectra"+to_string(j)+typestr[i]).c_str())),scfs[i-1]);
	    }
	}
    }

  TCanvas* c = new TCanvas("","",1000,1000);
  c->cd();
  for(int i=0; i<ntype; ++i)
    {
      hratio[i][0]->Divide(hspectra[i][7+toad],hspectra[i][6+toad],1,1,"B");
      hratio[i][1]->Divide(hspectra[i][9+toad],hspectra[i][8+toad],1,1,"B");
    }

  TLegend *tleg = new TLegend(0.13,0.13,0.47,0.35,"","brNDC");
  tleg->SetLineWidth(0);
  tleg->SetFillStyle(0);
  int markers[4] = {20,21,71,72};
  int colors[4] = {kRed+2,kSpring+2,kOrange+2,kMagenta+2};
  for(int i=0; i<ntype; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  hratio[i][j]->SetMarkerSize(1.5);
	  hratio[i][j]->SetMarkerStyle(markers[2*i+(j>1?j-2:j)]);
	  hratio[i][j]->SetMarkerColor(colors[2*i+(j>1?j-2:j)]);
	  hratio[i][j]->SetLineColor(colors[2*i+(j>1?j-2:j)]);
	  hratio[i][j]->GetYaxis()->SetRangeUser(0,1.05);
	  hratio[i][j]->GetXaxis()->SetRangeUser(0,60);
	  if(i==0 && j==0) hratio[i][j]->Draw("PE");
	  else hratio[i][j]->Draw("SAME PE");
	}
    }
  
  tleg->AddEntry(hratio[0][0],"Dijet\&t\&MBD / MBD (Data)","p");
  tleg->AddEntry(hratio[1][0],"Dijet\&t\&MBD / MBD (Sim)","p");
  tleg->AddEntry(hratio[0][1],"Dijet\&t\&MBD / Dijet\&t (Data)","p");
  tleg->AddEntry(hratio[1][1],"Dijet\&t\&MBD / Dijet\&t (Sim)","p");
  tleg->Draw();
  string toadstr = (toad?"frac":"nofrac");
  c->SaveAs(("ratios"+toadstr+".pdf").c_str());
  return 0;
}
