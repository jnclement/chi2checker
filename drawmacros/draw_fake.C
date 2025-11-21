int draw_fake(int toad = 0, int mbdtimereq = 0)
{
  string toadstr = (toad?"frac":"");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  TFile* f[4];
  string typestr[4] = {"dat","jet10","jet30","jet50"};
  string mbdtimereqstr = (mbdtimereq?"_mbdtimereq":"");
  for(int i=0; i<4; ++i)
    {
      typestr[i] += toadstr;
      cout << typestr[i] << endl;
      f[i] = TFile::Open(("../hists"+mbdtimereqstr+"_out_"+typestr[i]+".root").c_str(),"READ");
    }
  const int ntype = 2;
  
  TH1D* hspectra[ntype][20];
  TH1D* hleadtimeNOMBD[ntype];
  TH1D* hleadtimeYESMBD[ntype];
  TH1D* hleadtimeNOMBDwdijet[ntype];
  TH1D* hleadtimeYESMBDwdijet[ntype];
  TH1D* hleadtimeNOMBDwdijetP[ntype];
  TH1D* hleadtimeYESMBDwdijetP[ntype];
  TH1D* hratio[ntype][13];
  TH1D* hfrac[2];
  hfrac[0] = (TH1D*)(f[0]->Get(("hfracgood0dat")));
  hfrac[1] = (TH1D*)(f[0]->Get(("hfracbad0dat")));
  float scfs[3] = {3.997e6/7.2695,2.502e3/7.2695,1};
  
  for(int i=0; i<4; ++i)
    {
      cout << i << endl;
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
	      cout << j << endl;
	      hspectra[index][j] = (TH1D*)(f[i]->Get(("hspectra"+to_string(j)+typestr[i]).c_str()));
	      if(i==1) hspectra[index][j]->Scale(scfs[i-1]);
	      if(j<10)
		{
		  hratio[index][j] = (TH1D*)(f[i]->Get(("hratio"+to_string(j)+typestr[i]).c_str()));
		  if(i==1) hratio[index][j]->Scale(scfs[i-1]);
		}
	      else if(j<13)
		{
		  hratio[index][j] = (TH1D*)(hratio[index][9]->Clone(("hratio"+to_string(j)+typestr[i]).c_str()));
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
	      cout << j << endl;
	      hspectra[index][j]->Add((TH1D*)(f[i]->Get(("hspectra"+to_string(j)+typestr[i]).c_str())),scfs[i-1]);
	    }
	}
    }

  TCanvas* c = new TCanvas("","",1000,1000);
  c->cd();
  int todivide[13] = {0,0,2,2,2,0,6,6,8,8,10,10,10};
  for(int i=0; i<ntype; ++i)
    {
      for(int j=0; j<13; ++j)
	{
	  cout<< i << " " << j << endl;
	  hratio[i][j]->Divide(hspectra[i][j],hspectra[i][todivide[j]],1,1,"B");
	}
    }
  cout << "madeit" << endl;
  string ytitles[13] = {"Should not be drawn","Dijet-dt \& t \& MBD-t / Dijet dt \& t","Should not be drawn","MBD-coinc \& MBD-t \& Dijet dt \& t / MBD-coinc \& Dijet dt \& t","MBD-coinc \& MBD-t \& t \& Dijet dt / MBD-coinc \& Dijet dt \& t","Dijet-dt \& t \& MBD-t / Dijet dt \& t","MBD-coinc \& MBD-t \& Dijet \& t / MBD-coinc \& MBD-t","MBD-coinc \& MBD-t \& Dijet-dt \& t / MBD-coinc \& MBD-t","Should not be drawn","Dijet-dt \& t \& MBD-coinc / Dijet-dt \& t","Should not be drawn","t \& MBD-coinc / t","t \& Dijet-dt / t"};
  TLegend *tleg = new TLegend(0.11,0.13,0.23,0.23,"","brNDC");
  tleg->SetLineWidth(0);
  tleg->SetFillStyle(0);
  int markers[2] = {20,72};
  int colors[2] = {kSpring+2,kMagenta+2};
  for(int j=0; j<13; ++j)
    {
      if(j==0 || j==2 || j==8 || j==6 || j==10) continue;
      for(int i=0; i<ntype; ++i)
	{
	  cout << i << " " << j << endl;
	  hratio[i][j]->SetMarkerSize(1.5);
	  hratio[i][j]->SetMarkerStyle(markers[i]);
	  hratio[i][j]->SetMarkerColor(colors[i]);
	  hratio[i][j]->SetLineColor(colors[i]);
	  hratio[i][j]->GetYaxis()->SetRangeUser(0,1.05);
	  hratio[i][j]->GetXaxis()->SetRangeUser(0,60);
	  hratio[i][j]->GetYaxis()->SetTitle(("Ratio of "+ytitles[j]).c_str());
	  if(j==1 && i==0)
	    {
	      tleg->AddEntry(hratio[0][1],"Data","p");
	      tleg->AddEntry(hratio[1][1],"Sim","p");
	    }
	  if(i==0) hratio[i][j]->Draw("PE");
	  //else if(j!=1) hratio[i][j]->Draw("SAME PE");

	}
      //if(j!=1) tleg->Draw();
      c->SaveAs(("ratios"+mbdtimereqstr+"_"+toadstr+to_string(j)+".pdf").c_str());
    }
  

  //tleg->AddEntry(hratio[0][1],"Dijet\&t\&MBD / Dijet\&t (Data)","p");
  //tleg->AddEntry(hratio[1][1],"Dijet\&t\&MBD / Dijet\&t (Sim)","p");


  TH3D* htdtmbdt;
  htdtmbdt = (TH3D*)(f[0]->Get("htdtmbdtdat"));
  htdtmbdt->GetZaxis()->SetTitle("MBD - t_{jet} [ns]");

  TH1D* projz = htdtmbdt->ProjectionZ("pz",111,170,136,165);

  projz->SetMarkerStyle(20);
  projz->SetMarkerColor(kGreen+2);
  projz->SetLineColor(kGreen+2);
  projz->GetYaxis()->SetRangeUser(0,2e6);
  projz->GetYaxis()->SetTitle("Counts");
  c->SetLeftMargin(0.15);
  projz->Draw("PE");

  TLine* line1 = new TLine(-3,0,-3,2e6);
  TLine* line2 = new TLine(3,0,3,2e6);
  line1->Draw();
  line2->Draw();

  c->SaveAs("mbdt.pdf");
  cout << (projz->Integral(1,135) + projz->Integral(166,300))/projz->Integral(1,300) << endl;

  hfrac[1]->Divide(hfrac[1],hfrac[0],1,1,"B");
  hfrac[1]->SetMarkerStyle(20);
  hfrac[1]->SetMarkerSize(1);
  hfrac[1]->SetLineColor(kRed);
  hfrac[1]->GetYaxis()->SetTitle("Fraction in Negative Side Tail");
  hfrac[1]->GetXaxis()->SetTitle("Run Number");
  hfrac[1]->Draw("PE");
  hfrac[1]->GetYaxis()->SetTitleOffset(1.5);
  c->SaveAs("frachist.pdf");
  return 0;
}
