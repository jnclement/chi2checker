#include <../dlUtility.h>

int drawprettyeff(std::vector<TH3D*> nums3, TH3D* den3, std::vector<vector<int>> xbounds, std::vector<vector<int>> ybounds, std::vector<vector<int>> zbounds, int axis, std::vector<int> colors, std::vector<int> markers, std::vector<string> numlabels, std::vector<string> denlabels, string title)
{

  if(nums3.size() > colors.size() || nums3.size() > markers.size())
    {
      cout << "error! nums must be shorter vector than colors, markers!" << endl;
      return 1;
    }

  std::vector<TH1D*> nums1 = {};
  std::vector<TH1D*> dens1 = {};

  for(int i=0; i<nums3.size(); ++i)
    {
      if(axis==0)
	{
	  int ylo = ybounds.at(i).at(0);
	  int yhi = ybounds.at(i).at(1);
	  int zlo = zbounds.at(i).at(0);
	  int zhi = zbounds.at(i).at(1);
	  nums1.push_back(nums3.at(i)->ProjectionX((string(nums3.at(i)->GetName())+"_projx_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
	  dens1.push_back(den3->ProjectionX((string(den3->GetName())+"_projx_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
	}
      if(axis==1)
	{
	  int xlo = xbounds.at(i).at(0);
	  int xhi = xbounds.at(i).at(1);
	  int zlo = zbounds.at(i).at(0);
	  int zhi = zbounds.at(i).at(1);
	  nums1.push_back(nums3.at(i)->ProjectionX((string(nums3.at(i)->GetName())+"_projx_"+to_string(xlo)+"_"+to_string(xhi)+"_"+to_string(zlo)+"_"+to_string(xhi)).c_str(),xlo,xhi,zlo,zhi,"e"));
	  dens1.push_back(den3->ProjectionX((string(den3->GetName())+"_projx_"+to_string(xlo)+"_"+to_string(xhi)+"_"+to_string(zlo)+"_"+to_string(xhi)).c_str(),xlo,xhi,zlo,zhi,"e"));
	}
      if(axis==2)
	{
	  int xlo = xbounds.at(i).at(0);
	  int xhi = xbounds.at(i).at(1);
	  int ylo = ybounds.at(i).at(0);
	  int yhi = ybounds.at(i).at(1);
	  nums1.push_back(nums3.at(i)->ProjectionX((string(nums3.at(i)->GetName())+"_projx_"+to_string(xlo)+"_"+to_string(xhi)+"_"+to_string(ylo)+"_"+to_string(xhi)).c_str(),xlo,xhi,ylo,yhi,"e"));
	  dens1.push_back(den3->ProjectionX((string(den3->GetName())+"_projx_"+to_string(xlo)+"_"+to_string(xhi)+"_"+to_string(ylo)+"_"+to_string(xhi)).c_str(),xlo,xhi,ylo,yhi,"e"));
	}
    }

  std::vector<TH1D*> effs1 = {};

  for(int i=0; i<nums1.size(); ++i)
    {
      effs1.push_back(nums1.at(i)->Clone((string(nums1.at(i)->GetName())+"_eff").c_str()));
      effs1.at(i)->Divide(nums1.at(i),dens1.at(i),1,1,"B");
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  ratioPanelCanvas(can,0.3);
  can->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetRightMargin(0.05);

  TLegend* leg = new TLegend(0.6,0.5,0.9,0.8);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetNColumns(2);
  
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->SetMarkerColor(colors.at(i));
      nums1.at(i)->SetLineColor(colors.at(i));
      nums1.at(i)->SetMarkerStyle(markers.at(i));
      nums1.at(i)->SetMarkerSize(2);
      leg->AddEntry(nums1.at(i),numlabels.at(i),"p");

      effs1.at(i)->SetMarkerColor(colors.at(i));
      effs1.at(i)->SetLineColor(colors.at(i));
      effs1.at(i)->SetMarkerStyle(markers.at(i));
      effs1.at(i)->SetMarkerSize(2);

      dens1.at(i)->SetMarkerColor(colors.at(i)+2);
      dens1.at(i)->SetLineColor(colors.at(i)+2);
      dens1.at(i)->SetMarkerStyle(markers.at(i));
      dens1.at(i)->SetMarkerSize(2);
      leg->AddEntry(dens1.at(i),denlabels.at(i),"p");
    }

  for(int i=0; i<nums1.size(); ++i)
    {
      if(i==0) nums1.at(i)->Draw("PE");
      else nums1.at(i)->Draw("SAME PE");
      dens1.at(i)->Draw("SAME PE");
    }

  c->cd(2);
  for(int i=0; i<effs1.size(); ++i)
    {
      if(i==0) effs1.at(i)->Draw("PE");
      else effs1.at(i)->Draw("SAME PE");
    }

  leg->Draw();
  can->cd(0);

  maintexts(0.98,0.6,0,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;

  return 0;
}


int draw_cuteff()
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  
  
  return 0;
}
