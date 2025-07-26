void drawText(const char *text, float xp, float yp, bool isRightAlign=0, int textColor=kBlack, double textSize=0.04, int textFont = 42, bool isNDC=true){
  // when textfont 42, textSize=0.04                                                                                                                         
  // when textfont 43, textSize=18                                                                                                                           
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(textFont);
  //   if(bold)tex->SetTextFont(43);                                                                                                                         
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  if(isNDC) tex->SetNDC();
  if(isRightAlign) tex->SetTextAlign(31);
  tex->Draw();
  //delete tex;
}

void sqrt_snn_text(float xp = 0.7, float yp = 0.8, bool isRightAlign=0, double textsize = 0.04)
{
  drawText("#sqrt{S_{NN}} = 200 GeV",xp,yp,isRightAlign,kBlack,textsize);
}

void sphenixtext(float xpos = 0.7, float ypos = 0.96, int ra = 0, float textsize = 0.04)
{
  drawText("#bf{#it{sPHENIX}} Internal", xpos, ypos, ra, kBlack, textsize);
}

void drawCalo(float** towersem, float** towersih, float** towersoh, float* jet_e, float* jet_et, float* jet_ph, int jet_n, float zvtx, int failscut, int runnum, int evtnum, float frcoh, float frcem)
{
  std::stringstream full_stream;
  full_stream << std::fixed << std::setprecision(3) << "Run " << runnum << ", event " << evtnum << ", EM fraction: " << frcem << ", OH fraction: " << frcoh << ", z_{vtx} = " <<  std::setprecision(1) << zvtx << " ";
  std::string full_string = full_stream.str();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetNumberContours(50);
  int ncircle = 64;
  const int ncol = 50;
  const int nstp = 3;
  double red[nstp] = {0.,1.,1.};// = {1, 1, 1, 1, 1, 1, 1, 1};
  double grn[nstp] = {0.,1.,0.};// = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double blu[nstp] = {1.,1.,0.};// = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double stp[nstp] = {0,2./27,1};;// = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  

  
  TCanvas* c = new TCanvas("","",1900,1000);
  c->Divide(4,1,0,0.1);
  
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,-0.5,23.5,64,-0.5,63.5);
  TH2D* event_disrt[3];
  for(int i=0; i<3; ++i)
    {
      int nbinx = (i==0?96:24);
      int nbiny = (i==0?256:64);
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);

    }
  TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta Bin");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi Bin");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta Bin");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi Bin");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta Bin");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi Bin");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->GetZaxis()->SetTitleOffset(1);
      event_disrt[i]->GetYaxis()->SetTitleOffset(1);
      event_disrt[i]->GetZaxis()->SetTitle("Tower Energy [GeV]");
      event_disrt[i]->GetZaxis()->SetRangeUser(-2,25);
      //event_disrt[i]->GetXaxis()->SetNdivisions(4,kFALSE);
      event_disrt[i]->GetXaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetYaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetZaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetXaxis()->SetTitleOffset(1);
      event_disrt[i]->GetXaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetYaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetZaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetXaxis()->SetLabelOffset(0.02);
    }
  event_sum->GetXaxis()->SetTitle("Calo Sum #eta Bin");
  event_sum->GetYaxis()->SetTitle("Calo Sum #phi Bin");
  event_sum->GetZaxis()->SetTitleOffset(1);
  event_sum->GetYaxis()->SetTitleOffset(1);
  event_sum->GetZaxis()->SetTitle("Tower Energy [GeV]");
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  //event_sum->GetXaxis()->SetNdivisions(4,kFALSE);
  event_sum->GetXaxis()->SetTitleSize(0.04);
  event_sum->GetYaxis()->SetTitleSize(0.04);
  event_sum->GetZaxis()->SetTitleSize(0.04);
  event_sum->GetXaxis()->SetTitleOffset(1);
  event_sum->GetXaxis()->SetLabelSize(0.04);
  event_sum->GetYaxis()->SetLabelSize(0.04);
  event_sum->GetZaxis()->SetLabelSize(0.04);
  event_sum->GetXaxis()->SetLabelOffset(0.02);
  event_sum->Reset();
  for(int j=0; j<3; ++j)
    {
      event_disrt[j]->Reset();
      for(int eta=0; eta<(j==0?96:24); ++eta)
	{
	  for(int phi = 0; phi<(j==0?256:64); ++phi)
	    {
	      float energy = 0;
	      if(j==0) energy = towersem[eta][phi];
	      if(j==1) energy = towersih[eta][phi];
	      if(j==2) energy = towersoh[eta][phi];
	      event_disrt[j]->Fill(eta,phi,);
	      if(j==0)event_sum->Fill(eta/4,phi/4,tower->get_energy());
	      else event_sum->Fill(eta,phi,tower->get_energy());
	    }
	}

      c->cd(j+1);
      gPad->SetLogz(0);
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.2);
      gPad->SetTopMargin(0.05);
      event_disrt[j]->SetContour(ncol);
      event_disrt[j]->Draw("COLZ");

    }
  c->cd(4);
  gPad->SetLogz(0);            
  gPad->SetTopMargin(0.05);                                                       
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  event_sum->SetContour(ncol);
  event_sum->Draw("COLZ");
  c->cd(0);
  sphenixtext(0.96,0.96,1,0.04);
  
  std::string fails = "";
  std::string whichcut = "";
  if(failscut == 0)
    {
      fails = "Passes dijet cut";
      whichcut = "dijet";
    }
  else if(failscut == 1)
    {
      fails = "Passes frac. cut";
      whichcut = "frac";
    }
  else
    {
      fails = "Passes both";
      whichcut = "both";
    }
  full_string += fails;
  
  //drawText(fails.c_str(),0.7,0.95);
  drawText(full_string.c_str(),0.05,0.95,0,kBlack,0.02);
  c->cd(4);
  maxJetE = 0;
  for(int k=0; k<jet_n; ++k)
    {
      if(maxJetE < jet_e[k])
	{
	  maxJetE = jet_e[k];
	}
      std::stringstream e_stream;
      e_stream << std::fixed << std::setprecision(2) << jet_e[k];
      std::string e_string = e_stream.str();
      drawText((e_string+" GeV").c_str(),12,((jet_ph[k]/*+(jet_ph[k]+M_PI>3.84?-0.53:0.43)*/+(jet_ph[k]<0?2*M_PI:0))/(2*M_PI))*64,/*(jet_et[k]>0?1:*/0/*)*/,kBlack,0.04,42,false);

    }
  //c->Update();
  string dirstring = "";
  for(int i=46; i<101; ++i)
    {
      if(maxJetE < i)
	{
	  dirstring = to_string(i-1)+"to"+to_string(i);
	  break;
	}
    }
  if(maxJetE > 100) dirstring = "gr100";
      
  c->SaveAs(("/sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/candidate_"+dirstring+"_"+_name+"_"+whichcut+"_"+to_string(cancount)+".png").c_str());
  cout << "Saved" << endl;

  for(int i=0; i<3; ++i)
    {
      c->cd(i+1);
      gPad->SetLogz();
      event_disrt[i]->GetZaxis()->SetRangeUser(0.05,25);
      gPad->Update();
    }
  c->cd(4);
  gPad->SetLogz();
  event_sum->GetZaxis()->SetRangeUser(0.05,25);
  gPad->Update();
  c->SaveAs(("/sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/candidate_"+dirstring+"_"+_name+"_"+whichcut+"_"+to_string(cancount)+"_log.png").c_str());
  if(c) delete c;
  if(event_disrt[0]) delete event_disrt[0];
  if(event_disrt[1]) delete event_disrt[1];
  if(event_disrt[2]) delete event_disrt[2];
  if(event_sum) delete event_sum;
  
}
