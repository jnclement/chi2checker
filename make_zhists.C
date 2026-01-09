int make_zhists(string filename)
{
  TFile* f=TFile::Open(filename.c_str(),"READ");
  TTree* outt = (TTree*)f->Get("outt");

  float zvtx;
  int run;

  outt->SetBranchStatus("*",0);
  outt->SetBranchStatus("zvtx",1);
  outt->SetBranchStatus("run",1);

  outt->SetBranchAddress("run",&run);
  outt->SetBranchAddress("zvtx",&zvtx);

  TH1D* zzmrad = new TH1D("zzmrad",";z_{vtx} [cm];Counts",600,-300,300);
  TH1D* nzmrad = new TH1D("nzmrad",";z_{vtx} [cm];Counts",600,-300,300);

  for(int i=0; i<outt->GetEntries(); ++i)
    {
      outt->GetEntry(i);
      if(zvtx==0) continue;
      if(run < 51212) zzmrad->Fill(zvtx);
      else nzmrad->Fill(zvtx);
    }

  TFile* outf=TFile::Open("outzhists.root","RECREATE");
  outf->cd();

  zzmrad->Write();
  nzmrad->Write();

  outf->Close();

  return 0;
}
