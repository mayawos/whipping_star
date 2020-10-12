#include "peleeplotter.h"
using namespace sbn;


void DrawUniverses(std::vector<TH1D*> histos, std::string name, std::string title, std::string var ){
  
  TCanvas can(name.c_str(),name.c_str());
  histos[0]->SetTitle(title.c_str());
  histos[0]->GetXaxis()->SetTitle(var.c_str());
  histos[0]->GetYaxis()->SetTitle("Events / 0.1 GeV");
  histos[0]->SetMaximum(1.2*histos[0]->GetMaximum());
  histos[0]->SetLineColor(1);
  histos[1]->SetLineColor(2);
  histos[2]->SetLineColor(2);
  //histos[1]->SetLineStyle(10);
  histos[2]->SetLineStyle(2);
  for(int i=0; i < 3; i++){
    histos[i]->SetLineWidth(2);
    if(i==0)histos[i]->Draw("hist");
    histos[i]->Draw("histsame");
  }
  TLegend  *legend = new TLegend(0.70,0.70,0.98,0.95); // we need different positions for the legend to not 
  // get the plot titles for the legend
  std::vector<std::string> legends_str = {"CV","Univ 0","Univ 1"};
  for(int i = 0; i < 3 ; i++){
    legend->AddEntry(histos[i],legends_str[i].c_str(),"l"); 
  }
  legend->Draw("same");
  can.Print(Form("%s_universes.pdf",name.c_str()),"pdf");
  
}

//borrow this function from SBNchi
void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels){
  bool debug = true;
  //if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<std::endl;
  //if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<std::endl;
  
  std::vector<std::vector<TMatrixD>> Summed(num_channels, std::vector<TMatrixD>(num_channels) );	//Initialise a matrix of matricies, to ZERO.
  for(int ic = 0; ic < num_channels; ic++){
    for(int jc =0; jc < num_channels; jc++){
      Summed[ic][jc].ResizeTo(num_bins[jc],num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
      Summed[ic][jc] = 0.0;
    }
  }
  
  int mrow = 0.0;
  int mcol = 0.0;
  
  for(int ic = 0; ic < num_channels; ic++){ 	 //Loop over all rows
    for(int jc =0; jc < num_channels; jc++){ //Loop over all columns
      
      //if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;
      
      for(int m=0; m < num_subchannels[ic]; m++){
	for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing toGether
	  //std::cout << mrow << ", " << n << ", " << num_bins[jc] << ", " << ", " << mcol << ", " << m << ", " << m << ", " << num_bins[ic] << std::endl;
	  Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
	  //std::cout << "Summed[" << ic << "][" << jc << "] = " << std::endl;
	  //Summed[ic][jc].Print();
	}
      }
      mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
    }//end of column loop
    
    mrow = 0; // as we end this row, reSet row count, but jump down 1 column
    mcol += num_subchannels[ic]*num_bins[ic];
  }//end of row loop
  
  /// ********************************* And put them back toGether! ************************ //
  Mc.Zero();
  mrow = 0;
  mcol = 0;
  
  //Repeat again for Contracted matrix
  for(int ic = 0; ic < num_channels; ic++){
    for(int jc =0; jc < num_channels; jc++){
      Mc.SetSub(mrow,mcol,Summed[ic][jc]);
      mrow += num_bins[jc];
    }
    
    mrow = 0;
    mcol +=num_bins[ic];
  }
  
  return;
}

//Pretty up the Plots!
void Draw_Stacked(TObjArray histos, TH1D *data, 
		  TString samplename,
		  TPad *pad,
		  bool normalised,
		  std::string stacktitle,
		  std::string var )
{// this function draws histograms stacked and correctly takes into account the
  // stats boxes for each
  gStyle->SetOptStat(0); 
  if( histos.GetEntries() == 0 ) return; // nothing to do
  
  // Initial set up
  TObjArray statsboxes;    
  std::vector<std::string> legends_str;
  std::vector<int> colours;
  TPaveStats *st1;
  
  const int n = histos.GetEntries();
  
  // choose the colours
  colours.push_back(kRed-3);
  colours.push_back(kOrange-5);
  colours.push_back(kMagenta-6);
  
  // lets open and draw the canvas 
  TCanvas *canvas;
  TString cName = Form("%s",samplename.Data());
  if(pad == 0){
    canvas = new TCanvas(cName,"Stacked Histograms", 1200., 800.);
    pad = (TPad*)canvas->cd();
  }
  pad->cd();
  pad->SetTicks(0,0);
  pad->SetRightMargin(0.1);
  pad->SetBottomMargin(0.13);
  
  // lets take the first histoname and see if we can match a title to its which will be HS stack title
  if( stacktitle == "" ) stacktitle = ((TH1F*)histos[0])->GetTitle();
  THStack *Hs = new THStack("hs2",stacktitle.c_str());
  
  //Set Axis Units
  //Hs->GetXaxis()->SetTitle("Number of Hits in Slice");
  
  // Set up the LEGEND
  TLegend  *legend = new TLegend(0.70,0.5,0.98,0.7); // we need different positions for the legend to not 
  // get the plot titles for the legend
  for(int i = 0; i < n ; i++){
    
    TH1D *h =  (TH1D*)histos[i];
    std::cout << "histo = " << h->GetTitle() << std::endl;
    legends_str.push_back( h->GetTitle() );
    //h->Scale(1,"width");
    h->SetLineWidth(0);
    h->SetLineColor(colours[i]);
    h->SetFillColor(colours[i]);
    legend->AddEntry(h,legends_str[i].c_str(),"f"); 

    std::cout << "histo name =  " << h->GetTitle() << std::endl;
    std::cout << "integral =  " << h->Integral() << std::endl;
    Hs->Add(h,"sames");
  }
  
  float heightboxes;
  // the position of the top corner
  float top_corner,deltay;
  // don't draw errors if it is 1
  int no_error = 1;
  top_corner = 0.9;

  /* 
  // if the array has more than 0 histograms lets specify the stat boxes 
  if (histos.GetEntries() > 0) { 
  heightboxes = (float)0.5/(float)histos.GetEntries();
  }
  else
  heightboxes = (float)0.5;
  */
  
  //HACK -need to draw a smaller histogram first to set the x-axis range correctly 
  //for stacks with variable-width bins
  //TH1D* tmp_mnv = (TH1D*)Hs->GetHists()->At(0);
  Hs->SetMaximum(1.4*data->GetMaximum());
  
  data->SetMarkerStyle(20);
  data->SetMarkerSize(2);
  data->SetMarkerColor(kBlack);
  data->SetLineWidth(2);
  data->SetLineColor(kBlack);
  
  data->Draw("PE1"); 
  // DRAW not stacked to get correct stats boxes
  if (no_error == 1) // do not draw errors
    Hs->Draw("HIST");
  else
    Hs->Draw("HISTE");

  Hs->GetXaxis()->SetTitle( var.c_str() );
  Hs->GetYaxis()->SetTitle( "Events" );
  //Hs->GetYaxis()->SetNdivisions(220, kTRUE);
  //Hs->GetXaxis()->SetNdivisions(220, kTRUE);
  Hs->GetXaxis()->SetTitleFont(43);
  Hs->GetYaxis()->SetTitleFont(43);
  Hs->GetXaxis()->SetTitleSize(40);
  Hs->GetYaxis()->SetTitleSize(40);
  Hs->GetXaxis()->SetLabelFont(43);
  Hs->GetYaxis()->SetLabelFont(43);
  Hs->GetXaxis()->SetLabelSize(25);
  Hs->GetYaxis()->SetLabelSize(25);
  Hs->GetXaxis()->CenterTitle(kTRUE);

  data->DrawCopy("PE1 same");
  legend->Draw("");
  
  gPad->Modified(); // so it updates the pad with the new changes
 
  canvas->Modified(); 
  canvas->Print(Form("%s.pdf",cName.Data()),"pdf"); 
  
  delete canvas;
  
  return;
}

TH1D* GetTotalError(std::string name, TMatrixD errMatrix, TH1D *histo ){
  // Make a copy of this histogram as a TH1D and rename it
  TString hname = name+"_total_error";
  TH1D *err = (TH1D*)histo->Clone(hname);
  err->Reset();
  
  const int highBin = err->GetNbinsX() + 1;
  const int lowBin = 0;
  
  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    {
      double derr = errMatrix[iBin][iBin];
      err->SetBinContent( iBin+1, ( derr > 0 ) ? sqrt(derr): 0. );
    }
  
  return err;
}

void DrawMCAndSyst(TString cName, TH1D* MC, std::string var, std::string title){

  TCanvas can(cName,cName,1200,800);

  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
  TH1D* sysErr = (TH1D*)MC->Clone("tempSysCV");
  TH1D* tmpsys = (TH1D*)MC->Clone("tempSysRatio");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC->GetBinError(bin)/tmpMC->GetBinContent(bin));
    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
    sysErr->SetBinError(bin,tmpMC->GetBinError(bin));
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  tmpMC->SetLineColor(kRed);
  tmpMC->SetLineWidth(3);
  tmpMC->SetLineStyle(1);
  sysErr->SetTitleSize(80);
  sysErr->SetTitleFont(43);
  sysErr->SetTitle(title.c_str());
  sysErr->SetLineColor(kRed-10);
  sysErr->SetLineWidth(1);
  sysErr->SetFillColor(kRed-10);
  sysErr->SetFillStyle(1001);
  sysErr->SetStats(0);
  sysErr->GetYaxis()->SetTitleSize(25);
  sysErr->GetYaxis()->SetTitleFont(43);
  sysErr->GetYaxis()->SetTitle("Events");
  sysErr->GetXaxis()->SetLabelSize(0);
  sysErr->GetXaxis()->SetTitleFont(43);
  sysErr->GetXaxis()->SetTitle(var.c_str());
  sysErr->SetMaximum(1.2*sysErr->GetMaximum());
  sysErr->Draw("E2");
  tmpMC->Draw("HIST same");

  // lower plot will be in pad
  can.cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  tmpsys->SetTitle("");
  tmpsys->GetYaxis()->SetTitle("uncertainties");
  tmpsys->GetXaxis()->SetTitle(var.c_str());
  tmpsys->SetStats(0);
  tmpsys->SetLineColor(kRed);
  tmpsys->SetLineWidth(1);
  tmpsys->SetFillColor(kRed-10);
  tmpsys->SetFillStyle(1001);
  tmpsys->GetYaxis()->SetTitleSize(25);
  tmpsys->GetYaxis()->SetTitleFont(43);
  tmpsys->GetXaxis()->SetTitleSize(25);
  tmpsys->GetXaxis()->SetTitleOffset(3.5);
  tmpsys->GetYaxis()->SetTitleOffset(1.0);
  tmpsys->GetXaxis()->SetLabelSize(0.1);
  tmpsys->GetYaxis()->SetLabelSize(0.08);
  tmpsys->GetYaxis()->SetNdivisions(210, kTRUE);
  tmpsys->GetXaxis()->SetTitleFont(43);
  tmpsys->SetMaximum(1.4);
  tmpsys->SetMinimum(0.6);
  tmpsys->Draw("E2");

  const TAxis *axis = tmpsys->GetXaxis();
  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
  double highX = axis->GetBinUpEdge(  axis->GetLast() );

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(kRed);
  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad
 
  
  pad2->Modified(); // so it updates the pad with the new changes
  pad2->Draw("");

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}

void DrawMCAndSystOverlay(TString cName, TH1D* MC1, TH1D* MC2, std::string var, std::string title){

  TCanvas can(cName,cName,1200,800);

  TH1D* tmpMC = (TH1D*)MC1->Clone("tempCV");
  TH1D* tmpMC2 = (TH1D*)MC2->Clone("tempCV2");
  TH1D* sysErr = (TH1D*)MC1->Clone("tempSysCV");
  TH1D* sysErr2 = (TH1D*)MC2->Clone("tempSysCV2");
  TH1D* tmpsys = (TH1D*)MC1->Clone("tempSysRatio");
  TH1D* tmpsys2 = (TH1D*)MC2->Clone("tempSysRatio2");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC->GetBinError(bin)/tmpMC->GetBinContent(bin));
    tmpsys2->SetBinContent(bin,1.0);
    tmpsys2->SetBinError(bin,tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin));
    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
    sysErr->SetBinError(bin,tmpMC->GetBinError(bin));
    sysErr2->SetBinContent(bin,tmpMC2->GetBinContent(bin));
    sysErr2->SetBinError(bin,tmpMC2->GetBinError(bin));
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  tmpMC->SetLineColor(kRed);
  tmpMC->SetLineWidth(3);
  tmpMC->SetLineStyle(1);
  tmpMC2->SetLineColor(kMagenta);
  tmpMC2->SetLineWidth(3);
  tmpMC2->SetLineStyle(2);
  sysErr->SetTitleSize(80);
  sysErr->SetTitleFont(43);
  sysErr->SetTitle(title.c_str());
  sysErr->SetLineColor(kRed-10);
  sysErr->SetLineWidth(1);
  sysErr2->SetLineColor(kMagenta-10);
  sysErr2->SetLineWidth(1);
  sysErr2->SetFillColor(kMagenta-10);
  sysErr2->SetFillStyle(3144);
  sysErr->SetStats(0);
  sysErr->GetYaxis()->SetTitleSize(25);
  sysErr->GetYaxis()->SetTitleFont(43);
  sysErr->GetYaxis()->SetTitle("Events");
  sysErr->GetXaxis()->SetLabelSize(0);
  sysErr->GetXaxis()->SetTitleFont(43);
  sysErr->GetXaxis()->SetTitle(var.c_str());
  sysErr->SetMaximum(1.5*sysErr->GetMaximum());
  sysErr->Draw("E2");
  sysErr2->Draw("E2 same");
  tmpMC->Draw("HIST same");
  tmpMC2->Draw("HIST same");

  
  // lower plot will be in pad
  can.cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  tmpsys->SetTitle("");
  tmpsys->GetYaxis()->SetTitle("uncertainties");
  tmpsys->GetXaxis()->SetTitle(var.c_str());
  tmpsys->SetStats(0);
  tmpsys->SetLineColor(kRed);
  tmpsys->SetLineWidth(3);
  //tmpsys->SetFillColor(kRed-10);
  //tmpsys->SetFillStyle(1001);
  tmpsys2->SetLineColor(kBlue);
  tmpsys2->SetLineWidth(2);
  tmpsys2->SetMarkerStyle(20);
  tmpsys->SetMarkerStyle(21);
  tmpsys2->SetMarkerColor(kBlue);
  tmpsys->SetMarkerColor(kRed);
  //tmpsys2->SetFillColor(kMagenta-10);
  //tmpsys2->SetFillStyle(3144);
  tmpsys->GetYaxis()->SetTitleSize(25);
  tmpsys->GetYaxis()->SetTitleFont(43);
  tmpsys->GetXaxis()->SetTitleSize(25);
  tmpsys->GetXaxis()->SetTitleOffset(3.5);
  tmpsys->GetYaxis()->SetTitleOffset(1.0);
  tmpsys->GetXaxis()->SetLabelSize(0.1);
  tmpsys->GetYaxis()->SetLabelSize(0.08);
  tmpsys->GetYaxis()->SetNdivisions(210, kTRUE);
  tmpsys->GetXaxis()->SetTitleFont(43);
  tmpsys->SetMaximum(2.0);
  tmpsys->SetMinimum(0.0);
  tmpsys->Draw("PE1");
  tmpsys2->Draw("PE1 same");

  const TAxis *axis = tmpsys->GetXaxis();
  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
  double highX = axis->GetBinUpEdge(  axis->GetLast() );

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(kRed);
  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad
 
  
  pad2->Modified(); // so it updates the pad with the new changes
  //pad2->Draw("");

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}

void DrawDataMCAndSyst(TString cName, TH1D* MC1, TH1D* MC2, TH1D* data, std::string var, std::string title){

  TCanvas can(cName,cName,1200,800);

  TH1D* tmpMC1 = (TH1D*)MC1->Clone("tempCV");
  TH1D* tmpMC2 = (TH1D*)MC2->Clone("tempCVnolee");
  TH1D* sysErr1 = (TH1D*)MC1->Clone("tempSysCV");
  TH1D* sysErr2 = (TH1D*)MC2->Clone("tempSysCVnolee");
  TH1D* tmpData = (TH1D*)data->Clone("tempData");

  TH1D* tmpratio1 = (TH1D*)data->Clone("tempRatio1");
  TH1D* tmpMC1_staterr = (TH1D*)MC1->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC1_staterr->GetNbinsX()+1; bin++){ tmpMC1_staterr->SetBinError(bin,sqrt(tmpMC1_staterr->GetBinContent(bin))); }
  tmpratio1->Reset();
  tmpratio1->Divide(tmpData,tmpMC1_staterr,1.0,1.0,"B");
  std::cout << "nbins data, mc = " << tmpData->Integral() << ", " << tmpMC1->Integral() << std::endl;

  TH1D* tmpratio2 = (TH1D*)data->Clone("tempRatio2");
  TH1D* tmpMC2_staterr = (TH1D*)MC2->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC2_staterr->GetNbinsX()+1; bin++){ tmpMC2_staterr->SetBinError(bin,sqrt(tmpMC2_staterr->GetBinContent(bin))); }
  tmpratio2->Reset();
  tmpratio2->Divide(tmpData,tmpMC2_staterr,1.0,1.0,"B");
  std::cout << "nbins data, mc nolee = " << tmpData->Integral() << ", " << tmpMC2->Integral() << std::endl;

  TH1D* tmpsys = (TH1D*)MC2->Clone("tempSysRatio");
  TH1D* tmpsys2 = (TH1D*)MC1->Clone("tempSysRatio1");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin));
    std::cout << "tempSys MC2   " << tmpsys->GetBinError(bin) << std::endl;
    tmpsys2->SetBinContent(bin,1.0);
    tmpsys2->SetBinError(bin,tmpMC1->GetBinError(bin)/tmpMC1->GetBinContent(bin));
    if(bin==1) std::cout << "tempSys MC1   " << tmpMC1->GetBinError(bin) << "/" << tmpMC1->GetBinContent(bin) << " = " << tmpMC1->GetBinError(bin)/tmpMC1->GetBinContent(bin) << std::endl;
    if(bin==1) std::cout << "tempSys MC2   " << tmpMC2->GetBinError(bin) << "/" << tmpMC2->GetBinContent(bin) << " = " << tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin) << std::endl;
    std::cout << "tempSys MC1   " << tmpsys2->GetBinError(bin) << std::endl;
    //tmpratio->SetBinError(bin,0.01);
    sysErr1->SetBinContent(bin,tmpMC1->GetBinContent(bin));
    sysErr1->SetBinError(bin,tmpMC1->GetBinError(bin));
    sysErr2->SetBinContent(bin,tmpMC2->GetBinContent(bin));
    sysErr2->SetBinError(bin,tmpMC2->GetBinError(bin));
    //std::cout << "temp ratio = " << tmpratio->GetBinContent(bin) << std::endl; 
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1, 1.0);
  pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  double max = 1.1;
  double max1 = 1.1;
  double maxMC1 = tmpMC1->GetMaximum();
  double maxMC2 = tmpMC2->GetMaximum();
  double maxData = tmpData->GetMaximum();
  if(maxMC1 > maxData) max1 *= maxMC1;
  if(maxData > maxMC1 ) max1 *= maxData;
  if(maxMC2 > max1) max *= maxMC2;
  else max = max1;

  sysErr1->SetTitleSize(80);
  sysErr1->SetTitleFont(43);
  sysErr1->SetTitle(title.c_str());
  sysErr1->SetStats(0);
  sysErr1->GetYaxis()->SetTitleSize(25);
  sysErr1->GetYaxis()->SetTitleFont(43);
  sysErr1->GetYaxis()->SetTitle("Events");
  sysErr1->GetXaxis()->SetLabelSize(0);

  tmpMC1->SetLineColor(kBlue);
  tmpMC1->SetLineWidth(3);
  tmpMC1->SetLineStyle(1);
  tmpMC2->SetLineColor(kRed);
  tmpMC2->SetLineWidth(3);
  tmpMC2->SetLineStyle(1);
  tmpData->SetLineWidth(2);
  tmpData->SetMarkerStyle(20);
  tmpData->SetLineColor(kBlack);
  if(std::string(cName).find("nowgt") != std::string::npos && std::string(cName).find("nue") != std::string::npos ) tmpData->SetLineColor(kBlack);

  sysErr1->SetLineColor(kBlue-10);
  sysErr1->SetLineWidth(1);
  sysErr1->SetFillColor(kBlue-10);
  sysErr1->SetFillStyle(1001);

  sysErr2->SetLineColor(kRed-10);
  sysErr2->SetLineWidth(1);
  sysErr2->SetFillColor(kRed-10);
  sysErr2->SetFillStyle(1001);

  sysErr1->SetMaximum(1.2*max);
  sysErr1->SetMinimum(0.);

  sysErr1->Draw("E2");
  sysErr2->Draw("E2 same");
  tmpMC1->Draw("hist same"); 
  tmpMC2->Draw("hist same"); 
  tmpData->Draw("PE same"); 
  if(std::string(cName).find("nowgt") != std::string::npos && std::string(cName).find("nue") != std::string::npos ) tmpData->Draw("hist same"); 

  double x1=0.0; 
  double y1=0.0; 
  double x2=0.0; 
  double y2=0.0; 
  if(title == "#nu_{#mu} Selection"){
    x1=0.7;
    x2=0.9;
    y1=0.75;
    y2=0.9;
  }else{
    x1=0.6;
    x2=0.9;
    y1=0.6;
    y2=0.9;
  } 
  TLegend  *legend = new TLegend(x1,y1,x2,y2); // we need different positions for the legend to not 
  // get the plot titles for the legend
  if(title == "#nu_{#mu} Selection"){
    if(std::string(cName).find("nowgt") != std::string::npos ){
      legend->AddEntry(tmpMC1,"#nu_{#mu} CC, post-constraint","l"); 
      legend->AddEntry(tmpMC2,"#nu_{#mu} CC, no Genie wgt","l"); 
      legend->AddEntry(tmpData,"fake data","l");
      //legend->AddEntry(tmpData,"All_UBGenie_Univ_75","l");
    }else{
      legend->AddEntry(tmpMC2,"#nu_{#mu} CC","l"); 
      legend->AddEntry(tmpData,"fake data","l");
      //legend->AddEntry(tmpData,"All_UBGenie_Univ_75","l");
    }
  }else{
    if(std::string(cName).find("nowgt") != std::string::npos ){
      legend->AddEntry(tmpMC1,"#nu_e, Genie Wgt","l"); 
      legend->AddEntry(tmpMC2,"#nu_e, no Genie wgt","l"); 
      legend->AddEntry(tmpData,"#nu_e, post-constraint","l");
    }else{
      legend->AddEntry(tmpMC1,"#nu_e intrinsic + LEE","l"); 
      legend->AddEntry(tmpMC2,"#nu_e intrinsic","l"); 
      legend->AddEntry(tmpData,"fake data","l");
      //legend->AddEntry(tmpData,"All_UBGenie_Univ_75","l");
    }
  } 
  legend->Draw("same");

  // lower plot will be in pad
  can.cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.45);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  tmpsys2->SetTitle("");
  tmpsys2->GetYaxis()->SetTitle("(BNB/MC+EXT)");
  tmpsys2->GetXaxis()->SetTitle(var.c_str());
  tmpsys2->SetStats(0);
  tmpsys2->SetLineColor(kBlue-10);
  tmpsys2->SetLineWidth(1);
  tmpsys2->SetFillColor(kBlue-10);
  tmpsys2->SetFillStyle(1001);
  tmpsys->SetLineColor(kRed-10);
  tmpsys->SetLineWidth(1);
  tmpsys->SetFillColor(kRed-10);
  tmpsys->SetFillStyle(1001);
  tmpsys2->GetYaxis()->SetTitleSize(25);
  tmpsys2->GetYaxis()->SetTitleFont(43);
  tmpsys2->GetXaxis()->SetTitleSize(25);
  tmpsys2->GetXaxis()->SetTitleOffset(3.5);
  tmpsys2->GetYaxis()->SetTitleOffset(1.0);
  tmpsys2->GetXaxis()->SetLabelSize(0.1);
  tmpsys2->GetYaxis()->SetLabelSize(0.08);
  tmpsys2->GetYaxis()->SetNdivisions(210, kTRUE);
  tmpsys2->GetXaxis()->SetTitleFont(43);
  std::cout << "@@@@@@ tmpratio2, tmpratio1 = " << tmpratio2->GetMaximum() << ", " << tmpratio1->GetMaximum() << ", " << tmpratio2->GetMinimum() << ", " << tmpratio1->GetMinimum() << std::endl;
  if(tmpratio2->GetMaximum() >= tmpratio1->GetMaximum() ) tmpsys2->SetMaximum(1.1*tmpratio2->GetMaximum()); 
  if(tmpratio1->GetMaximum() >= tmpratio2->GetMaximum() ) tmpsys2->SetMaximum(1.1*tmpratio1->GetMaximum()); 
  if(tmpratio2->GetMinimum() <= tmpratio1->GetMinimum() ) tmpsys2->SetMinimum(0.8*tmpratio2->GetMinimum()); 
  if(tmpratio1->GetMinimum() <= tmpratio2->GetMinimum() ) tmpsys2->SetMinimum(0.8*tmpratio1->GetMinimum()); 
  //tmpsys2->SetMaximum(3.0);
  if(title == "#nu_{#mu} Selection") tmpsys2->SetMaximum(1.5);
  if(title == "#nu_{#mu} Selection") tmpsys2->SetMinimum(0.5);
  tmpsys2->SetMinimum(0.0);
  tmpsys2->SetMaximum(2.0);
  std::cout << "@@@@@@ tmpSys max, min = " << tmpsys2->GetMaximum() << ", " << tmpsys2->GetMinimum() << std::endl;
  tmpsys2->Draw("E2");
  tmpsys->Draw("E2 same");
  tmpsys2->Draw("E2 same");

  const TAxis *axis = tmpratio1->GetXaxis();
  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
  double highX = axis->GetBinUpEdge(  axis->GetLast() );

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(kRed);
  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad

  tmpratio1->SetMarkerStyle(20);
  tmpratio1->SetLineColor(kBlue);
  tmpratio1->SetLineWidth(2);
  tmpratio1->Draw("PE same");

  tmpratio2->SetMarkerStyle(20);
  tmpratio2->SetLineColor(kRed);
  tmpratio2->SetLineWidth(2);
  tmpratio2->Draw("PE same");

  pad2->Update();
  pad2->Modified(); // so it updates the pad with the new changes
  pad2->Draw("");

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}

void plot_one(TMatrixD matrix, TH1D *h_nue, TH1D *h_numu, std::string tag){

    std::string dir="/uboone/data/users/wospakrk/SBNFitPlots/";
    std::vector<std::string> channel_names = {"nu uBooNE nue intrinsic","nu uBooNE numu BNB"};
    int num_channels = 2;
    int num_bins_total = matrix.GetNrows();
    std::vector<int> num_bins = {h_nue->GetNbinsX(), h_numu->GetNbinsX()}; 
    TH2D h2_full(matrix);
    h2_full.SetName((tag+"_th2d").c_str());
    TCanvas *c_full = new TCanvas((tag+"_canvas").c_str());
    TPad *p_full = (TPad*)c_full->cd();
    c_full->SetFixedAspectRatio();
    h2_full.Draw("colz");
    h2_full.SetTitle(tag.c_str());
    h2_full.GetXaxis()->SetTitle("Global Bin Number");
    h2_full.GetYaxis()->SetTitle(" ");
    h2_full.GetYaxis()->SetLabelSize(0);

    c_full->SetFrameFillColor(kWhite);
    c_full->SetFillColor(kWhite);
    p_full->SetFillColor(kWhite);

    c_full->SetRightMargin(0.150);
    c_full->SetLeftMargin(0.250);
    c_full->SetTopMargin(0.10);
    int use_full =0;

    double percent_left = 0.15;
    double nice_shift = num_bins_total*0.02;

    for(int ic = 0; ic < num_channels; ic++){

        TText * tmd = new TText(-num_bins_total*percent_left*0.15, use_full+nice_shift*0.5, (channel_names[ic]).c_str() );
        tmd->SetTextColor(kBlack);
        tmd->SetTextSize(0.03);
        tmd->SetTextAlign(31);
        tmd->Draw();
       
        TLine *lv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
        TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
        lv->SetLineWidth(3);
        lh->SetLineWidth(3);
        lv->SetLineColor(kRed);
        lh->SetLineColor(kRed);
        use_full+=num_bins.at(ic);
        lv->Draw();
        lh->Draw();

    }
    gStyle->SetOptStat(0);
    c_full->Print((tag+"_collapsed.pdf").c_str(),"pdf");
}

double CalcChiNumu(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix, TH1D *h_nue){

  double tchi=0.0;

  for( int i = 0; i < inv_matrix.GetNcols(); i++ ){
    for( int j = 0; j < inv_matrix.GetNrows(); j++ ){
      if( i >= (h_nue->GetNbinsX()) && j >= (h_nue->GetNbinsX()) ){
        int numubin_i = i-h_nue->GetNbinsX();
        int numubin_j = j-h_nue->GetNbinsX();
        //std::cout << "i, j, numubin_i, numubin_j, h_data_i, h_mc_i, cov_matrix, h_data_j,  = " << i << ", " << j << ", " << numubin_i << ", " << numubin_j << ", " << h_data->GetBinContent(numubin_i+1)-h_mc->GetBinContent(numubin_j+1) << ", " << inv_matrix(i,j) << ", " << h_data->GetBinContent(numubin_j+1)-h_mc->GetBinContent(numubin_i+1) <<std::endl;
        double tchiperbin = (h_data->GetBinContent(numubin_i+1)-h_mc->GetBinContent(numubin_i+1))*inv_matrix(i,j)*(h_data->GetBinContent(numubin_j+1)-h_mc->GetBinContent(numubin_j+1));
        tchi += tchiperbin;
      }
    }
  }
  return tchi;

}

double CalcChiNue(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix){

  double tchi=0.0;

  for( int i = 0; i < inv_matrix.GetNcols(); i++ ){
    for( int j = 0; j < inv_matrix.GetNrows(); j++ ){
      //if(i==j) std::cout << "i,j = " << i << ", " << j << std::endl;
      if( i < (h_data->GetNbinsX()) && j < (h_data->GetNbinsX()) ){
        //std::cout << "i, j, i, j, h_data_i, h_mc_i, cov_matrix, h_data_j,  = " << i << ", " << j << ", " << i << ", " << j << ", " << h_data->GetBinContent(i+1)-h_mc->GetBinContent(j+1) << ", " << inv_matrix(i,j) << ", " << h_data->GetBinContent(j+1)-h_mc->GetBinContent(i+1) <<std::endl;
        double tchiperbin = (h_data->GetBinContent(i+1)-h_mc->GetBinContent(i+1))*inv_matrix(i,j)*(h_data->GetBinContent(j+1)-h_mc->GetBinContent(j+1));
        tchi += tchiperbin;
      }
    }
  }
  return tchi;

}

double CalcChi(std::vector<double> data, std::vector<double> mc, TMatrixD inv_matrix){

  double tchi=0.0;

  for( int i = 0; i < inv_matrix.GetNcols(); i++ ){
    for( int j = 0; j < inv_matrix.GetNrows(); j++ ){
        double tchiperbin = (data[i]-mc[i])*inv_matrix(i,j)*(data[j]-mc[j]);
        tchi += tchiperbin;
    }
  }
  return tchi;

}
