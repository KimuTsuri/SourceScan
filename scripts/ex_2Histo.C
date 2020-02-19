{
Int_t n = 4;
Double_t pt[4] = {0.32, 0.74, 1.2, 1.8};
Double_t v2[4] = {0.005, 0.0134, 0.0316, 0.0738};
Double_t errRel[4] = {285, 86, 49, 33};
Double_t zero[4] = {0., 0., 0., 0.}; 

TCanvas *c1 = new TCanvas("c1","transparent pad",200,10,700,500);
TPad *pad1 = new TPad("pad1","",0,0,1,1);
TPad *pad2 = new TPad("pad2","",0,0,1,1);
pad2->SetFillStyle(4000);

pad1->SetLeftMargin(0.2); // ======= THIS IS WHAT I HAVE ADDED
pad1->SetBottomMargin(0.2);
pad1->SetRightMargin(0.2);

pad1->Draw();
pad1->cd();

TH1F *hr = c1->DrawFrame(0., 0., 2.5, 0.1);

hr->SetXTitle("p_{T} [GeV / c]");
hr->GetXaxis()->CenterTitle();
hr->GetXaxis()->SetTitleSize(0.06);
hr->GetXaxis()->SetTitleOffset(1.2);
hr->GetXaxis()->SetLabelSize(0.06);
hr->GetXaxis()->SetLabelOffset(0.02);

hr->SetYTitle("v_{2}");
hr->GetYaxis()->CenterTitle();
hr->GetYaxis()->SetTitleSize(0.06);
hr->GetYaxis()->SetTitleOffset(1.3);
hr->GetYaxis()->SetLabelSize(0.06);
hr->GetYaxis()->SetLabelOffset(0.02);

gr1 = new TGraphErrors(n, pt, v2, zero, zero); // ======== I SIMPLIFIED HERE
gr1->SetMarkerColor(kBlack);
gr1->SetMarkerStyle(20);
gr1->SetMarkerSize(2.);
gr1->SetLineWidth(4);
gr1->Draw("P");

pad1->Update(); //this will force the generation of the "stats" box === I DON'T NEED STATS BOX
c1->cd();

//compute the pad range with suitable margins
Double_t ymin = 0;
Double_t ymax = 60;
Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
Double_t xmin = 0;
Double_t xmax = 2.5;
Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);

pad2->SetLeftMargin(0.2);
pad2->SetBottomMargin(0.2);
pad2->SetRightMargin(0.2);

pad2->Draw();
pad2->cd();

gr2 = new TGraphErrors(n,pt, errRel, zero, zero); // ======= I SIMPLIFIED HERE
gr2->SetMarkerColor(kRed);
gr2->SetMarkerStyle(20);
gr2->SetMarkerSize(2.);
gr2->SetLineWidth(4);
gr2->Draw("P"); //LP

// draw axis on the right side of the pad
TGaxis *axis = new TGaxis(xmax,ymin+3.8,xmax,ymax-3.8,ymin,ymax,50510,"+L");
axis->SetLabelColor(kRed);
axis->SetTitle("#delta v_{2} / v_{2} [%]");
axis->CenterTitle();
axis->SetTitleSize(0.06);
axis->SetTitleColor(kRed);
axis->SetTitleOffset(1.2);
axis->SetLabelSize(0.06);
axis->SetLabelOffset(0.02);
axis->Draw();
}
