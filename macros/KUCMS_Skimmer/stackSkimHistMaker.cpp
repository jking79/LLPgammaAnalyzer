//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "skimHistMaker.hh"

#include <TStyle.h>
#include <TColorWheel.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCanvas.h>

//#define DEBUG true
#define DEBUG false
#define doEBEEmaps false

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// HistMaker class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

void HistMaker::histMaker( std::string indir, std::string ifilelist, std::string ofilename, std::vector<std::vector<std::string>> deflist, std::vector<float> params ){

    //bool debug = true;
    bool debug = false;

	bkglist = deflist[0];
	siglist = deflist[1];
	datalist = deflist[2];
    bkgleg = deflist[3];
    sigleg = deflist[4];
    dataleg = deflist[5];
	title = deflist[6];
	varsel = deflist[7];

	lumi = params[0];
	maxy = params[1];
    miny = params[2];
    maxr = params[3];

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");

    preCutNPhotons = 0;
    preCut30NPhotons = 0;
    preCut100NPhotons = 0;
    postCutNPhotons = 0;
    postCut30NPhotons = 0;
    postCut100NPhotons = 0;

	std::cout << "Producing Histograms for : " << ofilename << std::endl;
    std::ifstream infile(listdir+ifilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    auto fConfigTree = new TChain(configtreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << ifilelist << " >> " << fInTree << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
        if( str[0] == '#' ) continue;
        if( str[0] == ' ' ) continue;
        auto tfilename = eosdir + indir + str;
        std::cout << "--  adding file: " << tfilename << std::endl;
        fInTree->Add(tfilename.c_str());
        fConfigTree->Add(tfilename.c_str());
		cnt++;
    }//<<>>while (std::getline(infile,str))


	Init(fInTree);
	initHists("");

    std::cout << "Filling Config Map." << std::endl;

    UInt_t          nEvents;
    UInt_t          nSelectedEvents;
    string          *sKey;
    Float_t         sCrossSection;
    Float_t         sGMSBGravMass;
    Float_t         sGMSBChi1Mass;
    Float_t         sMCWgt;
    Int_t           sMCType;

    TBranch        *b_nEvents;   //!
    TBranch        *b_nSelectedEvents;   //!
    TBranch        *b_sKey;   //!
    TBranch        *b_sCrossSection;   //!
    TBranch        *b_sGMSBGravMass;   //!
    TBranch        *b_sGMSBChi1Mass;   //!
    TBranch        *b_sMCWgt;   //!
    TBranch        *b_sMCType;   //!
    TBranch        *b_sumEvtGenWgt;

    sKey = 0;

    fConfigTree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
    fConfigTree->SetBranchAddress("nSelectedEvents", &nSelectedEvents, &b_nSelectedEvents);
    fConfigTree->SetBranchAddress("sKey", &sKey, &b_sKey);
    fConfigTree->SetBranchAddress("sCrossSection", &sCrossSection, &b_sCrossSection);
    fConfigTree->SetBranchAddress("sGMSBGravMass", &sGMSBGravMass, &b_sGMSBGravMass);
    fConfigTree->SetBranchAddress("sGMSBChi1Mass", &sGMSBChi1Mass, &b_sGMSBChi1Mass);
    fConfigTree->SetBranchAddress("sMCWgt", &sMCWgt, &b_sMCWgt);
    fConfigTree->SetBranchAddress("sMCType", &sMCType, &b_sMCType);
    fConfigTree->SetBranchAddress("sumEvtGenWgt", &sumEvtGenWgt, &b_sumEvtGenWgt);

    auto nConfigEntries = fConfigTree->GetEntries();
    std::cout << "Proccessing " << nConfigEntries << " config entries." << std::endl;
    for (Long64_t centry = 0; centry < nConfigEntries; centry++){
		
		auto entry = fConfigTree->LoadTree(centry);

        if(debug) std::cout << " - Getting Branches. " << std::endl;
    	b_nEvents->GetEntry(entry);   //!
    	b_nSelectedEvents->GetEntry(entry);   //!
    	b_sKey->GetEntry(entry);   //!  
    	b_sCrossSection->GetEntry(entry);   //!
    	b_sGMSBGravMass->GetEntry(entry);   //!
    	b_sGMSBChi1Mass->GetEntry(entry);   //!
    	b_sMCWgt->GetEntry(entry);   //!
    	b_sMCType->GetEntry(entry);   //!
        std::string configKey(*sKey);

		if( not configInfo.count(configKey) ){
        	if(debug) std::cout << " - Filling configValues. " << std::endl;
			std::map< std::string, float > configValues;
        	configValues["nEvents"] = nEvents;
        	configValues["nSelectedEvents"] = nSelectedEvents;
        	configValues["sCrossSection"] = sCrossSection;
        	configValues["sGMSBGravMass"] = sGMSBGravMass;
        	configValues["sGMSBChi1Mass"] = sGMSBChi1Mass;
        	configValues["sMCWgt"] = sMCWgt;
        	configValues["sMCType"] = sMCType;
        	if(debug) std::cout << " - Filling configInfo. " << std::endl;
        	configInfo[configKey] = configValues;
		} else {
			auto & configValues = configInfo[configKey];
			configValues["nEvents"] += nEvents;
			configValues["nSelectedEvents"] += nSelectedEvents;
		}//<<>>if( not configInfo.count(configKey) )

	}//<<>>for (Long64_t centry = 0; centry < nConfigEntries; centry++)

    for( auto item : configInfo ){ 
		std::cout << item.first << " ( ";  
		for( auto line : item.second ){ 
			std::cout << line.first <<  " " << line.second << " ";
		}//<<>>for( auto line : item.second )
		std::cout << ")" << std::endl;
	}//<<>>for( auto item : configInfo )

    std::cout << "<<<<<<<< Processing Event Loop <<<<<<<<<<<<<< " << std::endl;
	int loopCounter(100000);
    auto nEntries = fInTree->GetEntries();
    if(debug){ nEntries = 10; loopCounter = 1; }
    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
    for (Long64_t centry = 0; centry < nEntries; centry++){
        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
        if(debug) std::cout << "*****************************************************************************" << std::endl;
        auto entry = fInTree->LoadTree(centry);
		if(debug) std::cout << " - getBranches " << std::endl;
		getBranches(entry);
		if(debug) std::cout << " - eventLoop " << std::endl;
		eventLoop(entry);
    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
   
    if(debug) std::cout << " - Creating output file " << std::endl;
    TFile* fOutFile = new TFile( ofilename.c_str(), "RECREATE" );
    fOutFile->cd();

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

	endJobs();

	for( int it = 0; it < n1dHists; it++ ){ if(sigHist[it]){ sigHist[it]->Write(); delete sigHist[it]; } }
    for( int it = 0; it < n1dHists; it++ ){ if(bkgHist[it]){ bkgHist[it]->Write(); delete bkgHist[it]; } }
    for( int it = 0; it < n1dHists; it++ ){ if(dataHist[it]){ dataHist[it]->Write(); delete dataHist[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }

    fOutFile->Close();
    std::cout << "histMaker : Thats all Folks!!" << std::endl;
}//<<>>void kucmsSkimmer
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------- event loop -------------------------------------------------------------------------
void HistMaker::eventLoop( Long64_t entry ){

    //auto dskey  = *DataSetKey
    //float evtgwt = evtGenWgt;
    //float evtgwt = 1;
    //float scale = lumi;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    float segwt = (configInfo[configKey])["nEvents"];
	float evtgwt = evtGenWgt;
    auto fillwt = lumi * xsec * ( evtgwt / segwt );
    //std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;
    //auto fillwt = (configInfo[configKey])["sCrossSection"] * evtwt;
	// if (s1.find(s2) != std::string::npos) -> true ( s1 contians s2 )

/////////////////    set up varible to plot ////////////////////////////////////////

	//std::string varname = "photonTime_";
	//auto plotv1 = (*selPhoTime)[0];

    std::string varname = title[1]+"_";
	int nPhos = selPhoTime->size();
    //auto plotv3 = (*selPhoR9)[1];
    auto plotTime = (*selPhoTime)[0];
	//auto plotv4 = ( selPhoTime->size() > 1 ) ? (*SCosA)[0] : -10;
    //auto plotv5 = ( selPhoTime->size() == 1 ) ? selCMet : -10;

    auto plot2am1 = (*X2aMass)[0];
    auto plot2ac1 = (*X2aCosA)[0];
	auto plot2bm1 = (*X2bMass)[0];
	auto plot2bc1 = (*X2bCosA)[0];
    auto plot1am1 = (*X1aMass)[0];
    auto plot1ac1 = (*X1aCosA)[0];
    auto plot1bm1 = (*X1bMass)[0];
    auto plot1bc1 = (*X1bCosA)[0];
	auto plot1abm1 = std::sqrt((sq2(plot2am1)+sq2(plot2bm1))/2);

    auto plot2am2 = ( nPhos > 1 ) ? (*X2aMass)[1] : -9999;
    auto plot2ac2 = ( nPhos > 1 ) ? (*X2aCosA)[1] : -9999;
    auto plot2bm2 = ( nPhos > 1 ) ? (*X2bMass)[1] : -9999;
    auto plot2bc2 = ( nPhos > 1 ) ? (*X2bCosA)[1] : -9999;
    auto plot1am2 = ( nPhos > 1 ) ? (*X1aMass)[1] : -9999;
    auto plot1ac2 = ( nPhos > 1 ) ? (*X1aCosA)[1] : -9999;
    auto plot1bm2 = ( nPhos > 1 ) ? (*X1bMass)[1] : -9999;
    auto plot1bc2 = ( nPhos > 1 ) ? (*X1bCosA)[1] : -9999;

	auto plotsm1 = (*SMass)[0];
    auto plotsc1 = (*SCosA)[0];

	auto cmetx = selCMetPx + (*selPhoPt)[0] * std::cos( (*selPhoPhi)[0] );
    auto cmety = selCMetPy + (*selPhoPt)[0] * std::sin( (*selPhoPhi)[0] );
	if( nPhos > 1 ){
		cmetx = cmetx + (*selPhoPt)[1] * std::cos( (*selPhoPhi)[1] );
		cmety = cmety + (*selPhoPt)[1] * std::sin( (*selPhoPhi)[1] );
	}//<<>>if( selPhoTime->size() > 1 )
	auto met = hypo( selCMetPx, selCMetPy ); 
	auto cmet = hypo( cmetx, cmety );

	float plotv1 = -99999.0;
	if( varsel[0] == "x1am11" && nPhos == 1 ) plotv1 = plot1am1; 
    if( varsel[0] == "x2am11" && nPhos == 1 ) plotv1 = plot2am1;
    if( varsel[0] == "unmet1" && nPhos == 1 ) plotv1 = met;
    if( varsel[0] == "crmet1" && nPhos == 1 ) plotv1 = cmet;
    if( varsel[0] == "x2abm11" && nPhos == 1 ) plotv1 = plot1abm1;

	//std::cout << " X2bMass value : " <<  X2bMass->at(0) << std::endl;
    //std::cout << " nRjrPhotons value : " <<  nRjrPhotons->size() << std::endl;

///////////////////////////////////////////////////////////////////////////////////

	if ( configKey.find( "GJets" ) != std::string::npos && selPhoTime->size() == 1 && plotv1 > -99.0 ){

		hist2d[0]->Fill(plot2am1,plotTime);
    	hist2d[1]->Fill(plot2ac1,plotTime);
    	hist2d[2]->Fill(plot2bm1,plotTime);
    	hist2d[3]->Fill(plot2bc1,plotTime);

	}//<<>>if ( configKey.find( "GJets" ) != std::string::npos )


///////////////////////////////////////////////////////////////////////////////////

    int nBkgHists = bkglist.size();
	for( int it = 0; it < nBkgHists; it++ ){
		
		if ( configKey.find( bkglist[it] ) != std::string::npos ){ 

			if( plotv1 > -99.0 ) bkgHist[it]->Fill( plotv1, fillwt ); 
			//bkgHist[it]->Fill( plotv2, fillwt );

		}//<<>>if ( configKey.find( bkglist[it] ) != std::string::npos )

	}//<<>>for( int it = 0; it < nBkg; it++ )

	int nSigHists = siglist.size();
    for( int it = 0; it < nSigHists; it++ ){

        if ( configKey.find( siglist[it] ) != std::string::npos ){

            if( plotv1 > -99.0 ) sigHist[it]->Fill( plotv1, fillwt );

        }//<<>>if ( configKey.find( siglist[it] ) != std::string::npos )

    	if ( selPhoTime->size() == 1 && plotv1 > -99.0 ){

            hist2d[10]->Fill(plot2am1,plotTime);
            hist2d[11]->Fill(plot2ac1,plotTime);
            hist2d[12]->Fill(plot2bm1,plotTime);
            hist2d[13]->Fill(plot2bc1,plotTime);
			hist2d[14]->Fill(plotsm1,plotTime); 
			hist2d[15]->Fill(plotsc1,plotTime);
            hist2d[18]->Fill(plot2am1,plot2bm1);
            hist2d[19]->Fill(plot2ac1,plot2bc1);

		}//<<>>if ( selPhoTime->size() == 1 )

    }//<<>>for( int it = 0; it < nSigHists; it++ )

    int nDataHists = datalist.size();
    for( int it = 0; it < nDataHists; it++ ){

        if ( configKey.find( datalist[it] ) != std::string::npos ){

            if( plotv1 > -99.0 ) dataHist[it]->Fill( plotv1, 1 );

        }//<<>>if ( configKey.find( siglist[it] ) != std::string::npos )

    }//<<>>for( int it = 0; it < nSigHists; it++ )


}//<<>>void HistMaker::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){

    std::cout << " Filling Order Map " << std::endl;

	std::map<float,int> bkgorder;
    std::map<float,int> sigorder;
    std::map<float,int> dataorder;

	int cnt = 0;
	for( auto hist : bkgHist ){ 
		if( hist ){ float size = hist->Integral(); std::cout << " -- Filling : " << size << " cnt: " << cnt << std::endl; bkgorder[size] = cnt;} cnt++;
	}//<<>for( auto hist : bkgHist )
    std::cout << " Bkg Map : ";
    for( auto order : bkgorder ){ std::cout << order.first << " : " << order.second << " : "; }
    std::cout << " Done " << std::endl;

	cnt = bkgorder.size();
    std::cout << " Adding BkgHists - Make TotBkgHist for " << cnt << std::endl;
	cnt--;
	TH1D* totBkgHist = (TH1D*) bkgHist[cnt]->Clone("TOT_BKG");
	for( int iter1 = cnt; iter1 >= 0; iter1-- ){
		//std::cout << "  -- BkgHists iter1 : " << iter1 << std::endl;
		if( iter1 < cnt ) totBkgHist->Add(bkgHist[iter1]);
		for( int iter2 = iter1-1; iter2 >= 0; iter2-- ){
			//std::cout << "  --- Adding iter2 : " << iter2 << std::endl;
			bkgHist[iter1]->Add(bkgHist[iter2]);
		}//<<>>for( int iter2 = iter1-1; iter2 > -1; iter2-- )
	}//<<>>for( int iter1 = cnt-1; iter1 > -1; iter1-- )	
	auto intTotBkg = totBkgHist->Integral();
	std::cout << "Total BackGround : " << intTotBkg << std::endl;
	totBkgHist->Write();

	cnt = 0;
    for( auto hist : sigHist ){ if( hist ){ float size = hist->Integral(); sigorder[size] = cnt; } cnt++; }
	std::cout << " Sig Map : ";	
	for( auto order : sigorder ){ std::cout << order.first << " : " << order.second << " : "; } 
	std::cout << " Done " << std::endl;

    cnt = 0;
    for( auto hist : dataHist ){ if( hist ){ float size = hist->Integral(); dataorder[size] = cnt; } cnt++; }
    std::cout << " Data Map : ";
    for( auto order : dataorder ){ std::cout << order.first << " : " << order.second << " : "; }
    auto intDataHist0 = dataHist[0]->Integral();
    std::cout << " Done " << std::endl;

	std::cout << " Creating Canvas " << std::endl;
	
	float width = bkgHist[0]->GetBinWidth(1);
	std::stringstream convst;
	convst << std::setprecision(4) << width;
	auto widthstr = convst.str();
	std::string yAxisTitle = "Events / ";
	yAxisTitle = yAxisTitle + widthstr + " ";

	auto tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
    tdrStyle->cd();
	tdrStyle->SetOptTitle(0);
	tdrStyle->SetOptStat(0);
  	tdrStyle->SetOptFit(11111111);
    gROOT->ForceStyle();

  	auto can = new TCanvas("can","can",0,0,800,700);
  	can->SetGridx();
  	can->SetGridy();
	//can->SetLogy();
 	can->Draw();
  	can->cd();

    auto c1 = TPad("c1","c1", 0.0, 0.22, 1.0, 1.0);
    c1.SetTopMargin(0.08);
    c1.SetBottomMargin(0.05);
    //if layout['logx'] : c1.SetLogx()
    //if layout['logy'] : 
    c1.SetLogy();
    c1.SetGridx(1);
    c1.SetGridy(1);
    c1.Draw();
    c1.cd();

    std::cout << " Filling Stacked Plot  " << std::endl;

	totBkgHist->SetMaximum(maxy);
	totBkgHist->SetMinimum(miny);
	totBkgHist->GetXaxis()->CenterTitle(true);
	totBkgHist->GetXaxis()->SetTitle("");
    totBkgHist->GetYaxis()->SetTitleFont(42);
	totBkgHist->GetYaxis()->CenterTitle(true);
	totBkgHist->GetYaxis()->SetTitle(yAxisTitle.c_str());
	totBkgHist->SetLineWidth(3.0);
    totBkgHist->SetLineColor(kRed);
    totBkgHist->SetMarkerSize(0);
	totBkgHist->GetXaxis()->SetLabelSize(0);//#20
	totBkgHist->Draw("hist");

	int color = 1;
    std::vector<int> colorStyle{30,38,42,46,49,12};
	//std::vector<int> fillStyle{3305,3395,3345,3354,3325,3352};
	for( auto order = bkgorder.rbegin(); order != bkgorder.rend(); order++ ){

		int indx = order->second;
        bkgHist[indx]->SetLineColor(kBlack);
        bkgHist[indx]->SetLineWidth(1.0);
        bkgHist[indx]->SetFillColor(colorStyle[color]);
        //bkgHist[indx]->SetFillColorAlpha(colorStyle[color],0.85);
        //bkgHist[indx]->SetFillStyle(fillStyle[color]);
        bkgHist[indx]->SetFillStyle(1001);//solid
		bkgHist[indx]->Draw("hist same");
		color++;		

	}//<<>>for( auto order : bkgorder )

    color = 2;
    for( auto order = sigorder.rbegin(); order != sigorder.rend(); order++ ){

        int indx = order->second;
        sigHist[indx]->SetLineColor(color++);
        sigHist[indx]->SetLineWidth(3.0);
        sigHist[indx]->SetMarkerSize(0.);
        sigHist[indx]->SetMarkerColor(kBlack);
        sigHist[indx]->SetLineStyle(7);
        sigHist[indx]->Draw("hist same");

    }//<<>>for( auto order : bkgorder )

    color = 1;
    for( auto order = dataorder.rbegin(); order != dataorder.rend(); order++ ){

        int indx = order->second;
        //sigHist[indx]->SetLineColor(color++);
        //dataHist[indx]->SetLineWidth(1.0);
		dataHist[indx]->SetMarkerStyle(8);
        dataHist[indx]->SetMarkerSize(1.0);
        dataHist[indx]->SetMarkerColor(color++);
        dataHist[indx]->Draw("ep same");

    }//<<>>for( auto order : bkgorder )

	can->cd();
    auto c2 = TPad("c2","c2", 0.0, 0.0, 1.0, 0.24);
    c2.SetTopMargin(0.1);
    c2.SetBottomMargin(0.3);
    c2.SetGridx(1);
    c2.SetGridy(1);
    //if layout['logx'] : c2.SetLogx()
    c2.Draw();
    c2.cd();

    can->Modified();
    can->Update();
	
	TH1D* ratioHist = (TH1D*) totBkgHist->Clone("Ratio_Hist");
	ratioHist->Divide(dataHist[0]);
	
    ratioHist->SetLineColor(1);
    ratioHist->SetMarkerColor(1);
	ratioHist->SetLineWidth(1);
    ratioHist->SetMarkerStyle(8);
    ratioHist->SetMarkerSize(0.8);
    ratioHist->SetMinimum(0.0);  // Define Y ..
    ratioHist->SetMaximum(maxr); // .. range
    //ratioHist->GetXaxis().SetRangeUser(x[0],x[1])
    ratioHist->SetStats(0); // No statistics on lower plot

    ratioHist->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    ratioHist->GetXaxis()->SetLabelFont(43);
    ratioHist->GetXaxis()->SetLabelSize(18);//#20
    ratioHist->GetYaxis()->SetNdivisions(505);

    ratioHist->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    ratioHist->GetYaxis()->SetLabelFont(43);
    ratioHist->GetYaxis()->SetLabelSize(18);//#20
    ratioHist->GetYaxis()->SetNdivisions(505);
    
	ratioHist->GetXaxis()->CenterTitle(true);
    ratioHist->GetXaxis()->SetTitleFont(43);
    ratioHist->GetXaxis()->SetTitleSize(20);
    ratioHist->GetXaxis()->SetTitleOffset(4.0);
    ratioHist->GetXaxis()->SetTitle(title[0].c_str());

    ratioHist->GetYaxis()->CenterTitle(false);
    ratioHist->GetYaxis()->SetTitleFont(43);
    ratioHist->GetYaxis()->SetTitleSize(20);//#32
    ratioHist->GetYaxis()->SetTitleOffset(1.35);
    ratioHist->GetYaxis()->SetTitle("[Data/MC]");

	ratioHist->Draw("ep same");

    can->Modified();
    can->Update();

	//can->cd();
	c1.cd();

  	TLegend* leg = new TLegend(0.71,0.615,0.915,0.915);
  	leg->SetTextFont(132);
  	//leg->SetTextSize(0.045);
  	leg->SetFillColorAlpha(kWhite,0.1);
  	//leg->SetLineColor(kWhite);
  	//leg->SetShadowColor(kWhite);
	for( auto order = bkgorder.rbegin(); order != bkgorder.rend(); order++ ){ leg->AddEntry(bkgHist[order->second],bkgleg[order->second].c_str(),"F"); }	
    for( auto order = sigorder.rbegin(); order != sigorder.rend(); order++ ){ leg->AddEntry(sigHist[order->second],sigleg[order->second].c_str(),"L"); }
    for( auto order = dataorder.rbegin(); order != dataorder.rend(); order++ ){ leg->AddEntry(dataHist[order->second],dataleg[order->second].c_str(),"P"); }
  	//leg->SetLineColor(kWhite);
  	//leg->SetFillColor(kWhite);
  	//leg->SetShadowColor(kWhite);
  	leg->Draw("SAME");

    can->Modified();
    can->Update();

  	TLatex l;
  	l.SetTextFont(132);
  	l.SetNDC();
  	l.SetTextSize(0.05);
  	//l.SetTextFont(132);
  	l.DrawLatex(0.65,0.943,"");
  	//l.SetTextSize(0.05);
  	l.SetTextFont(42);
  	//l.DrawLatex(0.125,0.943,"#bf{#it{CMS}} Internal 13 TeV Simulation   #Chi_{0}->#gamma,G; #Lambda = 250 GeV");
    l.DrawLatex(0.125,0.943,"#bf{#it{CMS}} Internal 13 TeV Simulation   #Chi_{0}->#gamma,G; c#tau = 200 cm");
  	l.SetTextSize(0.05);
  	l.SetTextFont(132);
    std::stringstream lumiconvst;
    lumiconvst << std::setprecision(4) << lumi;
    auto lumistr = lumiconvst.str();
  	string s_lumi = "#scale[0.6]{#int} #it{L dt} = "+lumistr+" fb^{-1}";
  	l.DrawLatex(0.275,0.85,s_lumi.c_str());	

    can->Modified();
    can->Update();

	can->Write();
	can->Print( "stacked_plot_test.png" );
	can->Close();

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht ){

	for( int it = 0; it < n1dHists; it++ ){ bkgHist[it] = NULL; }
    for( int it = 0; it < n1dHists; it++ ){ sigHist[it] = NULL; }
    for( int it = 0; it < n1dHists; it++ ){ dataHist[it] = NULL; }
	for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

	// ---------: 0 - 100

	std::string name("hist_");
 	std::string title("stpl_");

	//int div = 12;
	//float strt = -15.0;
	//float end = 15.0;
	//std::string varname("photonTime_");

//////////////////////////////////////////   set up histograms //////////////////////////////////////////

	std::string varname = varsel[1];
	int div = std::stoi(varsel[2]);
	float strt = std::stof(varsel[3]);
    float end = std::stof(varsel[4]);

	// r9
    //int div = 20;
    //float strt = 0.0;
    //float end = 1.0;
    //std::string varname("photonR9_");
	// time
    //int div = 80;
    //float strt = -5.0;
    //float end = 15.0;
    //std::string varname("phoTime_");
	// Xmass
    //int div = 50;
    //float strt = 0.0;
	//float end = 1000.0;
    //int div = 70;
    //float strt = -0.00001;
    //float end = 0.00006;
    //int div = 100;
    //float strt = 0.0;
    //float end = 10000.0;
    //std::string varname("X2aMass_");
	// XCos
    //int div = 50;
    //float strt = -1.0;
    //float end = 1.0;
    //std::string varname("SCosA_");
	// MET
    //int div = 56;
    //float strt = 100.0;
    //float end = 1500.0;
    //std::string varname("MET_");

///////////////////////////////////////////////////////////////////////////////////////////////////////

    int nBkgHists = bkglist.size();
    for( int it = 0; it < nBkgHists; it++ ){
        auto ttitle = title + varname + bkglist[it];
        auto tname = name + varname + bkglist[it];;
        bkgHist[it] = new TH1D( tname.c_str(), ttitle.c_str(), div, strt, end);
	}//<<>>for( int it = 0; it < nBkgHists; it++ )
	for( int it = 0; it < n1dHists; it++ ){ if(bkgHist[it]) bkgHist[it]->Sumw2();}
    int nSigHists = siglist.size();
    for( int it = 0; it < nSigHists; it++ ){
        auto ttitle = title + varname + siglist[it];
        auto tname = name + varname + siglist[it];; 
        sigHist[it] = new TH1D( tname.c_str(), ttitle.c_str(), div, strt, end);
    }//<<>>for( int it = 0; it < nBkgHists; it++ )
	for( int it = 0; it < n1dHists; it++ ){ if(sigHist[it]) sigHist[it]->Sumw2();}
    int nDataHists = datalist.size();
    for( int it = 0; it < nDataHists; it++ ){
        auto ttitle = title + varname + datalist[it];
        auto tname = name + varname + datalist[it];;
        dataHist[it] = new TH1D( tname.c_str(), ttitle.c_str(), div, strt, end);
    }//<<>>for( int it = 0; it < nDataHists; it++ )
	for( int it = 0; it < n1dHists; it++ ){ if(dataHist[it]) dataHist[it]->Sumw2();}


    //hist1d[0] = new TH1D("phoSigMatchEff", "phoSigMatchEEff", 5, 0, 5); 
    //for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}


	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	hist2d[0] = new TH2D( "rjrX2aMassVtime","Rjr X2a Mass v leadPho Time Gjet;X2a Mass [GeV];time [ns]", 200, 0, 2000, 200, -5, 15 );
    hist2d[1] = new TH2D( "rjrX2aCosAVtime","Rjr X2a CosA v leadPho Time Gjet;X2a CosA;time [ns]", 30, -1.5, 1.5, 200, -5, 15 );
    hist2d[2] = new TH2D( "rjrX2bMassVtime","Rjr X2b Mass v leadPho Time Gjet;X2b Mass [GeV];time [ns]", 200, 0, 2000, 200, -5, 15 );
    hist2d[3] = new TH2D( "rjrX2bCosAVtime","Rjr X2b CosA v leadPho Time Gjet;X2b CosA;time [ns]", 30, -1.5, 1.5, 200, -5, 15 );

    hist2d[10] = new TH2D( "rjrX2aMassVtime1","Rjr X2a Mass v leadPho Time Sig;X2a Mass [GeV];time [ns]", 200, 0, 2000, 200, -5, 15 );
    hist2d[11] = new TH2D( "rjrX2aCosAVtime1","Rjr X2a CosA v leadPho Time Sig;X2a CosA;time [ns]", 30, -1.5, 1.5, 200, -5, 15 );
    hist2d[12] = new TH2D( "rjrX2bMassVtime1","Rjr X2b Mass v leadPho Time Sig;X2b Mass [GeV];time [ns]", 200, 0, 2000, 200, -5, 15 );
    hist2d[13] = new TH2D( "rjrX2bCosAVtime1","Rjr X2b CosA v leadPho Time Sig;X2b CosA;time [ns]", 30, -1.5, 1.5, 200, -5, 15 );

    hist2d[14] = new TH2D( "rjrSMassVtime1","Rjr S Mass v leadPho Time Sig;S Mass [GeV];time [ns]", 200, 0, 2000, 200, -5, 15 );
    hist2d[15] = new TH2D( "rjrSCosAVtime1","Rjr S CosA v leadPho Time Sig;S CosA;time [ns]", 30, -1.5, 1.5, 200, -5, 15 );

    hist2d[18] = new TH2D( "rjrX2aMassVX2b","Rjr X2a Mass v X2b Sig;X2a Mass [GeV];X2a Mass [GeV]", 200, 0, 2000, 200, 0, 2000 );
    hist2d[19] = new TH2D( "rjrX2aCosAVX2b","Rjr X2a CosA v X2b Sig;X2a CosA;X2b CosA", 30, -1.5, 1.5, 30, -1.5, 1.5 );



    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

	//------------------------------------------------------------------------------------

}//<<>>void HistMaker::initHists()

//void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "skim_files/";
				//const std::string listdir = "root://cmseos.fnal.gov//store/user/jaking/llpSkims/";				

				const std::string infilename = "KUCMS_Stack_Skim_List.txt";

                const std::string outfilename = "KUCMS_Skim_StackHists.root"; //9

				//auto htitle = "KUCMS_This_BaseHists_";

				std::vector<std::string> lbkgs{"GJets","QCD"};
				std::vector<std::string> lsigs{"Ct200_GMSBL150","Ct200_GMSBL250","Ct200_GMSBL350"};
                std::vector<std::string> ldata{"JetHT"};

                std::vector<std::string> bkgleg{" GJets"," QCD"};
                //std::vector<std::string> sigleg{" c#tau 0.1", " c#tau 400", " c#tau 1200" };
                std::vector<std::string> sigleg{" g/#chi 1207/212", " g/#chi 1915/358"," g/#chi 2599/503"};
				std::vector<std::string> dataleg{" Data"};

                //std::vector<std::string> title{"Lead Photon R9","leadPhoR9"};
                //std::vector<std::string> title{"SubLead Photon R9","subleadPhoR9"};
                //std::vector<std::string> title{"Lead Photon Time","leadPhoTime"};
                //std::vector<std::string> title{"SubLead Photon Time","subleadPhoTime"};
				//std::vector<std::string> title{"X1a Mass Lead Photon 1pho","pho1Smass"};
                //std::vector<std::string> title{"S CosA 2 Photon","pho2Scosa"};
                //std::vector<std::string> title{"Rjr Crtd Met 1 Photon","crtdmet"};

                //std::vector<std::string> title{"UnCor Met Lead Event 1pho","pho1UnMet"};
                //std::vector<std::string> varsel{"unmet1","UnMet_","50","0.0","1000"};
                //std::vector<std::string> title{"Cor Met Lead Event 1pho","pho1CrMet"};
                //std::vector<std::string> varsel{"crmet1","CrMet_","50","0.0","1000"};              
  
                //std::vector<std::string> title{"X1a Mass Lead Photon 1pho","pho11x1aMass"};
				//std::vector<std::string> varsel{"x1am11","X1a11pMass_","70","-0.00001","0.00006"};
                //std::vector<std::string> title{"X2a Mass Lead Photon 1pho","pho11x2aMass"};
                //std::vector<std::string> varsel{"x2am11","X2a11pMass_","60","0.0","1200"};
                std::vector<std::string> title{"X2ab Mass Lead Photon 1pho","pho11x2abMass"};
                std::vector<std::string> varsel{"x2abm11","X2a11pMass_","60","0.0","1200"};


				std::vector<std::vector<std::string>> deflist{lbkgs,lsigs,ldata,bkgleg,sigleg,dataleg,title,varsel};

				//float lumi = 0.248017205; // d
                //float lumi = 1.432791808; // c
                float lumi = 1.6808; // d+c
				float maxy = 1e6;
				float miny = 1e-6;
				float maxr = 5.0;

				std::vector<float> params{lumi,maxy,miny,maxr};

                HistMaker base;
                base.histMaker( listdir, infilename, outfilename, deflist, params );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

//
//#    legend = TLegend(0.25,0.20,0.52,0.525); # bottom left
//#    legend = TLegend(0.4,0.205,0.6,0.525);   # bottom middle
//#    legend = TLegend(0.4,0.60,0.6,0.90);   # top middle
//#    legend = TLegend(0.645,0.50,0.825,0.9);   # top mid right
//#    legend = TLegend(0.605,0.50,0.945,0.9);   # top right very wide
//#    legend = TLegend(0.705,0.50,0.945,0.9);   # top right wide 
//#    legend = TLegend(0.745,0.50,0.925,0.9);   # top right
//#    legend = TLegend(0.745,0.40,0.925,0.9);   # top right tall
//#    legend = TLegend(0.650,0.375,0.925,0.875);   # top mid right wide
//#    legend = TLegend(0.62,0.60,0.8,0.9);   # top right
//#    legend = TLegend(0.65,0.60,0.9,0.90);   # top right large
//
