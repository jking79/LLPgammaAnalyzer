head	1.1;
access;
symbols;
locks; strict;
comment	@ * @;


1.1
date	2022.09.30.18.05.33;	author jaking;	state Exp;
branches;
next	;


desc
@@


1.1
log
@Initial revision
@
text
@
        cljBc3dEx.push_back(cljBCEigen3D[0]);
        cljBc3dEy.push_back(cljBCEigen3D[1]);
        cljBc3dEz.push_back(cljBCEigen3D[2]);
        cljBc3dEv.push_back(cljBCEigen3D[3]);
        cljBc3dEslope.push_back(cljBCEigen3D[4]);
        cljBc3dEchisp.push_back(cljBCEigen3D[5]);

        cljBc2dEx.push_back(cljBCEigen2D[0]);
        cljBc2dEy.push_back(cljBCEigen2D[1]);
        cljBc2dEv.push_back(cljBCEigen2D[2]);
        cljBc2dEslope.push_back(cljBCEigen2D[3]);
        cljBc2dEchisp.push_back(cljBCEigen2D[4]);

        hist1d[146]->Fill( cljTimeStats[6] );//c mean
        hist1d[147]->Fill( cljSeedTOFTime );//lead time 
        hist1d[148]->Fill( cljTimeStats[6] - cljSeedTOFTime );//diff

        cljSeedTOFTime.push_back(seedTOFTime);
        cljCMeanTime.push_back(timeStats[6]);

		rhPosX.push_back(rhX);
		rhPosY.push_back(rhY);
		rhPosZ.push_back(rhZ);
		rhTOF.push_back(tof);
		rhPosEta.push_back(recHitPos.eta());
    	rhPosPhi.push_back(recHitPos.phi());
		rhTime.push_back(recHit.time());
		rhTimeErr.push_back(recHit.timeError());
		rhSubdet.push_back( idinfo.ecal == ECAL::EB ? 0 : idinfo.ecal == ECAL::EP ? 1 : 2 );
		rhXtalI1.push_back(idinfo.i1);
    	rhXtalI2.push_back(idinfo.i2);
		rhisOOT.push_back(recHit.checkFlag(EcalRecHit::kOutOfTime));
		rhEnergy.push_back(recHit.energy());
		energyError()
		rhisGS6.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain6));
		rhisGS1.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain1));
		 -- see line 2988 dispho code for below
    nRecHits = nRecHitCnt;

	{ fillTH1(recHit.time(),hist1d[131]); fillTH1(recHit.energy(),hist1d[132]); hist2d[122]->Fill(recHit.time(),recHit.energy()); 
	{ fillTH1(recHit.time(),hist1d[133]); fillTH1(recHit.energy(),hist1d[134]); hist2d[123]->Fill(recHit.time(),recHit.energy()); 
		fillTH1(recHit.checkFlag(EcalRecHit::kOutOfTime),hist1d[136]); }


        phoPt.push_back(photon.pt());
        phoEnergy.push_back(photon.energy());
        phoPhi.push_back(photon.phi());
        phoEta.push_back(photon.eta());
        phoPx.push_back(photon.px());
        phoPy.push_back(photon.py());
        phoPz.push_back(photon.pz());

		iGedPhos++;

        if( phoSCEigen3D[0] != -999 ){
            auto epanlge = getAngle( phoSCEigen3D[0], phoSCEigen3D[1] );
            //hist1d[63]->Fill(epanlge);//etaphi angle
            auto ephypo3D = hypo( phoSCEigen3D[0], phoSCEigen3D[1] );
            auto etanlge = getAngle( ephypo3D, phoSCEigen3D[2] );
            //hist1d[64]->Fill(etanlge);//etatim angle
            //hist2d[82]->Fill( phoSCEigen3D[0], phoSCEigen3D[1] );
            //hist2d[83]->Fill( phoSCEigen3D[0], phoSCEigen3D[2] );
            //hist1d[86]->Fill(phoSCEigen3D[3]);
        }//<<>>if( phoSCEigen3D[0] != -999 )
        if( phoSCEigen2D[0] != -999 ){
            auto sphanlge = getAngle( phoSCEigen2D[0], phoSCEigen2D[1] );
            //hist1d[65]->Fill(sphanlge);//eliptical angle
            //hist2d[81]->Fill( phoSCEigen2D[0], phoSCEigen2D[1] );
            //hist1d[81]->Fill(phoSCEigen2D[2]);
            //hist2d[88]->Fill(sphanlge, phoSCEigen2D[2]);
        }//<<>>if( phoSCEigen2D[0] != -999 )

        phoSc3dEx.push_back(phoSCEigen3D[0]);
    	phoSc3dEy.push_back(phoSCEigen3D[1]);
    	phoSc3dEz.push_back(phoSCEigen3D[2]);
        phoSc3dEv.push_back(phoSCEigen3D[3]);
        phoSc3dEslope.push_back(phoSCEigen3D[4]);
        phoSc3dEchisp.push_back(phoSCEigen3D[5]);

        phoSc2dEx.push_back(phoSCEigen2D[0]);
        phoSc2dEy.push_back(phoSCEigen2D[1]);
        phoSc2dEv.push_back(phoSCEigen2D[2]);
        phoSc2dEslope.push_back(phoSCEigen2D[3]);
        phoSc2dEchisp.push_back(phoSCEigen2D[4]);

        hist1d[146]->Fill( phoTimeStats[6] );//c mean
        hist1d[147]->Fill( phoSeedTOFTime );//lead time 
        hist1d[148]->Fill( phoTimeStats[6] - phoSeedTOFTime );//diff

		phoSeedTOFTime.push_back(seedTOFTime);
        phoCMeanTime.push_back(timeStats[6]);


        ootPhoPt.push_back(ootphoton.pt());
        ootPhoEnergy.push_back(ootphoton.energy());
        ootPhoPhi.push_back(ootphoton.phi());
        ootPhoEta.push_back(ootphoton.eta());
        ootPhoPx.push_back(ootphoton.px());
        ootPhoPy.push_back(ootphoton.py());
        ootPhoPz.push_back(ootphoton.pz());

        iOOTPhos++;


        hist1d[149]->Fill( ootPhoTimeStats[6] );//c mean
        hist1d[150]->Fill( ootPhoSeedTOFTime );//lead time 
        hist1d[151]->Fill( ootPhoTimeStats[6] - ootPhoSeedTOFTime );//diff

        ootPhoSc3dEx.push_back(ootPhoSCEigen3D[0]);
        ootPhoSc3dEy.push_back(ootPhoSCEigen3D[1]);
        ootPhoSc3dEz.push_back(ootPhoSCEigen3D[2]);
        ootPhoSc3dEv.push_back(ootPhoSCEigen3D[3]);
        ootPhoSc3dEslope.push_back(ootPhoSCEigen3D[4]);
        ootPhoSc3dEchisp.push_back(ootPhoSCEigen3D[5]);

        ootPhoSc2dEx.push_back(ootPhoSCEigen2D[0]);
        ootPhoSc2dEy.push_back(ootPhoSCEigen2D[1]);
        ootPhoSc2dEv.push_back(ootPhoSCEigen2D[2]);
        ootPhoSc2dEslope.push_back(ootPhoSCEigen2D[3]);
        ootPhoSc2dEchisp.push_back(ootPhoSCEigen2D[4]);

        ootPhoSeedTOFTime.push_back(seedTOFTime);
        ootPhoCMeanTime.push_back(timeStats[6]);


        elePt.push_back(electron.pt());
        eleEnergy.push_back(electron.energy());
        elePhi.push_back(electron.phi());
        eleEta.push_back(electron.eta());
        elePx.push_back(electron.px());
        elePy.push_back(electron.py());
        elePz.push_back(electron.pz());

        iElectros++;


        hist1d[152]->Fill( eleTimeStats[6] );//c mean
        hist1d[153]->Fill( eleSeedTOFTime );//lead time 
        hist1d[154]->Fill( eleTimeStats[6] - eleSeedTOFTime );//diff

        eleSc3dEx.push_back(eleSCEigen3D[0]);
        eleSc3dEy.push_back(eleSCEigen3D[1]);
        eleSc3dEz.push_back(eleSCEigen3D[2]);
        eleSc3dEv.push_back(eleSCEigen3D[3]);
        eleSc3dEslope.push_back(eleSCEigen3D[4]);
        eleSc3dEchisp.push_back(eleSCEigen3D[5]);

        eleSc2dEx.push_back(eleSCEigen2D[0]);
        eleSc2dEy.push_back(eleSCEigen2D[1]);
        eleSc2dEv.push_back(eleSCEigen2D[2]);
        eleSc2dEslope.push_back(eleSCEigen2D[3]);
        eleSc2dEchisp.push_back(eleSCEigen2D[4]);

        eleSeedTOFTime.push_back(seedTOFTime);
        eleCMeanTime.push_back(timeStats[6]);

	   	jetHt += jet.pt();

	   	jetE.push_back(jet.energy());
	   	jetPt.push_back(jet.pt());
	   	jetPhi.push_back(jet.phi());
	   	jetEta.push_back(jet.eta());
		jetEtaetaMmt.push_back(jet.etaetaMoment());
		jetPhiphiMnt.push_back(jet.phiphiMoment());
		jetEtaphiMnt.push_back(jet.etaphiMoment());
		jetMaxD.push_back(jet.maxDistance());
		jetConPtDis.push_back(jet.constituentPtDistribution());
		jetConEtaPhiSprd.push_back(jet.constituentEtaPhiSpread());
		jetArea.push_back(jet.jetArea());
		jetNCarry.push_back(jet.nCarrying(0.1));
		jetNConst.push_back(jet.nConstituents());

	   	jetID.push_back(jet.userInt("jetID"));
		jetID.push_back(jetid);
	   	jetNHF.push_back(jet.neutralHadronEnergyFraction());
	   	jetNEMF.push_back(jet.neutralEmEnergyFraction());
	   	jetCHF.push_back(jet.chargedHadronEnergyFraction());
	   	jetCEMF.push_back(jet.chargedEmEnergyFraction());
	   	jetMUF.push_back(jet.muonEnergyFraction());
	   	jetNHM.push_back(jet.neutralMultiplicity());
		jetCHM.push_back(jet.chargedMultiplicity());
		jetCharge.push_back(jet.jetCharge());

		jetPHE.push_back(jet.photonEnergy());
		jetPHEF.push_back(jet.photonEnergyFraction()); 
		jetELE.push_back(jet.electronEnergy());
		jetELEF.push_back(jet.electronEnergyFraction());
		jetMUE.push_back(jet.muonEnergy());
		jetPHM.push_back(jet.photonMultiplicity());
		jetELM.push_back(jet.electronMultiplicity());

   <<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
   
	   	if( DEBUG ) std::cout << "Fill jet pt/phi/eta Histograms" << std::endl;
 		
      	const auto jetepafrac 	= jet.photonEnergyFraction() + jet.electronEnergyFraction();
      	const auto jetepe 		= jet.photonEnergy() + jet.electronEnergy();
		const auto jeteme 		= jet.chargedEmEnergy() + jet.neutralEmEnergy();
      	const auto jetemfrac 	= jeteme/jet.energy();
      	const auto jetepfrac 	= jetepe/jet.energy();

		jetSumEPFrac.push_back(jetepafrac);
		jetEPEnergy.push_back(jetepe);
		jetEMEnergy.push_back(jeteme);
		jetEMEnrFrac.push_back(jetemfrac);
		jetEPEnrFrac.push_back(jetepfrac);

		hist2d[61]->Fill(jetepafrac,jetepfrac);
      	hist2d[62]->Fill(jetepfrac,jetemfrac);

	   	fillTH1(jet.pt(),hist1d[12]);//hist1d[12]->Fill(jet.pt());
	   	fillTH1(jet.phi(),hist1d[13]);//hist1d[13]->Fill(jet.phi());
	   	fillTH1(jet.eta(),hist1d[14]);//hist1d[14]->Fill(jet.eta());

		sJetDrRHEnergy.push_back(sumdrrhe);
		jetDrEMF.push_back(dremf);
        jetDrRhCnt.push_back(rhCount);

		   	const auto leadJetRhId = leadJetRh.detid();
		   	const auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
		   	auto sc_eta = leadJetRhIdPos.eta();
		   	auto sc_phi = leadJetRhIdPos.phi();
		   	auto sc_enr = leadJetRh.energy();
			  make jettime varible
	      	for( auto t : tofTimes ) hist1d[5]->Fill(t);
			auto jetTimeStats = getTimeDistStats( tofTimes, jetDrRhGroup ); 
			auto jmutime = jetTimeStats[0];
		   	auto jterr = jetTimeStats[1];
		   	auto jtrms = jetTimeStats[4];
			auto jmedtime = jetTimeStats[2];
	      	auto mederr = jetTimeStats[3];
			auto jcmutime = jetTimeStats[6];
	      	auto jcmedtime = jetTimeStats[10];

			jetDrTime = jcmutime;
		   	njetRecHits.push_back(rhCount);
		   	jetMuTime.push_back(jmutime);
		   	jetTimeError.push_back(jterr);
		   	jetTimeRMS.push_back(jtrms);
		   	jetMedTime.push_back(jmedtime);
	      	jetCMuTime.push_back(jcmutime);	
	      	jetCMedTime.push_back(jcmedtime);

			jetDrLeadEta.push_back(sc_eta);
            jetDrLeadPhi.push_back(sc_phi);
            jetDrLeadEnr.push_back(sc_enr);


		   	fillTH1(jmutime,hist1d[29]);//hist1d[29]->Fill(jmutime);
		   	fillTH1(rhCount,hist1d[1]);//hist1d[1]->Fill(rhCount);
		   	fillTH1(jterr,hist1d[2]);//hist1d[2]->Fill(jterr);
		   	fillTH1(jtrms,hist1d[3]);//hist1d[3]->Fill(jtrms);
		   	fillTH1(jmedtime,hist1d[4]);//hist1d[4]->Fill(jmedtime);
	      	fillTH1(jcmutime,hist1d[6]);//hist1d[6]->Fill(jcmutime);
	      	fillTH1(jcmedtime,hist1d[7]);//hist1d[7]->Fill(jcmedtime);

		   	hist2d[1]->Fill(jmutime,jet.pt());
		   	//hist2d[2]->Fill(jmutime,jet.userInt("jetID"));
            hist2d[2]->Fill(jmutime,jetid);
		   	hist2d[3]->Fill(jmutime,jet.neutralHadronEnergyFraction());
		   	hist2d[4]->Fill(jmutime,jet.chargedHadronEnergyFraction());
		   	hist2d[5]->Fill(jmutime,jet.neutralEmEnergyFraction());
		   	hist2d[6]->Fill(jmutime,jet.chargedEmEnergyFraction());
		   	hist2d[7]->Fill(jmutime,jet.muonEnergyFraction());
		   	hist2d[8]->Fill(jmutime,jet.neutralMultiplicity());
		   	hist2d[9]->Fill(jmutime,jet.chargedMultiplicity());
	
		   	hist2d[10]->Fill(jmutime,jmedtime);
		   	hist2d[24]->Fill(jmutime,rhCount);
		   	hist2d[25]->Fill(jmedtime,rhCount);
		   	hist2d[11]->Fill(jmutime,jtrms);
		   	hist2d[12]->Fill(jmutime,jterr);
		
		   	hist2d[32]->Fill(jmutime,sc_eta);
		   	hist2d[33]->Fill(jmutime,sc_phi);
		   	hist2d[34]->Fill(jmutime,sc_enr);
		   	hist2d[35]->Fill(jmedtime,sc_eta);
		   	hist2d[36]->Fill(jmedtime,sc_phi);
		   	hist2d[37]->Fill(jmedtime,sc_enr);
	
	
		   	hist2d[13]->Fill(jmedtime,jet.pt());
		   	//hist2d[14]->Fill(jmedtime,jet.userInt("jetID"));
            hist2d[14]->Fill(jmedtime,jetid);
		   	hist2d[15]->Fill(jmedtime,jet.neutralHadronEnergyFraction());
		   	hist2d[16]->Fill(jmedtime,jet.chargedHadronEnergyFraction());
		   	hist2d[17]->Fill(jmedtime,jet.neutralEmEnergyFraction());
		   	hist2d[18]->Fill(jmedtime,jet.chargedEmEnergyFraction());
		   	hist2d[19]->Fill(jmedtime,jet.muonEnergyFraction());
		   	hist2d[20]->Fill(jmedtime,jet.neutralMultiplicity());
		   	hist2d[21]->Fill(jmedtime,jet.chargedMultiplicity());


            	auto gjeta = genJet.eta();
            	auto gjphi = genJet.phi();
                auto jtgjdr = reco::deltaR2(gjeta, gjphi, jet.eta(), jet.phi() );
                if( jtgjdr <= 0.1 && rhCount > 0 ){

            		std::cout << " - genJet mothers : " << nMother << " daughters : " << nDaughter << " sources : " << nSources << std::endl;
            		if( DEBUG ) std::cout << " - genJet srcs : " << nSources << " PV (" << vtxX << "," << vtxY << "," << vtxZ << ")" << std::endl;
            		auto kids = genJet.daughterPtrVector();
            		std::cout << bigKidChase( kids, vtxX ) << std::endl;
            		kidChase( kids, vtxX, vtxY, vtxZ );
					auto leadJetRh = getLeadRh( jetDrRhGroup );
            		auto leadJetRhId = leadJetRh.detid();
            		auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
					auto cx = leadJetRhIdPos.x();
					auto cy = leadJetRhIdPos.y();
					auto cz = leadJetRhIdPos.z();
					auto tofcor = hypo( cx, cy, cz )/SOL;
					if( DEBUG ) kidChase( kids, vtxX, vtxY, vtxZ );
            		auto genKidInfo = kidTOFChain( kids, cx, cy, cz );
					if( DEBUG ) std::cout << " - genJet GenTime noTOF : " << genKidInfo[0] << " rhPos: " << cx << "," << cy << "," << cz << std::endl;
					genEta = genJet.eta();
					if( genKidInfo[0] > 25.0 ) genTime = -28.0;
					else if( genKidInfo[0] > -25.0 ) genTime = genKidInfo[0]-tofcor;
					else genTime = -27.0;
					genImpactAngle = genKidInfo[1];
					if( DEBUG ) std::cout << " - genJet GenTime : " << genTime << " Angle: " << genImpactAngle << std::endl;
					genPt = genJet.pt();
                	genEnergy = genJet.energy();
					genEMFrac = (genJet.chargedEmEnergy() + genJet.neutralEmEnergy())/genEnergy;
					genDrMatch = jtgjdr; std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), genJet.eta(), genJet.phi()));
					genTimeVar = genKidInfo[2];
                	genNextBX = genKidInfo[3];
                	genTimeLLP = genKidInfo[4];
                	genLLPPurity = genKidInfo[5];
					genNKids = genKidInfo[6];
					genTOF = tofcor;
					if( DEBUG ) std::cout << " -- Energy : " << genEnergy << " Pt : " << genPt << " EMfrac : " << genEMFrac << std::endl;
                	hist2d[109]->Fill(genTime[0],tofcor);
                	hist1d[110]->Fill(genTime[0]);
                	hist1d[111]->Fill(tofcor);
					
            jetGenTime.push_back(genTime);
            jetGenPt.push_back(genPt);
            jetGenEta.push_back(genEta);
            jetGenEnergy.push_back(genEnergy);
            jetGenEMFrac.push_back(genEMFrac);
            jetGenDrMatch.push_back(genDrMatch);
            jetGenTimeVar.push_back(genTimeVar);
            jetGenTimeLLP.push_back(genTimeLLP);
            jetGenLLPPurity.push_back(genLLPPurity);
            jetGenNextBX.push_back(genNextBX);
            jetGenNKids.push_back(genNKids);
            jetGenTOF.push_back(genTOF);


	   	hist1d[43]->Fill( nMatched ); // # of SC matched to a jet
	    hist2d[40]->Fill( jet.energy(), sum_sce );
	    hist2d[41]->Fill( jet.energy(), sum_phe );
	    hist2d[42]->Fill( sum_phe, sum_sce );

		nJetScMatch.push_back(nMatched);
		sJetScEnergy.push_back(sum_sce);
		sJetScPhEnergy.push_back(sum_phe);

			jetImpactAngle.push_back(impangle);
			jetSc3dEx.push_back(jetSCEigen3D[0]);
            jetSc3dEy.push_back(jetSCEigen3D[1]);
            jetSc3dEz.push_back(jetSCEigen3D[2]);
            jetSc3dEv.push_back(jetSCEigen3D[3]);
            jetSc3dEslope.push_back(jetSCEigen3D[4]);
            jetSc3dEchisp.push_back(jetSCEigen3D[5]);	

            jetSc2dEx.push_back(jetSCEigen2D[0]);
            jetSc2dEy.push_back(jetSCEigen2D[1]);
            jetSc2dEv.push_back(jetSCEigen2D[2]);
            jetSc2dEslope.push_back(jetSCEigen2D[3]);
            jetSc2dEchisp.push_back(jetSCEigen2D[4]);
            jetSc2dEslope2.push_back(jetSCEigen2D[5]);
            jetSc2dEchisp2.push_back(jetSCEigen2D[6]);
            jetSc2dErangle.push_back(jetSCEigen2D[7]);
            jetSc2dEnxsum.push_back(jetSCEigen2D[8]);

			if( jetSCEigen3D[0] != -999 ){
					auto epanlge = getAngle( jetSCEigen3D[0], jetSCEigen3D[1] );
            		hist1d[63]->Fill(epanlge);//etaphi angle
					auto ephypo3D = hypo( jetSCEigen3D[0], jetSCEigen3D[1] );
            		auto etanlge = getAngle( ephypo3D, jetSCEigen3D[2] );
            		fillTH1(etanlge,hist1d[64]);//hist1d[64]->Fill(etanlge);//etatim angle
            		hist2d[82]->Fill( jetSCEigen3D[0], jetSCEigen3D[1] );
            		hist2d[83]->Fill( jetSCEigen3D[0], jetSCEigen3D[2] );
					hist1d[86]->Fill(jetSCEigen3D[3]);
                	if( jetSCEigen3D[5] > 0.95 && jetSCEigen3D[3] < 0.9 && jetSCEigen3D[3] > 0.7 ){
						hist2d[96]->Fill( impangle, jetSCEigen3D[4] );
                        hist2d[126]->Fill( jet.eta(), jetSCEigen3D[4] );
                        hist2d[127]->Fill( jet.phi(), jetSCEigen3D[4] );
					}//<<>>if( jetSCEigen3D[5] > 0.95 ):
				}//<<>>if( jetSCEigen3D[0] != -999 )
				if( jetSCEigen2D[0] != -999 ){
					auto sphanlge = getAngle( jetSCEigen2D[0], jetSCEigen2D[1] );
		    		hist1d[65]->Fill(sphanlge);//eliptical angle
            		hist2d[81]->Fill( jetSCEigen2D[0], jetSCEigen2D[1] );
            		hist1d[81]->Fill(jetSCEigen2D[2]);
					hist2d[88]->Fill( sphanlge, jetSCEigen2D[2] );
					if( jetSCEigen2D[4] > 0.95 && jetSCEigen2D[2] < 0.9 && jetSCEigen2D[2] > 0.7 ){ 
                    	if( jet.eta() > 1.0 ) hist1d[117]->Fill( jetSCEigen2D[3] );
                    	if( jet.eta() < 0.5 ) hist1d[118]->Fill( jetSCEigen2D[3] );
						hist2d[89]->Fill( jet.eta(), jetSCEigen2D[3] );
                    	hist2d[121]->Fill( jet.phi(), jetSCEigen2D[3] );
                		hist2d[90]->Fill( impangle, jetSCEigen2D[3] );
					}//<<>>if( jetSCEigen2D[4] < 0.1 )
					//hist2d[91]->Fill( jetSCEigen2D[3], jetSCEigen2D[4]);
            	}//<<>>if( jetSCEigen2D[0] != -999 )
			}//if( true/false ) on/off switch for jetSC or PhotonSC

    <<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
    	<<<<if( jetScRhGroup.size() >= minRHcnt && scemf > minEmf )

			jetScRhCnt.push_back(jetScRhGroup.size());

			if( jetSCTimeStats[6] > -28.9 ) nGoodScJets++;
            std::cout << " - fill hists " << std::endl;
            for( auto t : jetSCtofTimes ) hist1d[47]->Fill(t);
            if( rhCount >= minRHcnt ) 
            hist2d[39]->Fill( rhCount, jetScRhGroup.size() );
            hist1d[44]->Fill(jetSCTimeStats[2]);//median
            hist1d[45]->Fill(jetSCTimeStats[0]);//mean
            hist1d[46]->Fill(jetSCTimeStats[4]);//rms
            hist1d[50]->Fill(jetSCTimeStats[5]);//skew
            hist1d[8]->Fill(jetSCTimeStats[6]);//c mean
            hist1d[9]->Fill(jetSCTimeStats[10]);//c med   

			hist1d[89]->Fill(jetSCTimeStats[6]-jetGenTime);

        	hist2d[53]->Fill( scemf, jetemfrac );
			hist2d[56]->Fill( sumscrhe, jeteme );
			hist2d[59]->Fill( dremf, scemf );

			hist2d[63]->Fill( jet.eta(), jetSCTimeStats[7] );  
            hist2d[64]->Fill( jet.etaetaMoment(), jetSCTimeStats[7] );
            hist2d[65]->Fill( jet.phiphiMoment(), jetSCTimeStats[7] );
            hist2d[66]->Fill( jet.etaphiMoment(), jetSCTimeStats[7] );
            hist2d[67]->Fill( jet.maxDistance(), jetSCTimeStats[7] );
            hist2d[68]->Fill( jet.constituentPtDistribution(), jetSCTimeStats[7] );
            hist2d[69]->Fill( jet.constituentEtaPhiSpread(), jetSCTimeStats[7] );
            hist2d[70]->Fill( jet.jetArea(), jetSCTimeStats[7] );
            hist2d[71]->Fill( jet.nCarrying(0.1), jetSCTimeStats[7] );
            hist2d[72]->Fill( jet.nConstituents(), jetSCTimeStats[7] );

        	if( DEBUG ) std::cout << " - fill vars " << std::endl;
            jetSCMedTime.push_back(jetSCTimeStats[2]);
            jetSCMuTime.push_back(jetSCTimeStats[0]);
            std::cout << "fill phCSCMuTimeTemp : " << jetSCTimeStats[6] << std::endl;
            jetCSCMuTime.push_back(jetSCTimeStats[6]);
            jetCSCMedTime.push_back(jetSCTimeStats[10]);

			jetBcTimesCnt.push_back(nBCTimes);
			jetBcSumRHEnr.push_back(sumbcrhe);
			jetBcEMFr.push_back(bcemf);
			jetBcRhCnt.push_back(bcRhCnt);
			jetBcGrpCnt.push_back(bcRhGroups.size());


			if( nBCTimes != 0 && bcemf > minEmf ){
				if( rhCount >= minRHcnt ) 
				hist2d[51]->Fill(rhCount,bcRhCnt);
				hist1d[55]->Fill(bcRhGroups.size());
				auto jetBCTimeStats = getTimeDistStats( bcTimes, bcRhGrpEnergy );
            	auto jetBCTimeStats = getTimeDistStats( bcTimes, bcEnergies );
				auto jetBCTimeStats = getTimeDistStats( bcTimes );
				if( jetBCTimeStats[6] > -28.9 ) nGoodBcJets++;
				jetCBCMedTime.push_back(jetBCTimeStats[10]);c med
	        	jetCBCMuTime.push_back(jetBCTimeStats[6]);c mu
				hist1d[10]->Fill(jetBCTimeStats[10]);//c med
	        	hist1d[19]->Fill(jetBCTimeStats[6]);//c mu
				std::cout << " - fill dbct hist " << std::endl;
				if( nBCTimes == 1 ){ hist1d[11]->Fill(-3.5); }
				else {<<>>if( nBCTimes == 1 )
					for( uInt ita = 0; ita < nBCTimes; ita++ ){
                 		hist2d[50]->Fill(bcTimes[ita],bcRhGrpEnergy[ita]);
						hist2d[50]->Fill(bcTimes[ita],bcEnergies[ita]);
						for( uInt itb = ita+1; itb < nBCTimes; itb++ ){
							auto dt = getdt(bcTimes[ita],bcTimes[itb]);
							hist1d[11]->Fill(dt);
							hist1d[23]->Fill(dt);
							auto effe = effMean(bcRhGrpEnergy[ita],bcRhGrpEnergy[itb]);
							hist2d[45]->Fill(dt, effe);
							hist2d[46]->Fill(dt, effe);
				}	}	}<<>>if( nBCTimes == 1 ) : else	
				hist2d[54]->Fill( bcemf, jetemfrac );
          		hist2d[57]->Fill( sumbcrhe, jeteme );
           		hist2d[58]->Fill( sumbcrhe, sumscrhe );
            	hist2d[60]->Fill( scemf, bcemf );

			
		//****************************	photon/electron to kid pfcand -> SC matcher ********************************

        auto jetERatio = jet.energy()/jetGenEnergy;
    	auto difSCTime = std::abs(jetGenTime-jetSCTime);
        auto difDrTime = std::abs(jetGenTime-jetDrTime);
        //if( jetSCTime < -27.9 ) jetSCTime = -25.0;
        //if( jetDrTime < -27.9 ) jetDrTime = -25.0;
        if( jetSCTime < -25.0 || jetGenTime < -25.0 ) difSCTime = 100.0;
        if( jetDrTime < -25.0 || jetGenTime < -25.0 ) difDrTime = 100.0;
        auto hasGoodGenTime = jetGenTime > -25.0;
        auto etaCut = std::abs(jetGenEta) < 1.5 && std::abs(jet.eta()) < 1.5;
        auto genEnergyCut = jetGenEnergy > 0.0;
        auto genVarCut = jetGenLLPPurity > 0.88 && jetGenTimeVar < 1;
        auto hasGoodGenSCMatch = difSCTime < 0.8;
        //auto hasGoodGenSCMatch = difSCTime < 20.0;
		auto genNoCut = true;// no cut
        //auto genCutTime = 9.0 - 4.0 * jetERatio;
        //auto genPhSpaceCut = jetGenTime > genCutTime;
        //auto genCutDr = jetGenDrMatch > 0.04;
        //auto genCutSCdiff = difSCTime > 4;
		//auto genDrPhSpaceCut = genCutSCdiff && genCutDr;
		auto hasGoodSCTime = jetSCTime > -14.0;

		if( genNoCut ){

            hist1d[112]->Fill(jetGenTimeLLP);
            hist1d[113]->Fill(jetGenLLPPurity);
            hist2d[114]->Fill(jetGenLLPPurity,jetGenTimeVar);
            hist1d[116]->Fill(jetGenNKids);
            hist2d[115]->Fill(jetGenTimeVar,jetGenNKids);
            hist2d[116]->Fill(jetGenLLPPurity,jetGenNKids);
            hist1d[93]->Fill(jetGenTime);
            hist1d[90]->Fill(jetGenImpactAngle);
            hist1d[105]->Fill(jetGenDrMatch);
            hist1d[108]->Fill(jetGenTimeVar);
            hist1d[109]->Fill(jetGenNextBX);
            hist1d[106]->Fill( difSCTime );
            hist1d[107]->Fill( difDrTime );
            hist2d[101]->Fill( jetERatio, jetGenTime );
            hist2d[117]->Fill( difSCTime, jetGenDrMatch );
            hist2d[118]->Fill( jetGenTime, jetemfrac );
			hist2d[124]->Fill( jetemfrac, jetGenDrMatch );
            hist2d[125]->Fill( jetERatio, jetGenDrMatch );

		}//<<>>if( genSpaceCut )

	    if( hasGoodGenSCMatch && etaCut && genEnergyCut && genVarCut && hasGoodGenTime && hasGoodSCTime ){

            hist2d[98]->Fill( jet.energy(), jetGenEnergy );
            hist2d[99]->Fill( jetERatio, jetGenTime );
            hist2d[110]->Fill( jetERatio, jetSCTime );
            hist2d[111]->Fill( jetERatio, jetDrTime );
            hist2d[100]->Fill( jetemfrac, jetSCTime );
            hist2d[103]->Fill( jetERatio, difSCTime );
            hist2d[104]->Fill( jetDrTime, jetSCTime );
            hist2d[102]->Fill( jetGenTime, jetDrTime );
            hist2d[105]->Fill( jetGenTime, jetSCTime );
            hist2d[106]->Fill( jetGenTime, jetGenEnergy );
            hist2d[107]->Fill( jetGenDrMatch, difSCTime );
            hist2d[108]->Fill( jetGenDrMatch, difDrTime );
            hist2d[112]->Fill( jetGenTimeVar, difSCTime );
            hist2d[113]->Fill( jetGenLLPPurity, difSCTime );

        }//<<>>if( jetSCTimeStats[0] > -28.0 && jetGenTime > -28.0 )

		if( DEBUG ) std::cout << "Next Jet .......................... " << std::endl; 	
	}<<>>for ( uInt ijet = 0; ijet < nJets; ijet++ )
	 ** end of jets	***************************************************************************************************

	


	hist1d[17]->Fill(jetHt);
   	hist1d[30]->Fill(nGoodDrJets);
   	hist1d[31]->Fill(nGoodScJets);
   	hist1d[32]->Fill(nGoodBcJets);
   	hist1d[33]->Fill(nUnJets);
   	if( nUnJets != 0 ) hist1d[34]->Fill(float(nJets)/nUnJets);
   	if( nGoodScJets != 0 ) hist1d[38]->Fill(float(nGoodBcJets)/nGoodScJets);
	if( nJets != 0 ){
   		hist1d[35]->Fill(float(nGoodDrJets)/nJets);
   		hist1d[36]->Fill(float(nGoodScJets)/nJets);
   		hist1d[37]->Fill(float(nGoodBcJets)/nJets);
	}//<<>>if( nUnJets != 0 )
	
	-----------------------------------------------------------------------------------------------------
	 ***************************** d jetTime for back-to-back high pt jets	*****************************
	auto dijetIdCut = 1;
	auto dijetPtMin = 200.0;
	auto difPtLmt = 0.8;
	auto htPctLmt = 0.8;
	auto dPhiLmt = 2.8;
	

  	if( jetCSCMuTime.size() != jetCBCMuTime.size() ) hist1d[53]->Fill(-3.25);
	else if( jetCSCMuTime.size() != nJets ) hist1d[53]->Fill(-3.0);
	else for( uInt q = 0; q < nJets; q++ ) hist1d[53]->Fill( getdt(jetCSCMuTime[q],jetCBCMuTime[q]) ); 

	if( DEBUG ) std::cout << "Finding jet dt pairs" << std::endl;
	for ( uInt q = 0; q < nJets; q++ ){
		for ( uInt p = q+1; p < nJets; p++ ){
	
			if( DEBUG ) std::cout << " - filter jet pairs" << std::endl;
	      	const auto & qjet = fjets[q];
	      	const auto & pjet = fjets[p];
	      	if( qjet.pt() < dijetPtMin ) continue;
	      	auto diffPt = pjet.pt()/qjet.pt();
	      	hist1d[24]->Fill(diffPt);
	      	if( diffPt < difPtLmt ) continue;
	      	auto htPct= (qjet.pt()+pjet.pt())/jetHt;
			hist1d[25]->Fill(htPct);
	      	if( htPct < htPctLmt ) continue;
	      	auto dPhi = reco::deltaPhi(qjet.phi(),pjet.phi());
	      	hist1d[26]->Fill(dPhi);
	      	if( dPhi < dPhiLmt ) continue;

			if( DEBUG ) std::cout << " - get jet pair dt" << std::endl;
			auto dTmu = getdt( jetMuTime[q], jetMuTime[p] );
            auto dTmed = getdt( jetMedTime[q], jetMedTime[p] );
            auto dTcmu = getdt( jetCMuTime[q], jetCMuTime[p] );
            auto dTcmed = getdt( jetCMedTime[q], jetCMedTime[p] );
            if( DEBUG ) std::cout << "dT dR      : " << dTmu <<  " " << jetMuTime[q] << " " << jetMuTime[p] << std::endl;
	      	auto dTmusc = getdt( jetSCMuTime[q], jetSCMuTime[p] );
	      	auto dTmedsc = getdt( jetSCMedTime[q], jetSCMedTime[p] );
	      	auto dTcmusc = getdt( jetCSCMuTime[q], jetCSCMuTime[p] );
	      	auto dTcmedsc = getdt( jetCSCMedTime[q], jetCSCMedTime[p] );
            if( DEBUG ) std::cout << "dT SC      : " << dTmusc <<  " " << jetSCMuTime[q] << " " << jetSCMuTime[p] << std::endl;
	      	auto dTcmubc = getdt( jetCBCMuTime[q], jetCBCMuTime[p] );
            if( DEBUG ) std::cout << "dT cMu BC  : " << dTcmubc <<  " " << jetCBCMuTime[q] << " " << jetCBCMuTime[p] << std::endl;
	      	auto dTcmedbc = getdt( jetCBCMedTime[q], jetCBCMedTime[p] );
            if( DEBUG ) std::cout << "dT cMed BC : " << dTcmedbc <<  " " << jetCBCMedTime[q] << " " << jetCBCMedTime[p] << std::endl;
            auto dTmuph = getdt( jetPhMuTime[q], jetPhMuTime[p] );
            if( DEBUG ) std::cout << "dT Ph      : " << dTmuph <<  " " << jetPhMuTime[q] << " " << jetPhMuTime[p] << std::endl;
            auto dTmuel = getdt( jetEleMuTime[q], jetEleMuTime[p] );
            if( DEBUG ) std::cout << "dT Ele     : " << dTmuel <<  " " << jetEleMuTime[q] << " " << jetEleMuTime[p] << std::endl;

   	//<<<<for ( uInt q = 0; q < nJets; q++ ){
      	//<<<<for ( uInt p = q+1; p < nJets; p++ ){

			if( DEBUG ) std::cout << " - fill hists" << std::endl;

			auto dtThrs = -2.5;//removes default dt values from getdt
			if( dTmu > dtThrs ) hist1d[15]->Fill(dTmu);
	      	if( dTmed > dtThrs ) hist1d[16]->Fill(dTmed);
            if( dTcmu > dtThrs ) hist1d[39]->Fill(dTcmu);
            if( dTcmed > dtThrs ) hist1d[40]->Fill(dTcmed);

	      	if( dTmedsc > dtThrs ) hist1d[48]->Fill(dTmedsc);
	      	if( dTmusc > dtThrs ) hist1d[49]->Fill(dTmusc);
	      	if( dTcmusc > dtThrs ) hist1d[41]->Fill(dTcmusc);
	      	if( dTcmedsc > dtThrs ) hist1d[42]->Fill(dTcmedsc);

	      	if( dTcmubc > dtThrs ) hist1d[27]->Fill(dTcmubc);
	      	if( dTcmedbc > dtThrs ) hist1d[28]->Fill(dTcmedbc);

            if( dTmuph > dtThrs ) hist1d[59]->Fill(dTmuph);
            if( dTmuel > dtThrs*2 ) hist1d[60]->Fill(dTmuel);

			hist2d[22]->Fill(dTmu,nJets);
	      	hist2d[26]->Fill(dTmu,diffPt);
	      	hist2d[27]->Fill(dTmu,htPct);
	      	hist2d[28]->Fill(dTmu,dPhi);
	     	hist2d[23]->Fill(dTmed,nJets);
	      	hist2d[29]->Fill(dTmu,diffPt);
	      	hist2d[30]->Fill(dTmu,htPct);
	      	hist2d[31]->Fill(dTmu,dPhi);

			if( DEBUG ) std::cout << " - fill dt vs eff e hists" << std::endl;
	      	auto effje = effMean(jetPHE[p],jetPHE[q]);
	      	hist2d[43]->Fill(dTmusc,effje);
	      	hist2d[44]->Fill(dTmu,effje);

		}//<<>>for ( uInt p = q+1; p < nJets; p++ )
	}//<<>>for ( uInt q = 0; q < nJets; q++ )
	-------------------------------------------------------------------------------

	if( goodJetEvent ) nGoodJetEvents++;

	 -- Fill output trees ------------------------------------------
void LLPgammaAnalyzer_AOD::beginJob(){

	 Global Varibles
	nGoodJetEvents = 0;

   	 Set up DetIdMap
   	DetIDMap = SetupDetIDs();

	 Book output files and trees
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("llpgtree","llpgtree");

	 Book histograms
	
	int jtdiv(400);
	float jtran(8);
    int jtdiv(625);
    float jtran(25);
	int jdtdiv(200);
	float jdtran(4);
   	int jztdiv(100);
   	float jztran(2);
	int rhcnt(80);

   	------ 1D Hists --------------------------------------------------------------------------

	hist1d[0] = fs->make<TH1D>("jetRHTime", "jetRHTime", 2000, -100, 100);
	hist1d[29] = fs->make<TH1D>("jetMuTime", "jetMuTime", jtdiv, -1*jtran, jtran);
	hist1d[1] = fs->make<TH1D>("jetRHMulti", "jetRHMulti", rhcnt, 0, rhcnt);
	hist1d[2] = fs->make<TH1D>("jetTimeError", "jetTimeError", 300, 0, 3);
	hist1d[3] = fs->make<TH1D>("jetTimeRMS", "jetTimeRMS", 200, 0, 20);
	hist1d[4] = fs->make<TH1D>("jetMedTime", "jetMedTime", jtdiv, -1*jtran, jtran);
	hist1d[5] = fs->make<TH1D>("jetRawTime", "jetRawTime", jtdiv, -1*jtran, jtran);
	hist1d[6] = fs->make<TH1D>("jetCMuTime", "jetCMuTime", jtdiv, -1*jtran, jtran);
	hist1d[7] = fs->make<TH1D>("jetCMedTime", "jetCMedTime", jtdiv, -1*jtran, jtran);
	hist1d[8] = fs->make<TH1D>("jetCSCMuTime", "jetCSCMuTime", jtdiv, -1*jtran, jtran);
	hist1d[9] = fs->make<TH1D>("jetCSCMedTime", "jetCSCMedTime", jtdiv, -1*jtran, jtran);
	hist1d[10] = fs->make<TH1D>("jetCBCMedTime", "jetCBCMedTime", jtdiv, -1*jtran, jtran);

	hist1d[12] = fs->make<TH1D>("jetPt", "jetPt", 500, 0, 5000);
	hist1d[13] = fs->make<TH1D>("jetPhi", "jetPhi", 700, -3.5, 3.5);
	hist1d[14] = fs->make<TH1D>("jetEta", "jetEta", 700, -3.5, 3.5);
	hist1d[15] = fs->make<TH1D>("jetdtmu", "jetdtmu", jdtdiv, -1*jdtran, jdtran);
	hist1d[16] = fs->make<TH1D>("jetdtmed", "jetdtmed", jdtdiv, -1*jdtran, jdtran);
	hist1d[17] = fs->make<TH1D>("jetHt", "jetHt", 1000, 0, 5000);
	hist1d[18] = fs->make<TH1D>("nJet", "nJets", 21, -0.5, 20.5);
	hist1d[19] = fs->make<TH1D>("jetCBCMuTime", "jetCBCMuTime", jtdiv, -1*jtran, jtran);

	hist1d[24] = fs->make<TH1D>("diffPt", "diffPt", 1000, 0, 10);
	hist1d[25] = fs->make<TH1D>("htPct", "htPct", 100, 0, 1);
	hist1d[26] = fs->make<TH1D>("dPhi", "dPhi", 70, -3.5, 3.5);

	hist1d[27] = fs->make<TH1D>("jetcmudtbc", "jetcmudtbc", jdtdiv, -1*jdtran, jdtran);
	hist1d[28] = fs->make<TH1D>("jetcmeddtbc", "jetcmeddtbc", jdtdiv, -1*jdtran, jdtran);
	 ----  moved to 2nd from top of list : hist1d[29] ( after hist1d[0] )
	hist1d[30] = fs->make<TH1D>("nGoodDrJets", "nGoodDrJets", 21, -0.5, 20.5);
	hist1d[31] = fs->make<TH1D>("nGoodScJets", "nGoodScJets", 21, -0.5, 20.5);
	hist1d[32] = fs->make<TH1D>("nGoodBcJets", "nGoodBcJets", 21, -0.5, 20.5);
	hist1d[33] = fs->make<TH1D>("nUnJets", "nUnJets", 101, -0.5, 100.5);
	hist1d[34] = fs->make<TH1D>("pJets", "pJets", 110, 0, 1.1);
	hist1d[35] = fs->make<TH1D>("pGoodDrJets", "pGoodDrJets", 110, 0, 1.1);
	hist1d[36] = fs->make<TH1D>("pGoodScJets", "pGoodScJets", 110, 0, 1.1);
	hist1d[37] = fs->make<TH1D>("pGoodBcJets", "pGoodBcJets", 110, 0, 1.1);
	hist1d[38] = fs->make<TH1D>("pGoodBcToScJets", "pGoodBcToScJets", 110, 0, 1.1);

	hist1d[39] = fs->make<TH1D>("jetcmudt", "jetcmudt", jdtdiv, -1*jdtran, jdtran);
	hist1d[40] = fs->make<TH1D>("jetcmeddt", "jetcmeddt", jdtdiv, -1*jdtran, jdtran);
	hist1d[41] = fs->make<TH1D>("jetcmudtsc", "jetcmudtsc", jdtdiv, -1*jdtran, jdtran);
	hist1d[42] = fs->make<TH1D>("jetcmeddtsc", "jetcmeddtsc", jdtdiv, -1*jdtran, jdtran);

	hist1d[43] = fs->make<TH1D>("nPhotonsPerJet","nPhotonsPerJet", 21, -0.5, 20.5);

	hist1d[44] = fs->make<TH1D>("jetSCmedTime", "jetSCmedTime", jtdiv, -1*jtran, jtran);
	hist1d[45] = fs->make<TH1D>("jetSCmuTime", "jetSCmuTime", jtdiv, -1*jtran, jtran);
	hist1d[46] = fs->make<TH1D>("jetSCTimeRms", "jetSCTimeRms", 60, 0, 6);
	hist1d[47] = fs->make<TH1D>("jetSCrawTime", "jetSCrawTime", jtdiv, -1*jtran, jtran);

	hist1d[48] = fs->make<TH1D>("jetmeddtsc", "jetmeddtsc", jdtdiv, -1*jdtran, jdtran);
	hist1d[49] = fs->make<TH1D>("jetmudtsc", "jetmudtsc", jdtdiv, -1*jdtran, jdtran);

	hist1d[50] = fs->make<TH1D>("jetSCTimeSkew", "jetSCTimeSkew", 80, -4.0, 4.0);
	hist1d[51] = fs->make<TH1D>("jetPhotons", "jetPhotons", 21, -0.5, 20.5);
	hist1d[52] = fs->make<TH1D>("jetElectrons", "jetElectrons", 21, -0.5, 20.5);

   	hist1d[53] = fs->make<TH1D>("scbcdt", "scbcdt", jdtdiv, -1*jdtran, jdtran);
   	hist1d[54] = fs->make<TH1D>("bc1rhef", "bc1rhef", 110, 0, 1.1);
   	hist1d[55] = fs->make<TH1D>("nBCinJet", "nBCinJet", 11, -0.5, 10.5);
   	hist1d[56] = fs->make<TH1D>("bcMrhef", "bcMrhef", 110, 0, 1.1);

    hist1d[57] = fs->make<TH1D>("jetCPhMuTime", "jetCPhMuTime", jtdiv, -1*jtran, jtran);
    hist1d[58] = fs->make<TH1D>("jetCPhMedTime", "jetCPhMedTime", jtdiv, -1*jtran, jtran);
    hist1d[59] = fs->make<TH1D>("jetmudtph", "jetmudtph", jdtdiv, -1*jdtran, jdtran);
    hist1d[60] = fs->make<TH1D>("jetmudtel", "jetmudtel", jdtdiv, -1*jdtran, jdtran);
    hist1d[61] = fs->make<TH1D>("jetCEleMuTime", "jetCEleMuTime", jtdiv, -1*jtran, jtran);
    hist1d[62] = fs->make<TH1D>("jetCEleMedTime", "jetCEleMedTime", jtdiv, -1*jtran, jtran);

    hist1d[63] = fs->make<TH1D>("scEtaPhiAngle3D", "scEtaPhiAngle3D", 660, -0.2, 6.4);
    hist1d[64] = fs->make<TH1D>("scEtaTimAngle3D", "scEtaTimAngle3D", 660, -0.2, 6.4);
    hist1d[65] = fs->make<TH1D>("scEtaPhiAngle2D", "scEtaTimAngle2D", 660, -0.2, 6.4);

    hist1d[66] = fs->make<TH1D>("sciEta3D", "sciEta", 171, -85, 85);
    hist1d[67] = fs->make<TH1D>("sciPhi3D", "sciPhi", 361, 0, 360);
    hist1d[68] = fs->make<TH1D>("scTime3D", "scTime", 5000, -25, 25);
    hist1d[69] = fs->make<TH1D>("sciEta2D", "sciEtaDiff", 201, -100, 100 );
    hist1d[70] = fs->make<TH1D>("sciPhi2D", "sciPhiDiff", 201, -100, 100);
    hist1d[71] = fs->make<TH1D>("sciTim2D", "sciTimDiff", 400, -10, 10);
    hist1d[72] = fs->make<TH1D>("sciAngle2D", "sciAngleSph", 660, -0.2, 6.4);
    hist1d[73] = fs->make<TH1D>("scRotAngle", "scRotAngle", 660, -0.2, 6.4);
    hist1d[74] = fs->make<TH1D>("scAngleTest", "scAngleTest", 660, -0.2, 6.4);
    hist1d[75] = fs->make<TH1D>("sciEta3Diff", "sciEta3Diff", 400, -10, 10 );
    hist1d[76] = fs->make<TH1D>("sciPhi3Diff", "sciPhi3Diff", 400, -10, 10);
    hist1d[77] = fs->make<TH1D>("sciTim3Diff", "sciTim3Diff", 400, -10, 10);
    hist1d[78] = fs->make<TH1D>("scSphEgn0", "scSphEgn0", 400, -10, 10);
    hist1d[79] = fs->make<TH1D>("scSphEgn1", "scSphEgn1", 400, -10, 10);
    hist1d[80] = fs->make<TH1D>("scSinTest", "scSinTest", 200, -1, 1);
    hist1d[81] = fs->make<TH1D>("eginValueSph", "eginValueSph", 150, 0.4, 1.1);
    hist1d[82] = fs->make<TH1D>("dIPhiTest", "dIPhiTest", 1440, -720, 720);
    hist1d[83] = fs->make<TH1D>("meanIPhiTest", "meanIPhiTest", 1440, -720, 720);
    hist1d[84] = fs->make<TH1D>("meanIEtaTest", "meanIEtaTest", 200, -100, 100);
    hist1d[85] = fs->make<TH1D>("meanTimeTest", "meanTimeTest", 2000, -25, 25);
    hist1d[86] = fs->make<TH1D>("eginValue3D", "eginValue3D", 110, 0, 1.1);

     see below in 2d hists for declaration, commeted here for clarity
    hist1d[87] = fs->make<TH1D>("cluster_etprofile", "Cluster Eta Time Profile Sph", cwdiv, -1*cwtrn, cwtrn);
    hist1d[88] = fs->make<TH1D>("cluster_et3Dprofl", "Cluster Eta Time Profile 3D", cl3ddiv1, -1*cl3dtrn1, cl3dtrn1);

    hist1d[89] = fs->make<TH1D>("jetmudtgen", "jetmudtgen", jdtdiv, -1*jdtran, jdtran);
    hist1d[90] = fs->make<TH1D>("genJetImpactAngle", "genJetImpactAngle", 660, -0.2, 6.4);

    hist1d[91] = fs->make<TH1D>("clEtaTimeSlope", "Cluster Eta Time Slope", 240, -120, 120);
    hist1d[92] = fs->make<TH1D>("clEtaTimeSlopeChi", "Cluster Eta Time Slope Chi2", 100, 0, 1);

    hist1d[93] = fs->make<TH1D>("jetGenTime", "jetGenTime", jtdiv, -1*jtran, jtran);

	hist1d[94] = fs->make<TH1D>("clBinProfile_m1", "spCluster Bin -1 Profile", 200, -10, 10);
    hist1d[95] = fs->make<TH1D>("clBinProfile_p1", "spCluster Bin +1 Profile", 200, -10, 10);
    hist1d[96] = fs->make<TH1D>("clBinProfile_m2", "spCluster Bin -2 Profile", 200, -10, 10);
    hist1d[97] = fs->make<TH1D>("clBinProfile_p2", "spCluster Bin +2 Profile", 200, -10, 10);
    hist1d[98] = fs->make<TH1D>("clBinProfile_m3", "spCluster Bin -3 Profile", 200, -10, 10);
    hist1d[99] = fs->make<TH1D>("clBinProfile_p3", "spCluster Bin +3 Profile", 200, -10, 10);

    hist1d[100] = fs->make<TH1D>("clETSlopes", "Cluster ET Slopes", 240, -120, 120);
	101 used below
	102 used blow

    hist1d[114] = fs->make<TH1D>("clETSlope3D", "Cluster EtaTimeSlope 3D", 500, -250, 250);
    hist1d[115] = fs->make<TH1D>("clETSlopeChi3D", "Cluster EtaTimeSlope Chi2 3D", 100, 0, 1);	
	
    hist1d[103] = fs->make<TH1D>("clEtaTimeSlopeInv", "Cluster Eta Time SlopeInv", 350, -0.1, 34.9);
    hist1d[104] = fs->make<TH1D>("clEtaTimeSlopeInv3D", "Cluster Eta Time SlopeInv 3D", 350, -0.1, 34.9);

    hist1d[105] = fs->make<TH1D>("genJetDrMatchJet", "genJetDrMatchJet", 100, 0, 0.1);
    hist1d[106] = fs->make<TH1D>("genJetSCTimeDiff", "genJetSCTimeDiff", 300, 0, 30);
    hist1d[107] = fs->make<TH1D>("genJetDrTimeDiff", "genJetSCTimeDiff", 300, 0, 30);

    hist1d[108] = fs->make<TH1D>("jetGenTimeVar", "jetGenTimeVar", 408, -2, 100);
    hist1d[109] = fs->make<TH1D>("jetGenTimeNextBX", "jetGenTimeNextBX", 3, -1, 2);
    hist1d[110] = fs->make<TH1D>("jetGenTimeNoTOF", "jetGenTimeNoTOF", 300, 0, 30);
    hist1d[111] = fs->make<TH1D>("jetGenTOF", "jetGenTOF", 300, 0, 30);
    hist1d[112] = fs->make<TH1D>("jetGenTimeIsLLP", "jetGenTimeIsLLP", 3, -1, 2);
    hist1d[113] = fs->make<TH1D>("jetGenTimeLLPPurity", "jetGenTimeLLPPurity", 100, 0, 1);
	hist1d[114] used above
	hist1d[115] usded above
	hist1d[116] = fs->make<TH1D>("jetGenNKids", "jetGenNKids", 100, 0, 100);
    hist1d[117] = fs->make<TH1D>("clEtaTimeSlopeRangeA", "Cluster Eta Time Slope Eta > 1.0", 480, -240, 240);
    hist1d[118] = fs->make<TH1D>("clEtaTimeSlopeRangeB", "Cluster Eta Time Slope Eta < 0.5", 480, -240, 240);

    hist1d[119] = fs->make<TH1D>("ootPhotonTime", "ootPhotonTime", jtdiv, -1*jtran, jtran);

    hist1d[120] = fs->make<TH1D>("jetCOOTPhMuTime", "jetCOOTPhMuTime", jtdiv, -1*jtran, jtran);
    hist1d[121] = fs->make<TH1D>("jetCOOTPhMedTime", "jetCOOTPhMedTime", jtdiv, -1*jtran, jtran);

	hist1d[122] = fs->make<TH1D>("jetPhClRhTime", "phClRhTime", jtdiv, -1*jtran, jtran);
    hist1d[123] = fs->make<TH1D>("jetPhClRhPkOOT", "phClRhPkOOT", 120, -0.1, 1.1);
    hist1d[124] = fs->make<TH1D>("jetPhClRhPMatched", "phClRhPMatched", 120, -0.1, 1.1);
    hist1d[125] = fs->make<TH1D>("jetOOTPhClRhTime", "ootPhClRhTime", jtdiv, -1*jtran, jtran);
    hist1d[126] = fs->make<TH1D>("jetOOTPhClRhPkOOT", "ootPhClRhPkOOT", 120, -0.1, 1.1);
    hist1d[127] = fs->make<TH1D>("jetOOTPhClRhPMatched", "ootPhClRhPMatched", 120, -0.1, 1.1);
    hist1d[128] = fs->make<TH1D>("jetEleClRhTime", "eleClRhTime", jtdiv, -1*jtran, jtran);
    hist1d[129] = fs->make<TH1D>("jetEleClRhPkOOT", "eleClRhPkOOT", 120, -0.1, 1.1);
    hist1d[130] = fs->make<TH1D>("jetEleClRhPMatched", "eleClRhPMatched", 120, -0.1, 1.1);

    hist1d[131] = fs->make<TH1D>("ebRhTime", "ebRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[132] = fs->make<TH1D>("ebRhEnergy", "ebRhEnergy", 1000, 0, 1000);
    hist1d[133] = fs->make<TH1D>("eeRhTime", "eeRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[134] = fs->make<TH1D>("eeRhEnergy", "eeRhEnergy", 1000, 0, 1000);
    hist1d[135] = fs->make<TH1D>("ebRhkOOT", "ebRhkOOT", 3, 0, 1);
    hist1d[136] = fs->make<TH1D>("eeRhkOOT", "eeRhkOOT", 3, 0, 1);

    hist1d[137] = fs->make<TH1D>("phClRhTime", "phClRhTime", jtdiv, -1*jtran, jtran);
    hist1d[138] = fs->make<TH1D>("phClRhPkOOT", "phClRhPkOOT", 120, -0.1, 1.1);
    hist1d[139] = fs->make<TH1D>("phClRhPMatched", "phClRhPMatched", 120, -0.1, 1.1);
    hist1d[140] = fs->make<TH1D>("ootPhClRhTime", "ootPhClRhTime", jtdiv, -1*jtran, jtran);
    hist1d[141] = fs->make<TH1D>("ootPhClRhPkOOT", "ootPhClRhPkOOT", 120, -0.1, 1.1);
    hist1d[142] = fs->make<TH1D>("ootPhClRhPMatched", "ootPhClRhPMatched", 120, -0.1, 1.1);
    hist1d[143] = fs->make<TH1D>("eleClRhTime", "eleClRhTime", jtdiv, -1*jtran, jtran);
    hist1d[144] = fs->make<TH1D>("eleClRhPkOOT", "eleClRhPkOOT", 120, -0.1, 1.1);
    hist1d[145] = fs->make<TH1D>("eleClRhPMatched", "eleClRhPMatched", 120, -0.1, 1.1);

    hist1d[146] = fs->make<TH1D>("phClTime", "phClTime", jtdiv, -1*jtran, jtran);
    hist1d[147] = fs->make<TH1D>("phSeedRhTime", "phLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[148] = fs->make<TH1D>("phClSeedTimeDiff", "phClLeadTimeDiff", jtdiv, -1*jtran, jtran);
    hist1d[149] = fs->make<TH1D>("ootPhClTime", "ootPhClTime", jtdiv, -1*jtran, jtran);
    hist1d[150] = fs->make<TH1D>("ootPhSeedRhTime", "ootPhLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[151] = fs->make<TH1D>("ootPhClSeedTimeDiff", "ootPhClLeadTimeDiff", jtdiv, -1*jtran, jtran);
    hist1d[152] = fs->make<TH1D>("eleClTime", "eleClTime", jtdiv, -1*jtran, jtran);
    hist1d[153] = fs->make<TH1D>("eleSeedRhTime", "eleLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[154] = fs->make<TH1D>("eleClSeedTimeDiff", "eleClLeadTimeDiff", jtdiv, -1*jtran, jtran);

	------ 2D Hists --------------------------------------------------------------------------

   hist2d[0] = fs->make<TH2D>("jt_pt", "jt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[1] = fs->make<TH2D>("jt_pt", "jt_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[2] = fs->make<TH2D>("jt_id", "jt_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
	hist2d[3] = fs->make<TH2D>("jt_nhf", "jt_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[4] = fs->make<TH2D>("jt_chf", "jt_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[5] = fs->make<TH2D>("jt_nemf", "jt_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[6] = fs->make<TH2D>("jt_cemf", "jt_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[7] = fs->make<TH2D>("jt_muf", "jt_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[8] = fs->make<TH2D>("jt_nhm", "jt_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
	hist2d[9] = fs->make<TH2D>("jt_chm", "jt_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

	hist2d[10] = fs->make<TH2D>("jt_medt", "jt_medt", jtdiv, -1*jtran, jtran, 200, -10, 10);
	hist2d[11] = fs->make<TH2D>("jt_rms", "jt_rms", jtdiv, -1*jtran, jtran, 200, 0, 20);
	hist2d[12] = fs->make<TH2D>("jt_err", "jt_err", jtdiv, -1*jtran, jtran, 300, 0, 3);

	hist2d[13] = fs->make<TH2D>("medt_pt", "medt_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[14] = fs->make<TH2D>("medt_id", "medt_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
	hist2d[15] = fs->make<TH2D>("medt_nhf", "medt_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[16] = fs->make<TH2D>("medt_chf", "medt_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[17] = fs->make<TH2D>("medt_nemf", "medt_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[18] = fs->make<TH2D>("medt_cemf", "medt_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[19] = fs->make<TH2D>("medt_muf", "medt_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[20] = fs->make<TH2D>("medt_nhm", "medt_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
	hist2d[21] = fs->make<TH2D>("medt_chm", "medt_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

	hist2d[22] = fs->make<TH2D>("jdtmu_nJets", "jdtmu_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);
	hist2d[23] = fs->make<TH2D>("jdtmed_nJets", "jdtmed_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);

	hist2d[24] = fs->make<TH2D>("jt_nrh", "jt_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);
	hist2d[25] = fs->make<TH2D>("medt_nrh", "medt_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);

	hist2d[26] = fs->make<TH2D>("jdtmu_diffPt", "jdtmu_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[27] = fs->make<TH2D>("jdtmu_htPct", "jdtmu_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[28] = fs->make<TH2D>("jdtmu_dPhi", "jdtmu_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);
	hist2d[29] = fs->make<TH2D>("jdtmed_diffPt", "jdtmed_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[30] = fs->make<TH2D>("jdtmed_htPct", "jdtmed_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[31] = fs->make<TH2D>("jdtmed_dPhi", "jdtmed_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);

	hist2d[32] = fs->make<TH2D>("jt_sceta", "jt_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[33] = fs->make<TH2D>("jt_scphi","jt_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[34] = fs->make<TH2D>("jt_scenr", "jt_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

	hist2d[35] = fs->make<TH2D>("medt_sceta", "medt_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[36] = fs->make<TH2D>("medt_scphi","medt_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[37] = fs->make<TH2D>("medt_scenr", "medt_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

	hist2d[38] = fs->make<TH2D>("rht_rhe", "rht_rhe", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

	hist2d[39] = fs->make<TH2D>("njrh_nscrh", "njrh_nscrh", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);
	hist2d[40] = fs->make<TH2D>("je_sce", "je_sce", 500, 0, 500, 500, 0, 500);
	hist2d[41] = fs->make<TH2D>("je_ege", "je_ege", 500, 0, 500, 500, 0, 500);
	hist2d[42] = fs->make<TH2D>("ege_sce", "ege_sce", 500, 0, 500, 500, 0, 500);
	hist2d[43] = fs->make<TH2D>("scdt_effje", "scdt_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);
	hist2d[44] = fs->make<TH2D>("jdt_effje", "jdt_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);

   	hist2d[51] = fs->make<TH2D>("njrh_nsbcrh", "njrh_nsbcrh", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);

   	hist2d[52] = fs->make<TH2D>("dremf_emf", "dremf_emf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[53] = fs->make<TH2D>("scemf_emf", "scemf_emf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[54] = fs->make<TH2D>("bcemf_emf", "bcemf_emf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[55] = fs->make<TH2D>("dreme_eme", "dreme_eme", 500, 0, 500, 500, 0, 500);
   	hist2d[56] = fs->make<TH2D>("sceme_eme", "sceme_eme", 500, 0, 500, 500, 0, 500);
   	hist2d[57] = fs->make<TH2D>("bceme_eme", "bceme_eme", 500, 0, 500, 500, 0, 500);
   	hist2d[58] = fs->make<TH2D>("sce_bce", "sce_bce", 500, 0, 500, 500, 0, 500);
   	hist2d[59] = fs->make<TH2D>("dremf_scemf", "dremf_scemf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[60] = fs->make<TH2D>("scemf_bcemf", "scemf_bcemf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[61] = fs->make<TH2D>("epaf_epf", "epaf_epf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[62] = fs->make<TH2D>("epf_emf", "epf_emf", 110, 0, 1.1, 110, 0, 1.1);

	auto stddiv = 120;
    auto stdtran = 3;
    hist2d[63] = fs->make<TH2D>("jetEta_stdSCt", "jetEta_stdSCt", 200, -2.5, 2.5, stddiv, 0, stdtran);
    hist2d[64] = fs->make<TH2D>("jetEtaetaMmt_stdSCt", "jetEtaetaMmt_stdSCt", 800, 0, 0.02, stddiv, 0, stdtran);
    hist2d[65] = fs->make<TH2D>("jetPhiphiMnt_stdSCt", "jetPhiphiMnt_stdSCt", 800, 0, 0.02, stddiv, 0, stdtran);
    hist2d[66] = fs->make<TH2D>("jetEtaphiMnt_stdSCt", "jetEtaphiMnt_stdSCt", 1600, -0.02, 0.02, stddiv, 0, stdtran);
    hist2d[67] = fs->make<TH2D>("jetMaxD_stdSCt", "jetMaxD_stdSCt", 400, 0, 1, stddiv, 0, stdtran);
    hist2d[68] = fs->make<TH2D>("jetConPtDis_stdSCt", "jetConPtDis_stdSCt", 400, 0, 1, stddiv, 0, stdtran);
    hist2d[69] = fs->make<TH2D>("jetConEtaPhiSprd_stdSCt", "jetConEtaPhiSprd_stdSCt", 600, -0.005, 0.01, stddiv, 0, stdtran);
    hist2d[70] = fs->make<TH2D>("jetArea_stdSCt", "jetArea_stdSCt", 200, 0, 1, stddiv, 0, stdtran);
    hist2d[71] = fs->make<TH2D>("jetNCarry_stdSCt", "jetNCarry_stdSCt", 40, 0, 40, stddiv, 0, stdtran);
    hist2d[72] = fs->make<TH2D>("jetNConst_stdSCt", "jetNConst_stdSCt", 100, 0, 100, stddiv, 0, stdtran);

	auto cldiv = 80;
	auto cltrn = 20;

    auto cl3ddiv = 200;
    auto cl3dtrn = 4;
    auto cl3ddiv1 = 200;
    auto cl3dtrn1 = 4;

    auto clsphdiv = 400;
    auto clsphtrn = 4;
    auto cwdiv = 80;
    auto cwtrn = 40;
	auto slmax = 120;
	auto slmin = -120;
	auto sldiv = 120;
    auto sl3max = 12;
    auto sl3min = -12;
    auto sl3div = 2400;
    auto sldiv = 320;
	auto chimax = 1.01;
    auto chimin = 0.91;
	auto chidiv = 400;

    hist2d[89] = fs->make<TH2D>("jetEtavSlope", "Jet Eta v Slope;Eta;Slope ps/cm", 60, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[90] = fs->make<TH2D>("jetImpAnglevSlope", "Jet ImpactAngle v Slope;ImpactAngle;Slope ps/cm", 60, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[91] = fs->make<TH2D>("clEtaTimeSlvChi", "Cluster EtaTime Slope v Chi2;Slope;Chi2", sldiv, slmin, slmax, chidiv, chimin, chimax);
    hist2d[92] = fs->make<TH2D>("clEtaTimeSlvEVal", "Cluster EtaTime Slope v EigenValue;Slope;EigenValue", sldiv, slmin, slmax, 180, 0.55, 1.0);
    hist2d[93] = fs->make<TH2D>("clEtaTimeSlvRotAng", "Cluster EtaTime Slope v Rotation Angle;Slope;rotAngle", sldiv, slmin, slmax, 660, -0.2, 6.4);
    hist2d[94] = fs->make<TH2D>("clEtaTimeSlvNumClRHs", "Cluster EtaTime Slope v nClRecHits;Slope;nClRecHits", sldiv, slmin, slmax, 60, 0, 60);
    hist2d[95] = fs->make<TH2D>("clEtaTimeChi2vEVal", "Cluster EtaTime Chi2 v EigenValue;Chi2;EigenValue", chidiv, chimin, chimax, 60, 0.45, 1.05);
    hist2d[96] = fs->make<TH2D>("jetImpAnglevSlope3D", "Jet ImpactAngle v Slope 3D;ImpactAngle;Slope", 150, 0, 1.5, sl3div, sl3min, sl3max);
    hist2d[97] = fs->make<TH2D>("clEtaTimeChi2vNumClRHs", "Cluster EtaTime Chi2 v nClRecHits;Chi2;nClRecHits", chidiv, chimin, chimax, 60, 0, 60);

    hist2d[98] = fs->make<TH2D>("jetEvGenE", "Jet Energy v GenEnergy;JetEnergy;GenEnergy", 100, 0, 1000, 100, 0, 1000 );
    hist2d[99] = fs->make<TH2D>("jetEGenERatiovGenTime", "Jet E/GenE v GenTime;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[100] = fs->make<TH2D>("jetEMFracvSCTime", "Jet EMFrac v SCTime;EMFrac;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[101] = fs->make<TH2D>("jetGenRatiovGenTimePre", "Jet E/GenE v GenTime Pre;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[102] = fs->make<TH2D>("jetGenTimevDrJetTime", "GenTime v DrJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[103] = fs->make<TH2D>("jetEGenERatiovSCTimeDiff", "Jet E/GenE v JetSC GenJet TimeDif;E/GenE;SCTimeDif", 80, 0, 2, 300, 0, 30.0 );
    hist2d[104] = fs->make<TH2D>("jetSCTimevDrTime", "JetSCTime v JetDrTime;JetSCTime;JetDrTime", 280, -15, 25, 280, -15, 25 );
    hist2d[105] = fs->make<TH2D>("jetGenTimevSCJetTime", "GenTime v SCJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[106] = fs->make<TH2D>("jetGenTimevGenEnergy", "GenTime v GenEnergy;GenTime;GenEnergy", 280, -15, 25, 100, 0, 1000 );
    hist2d[107] = fs->make<TH2D>("jetGenjetDrvSCTimeDiff", "Jet GenJet Dr v SCJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[108] = fs->make<TH2D>("jetGenjetDrvDRTimeDiff", "Jet GenJet Dr v DRJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[109] = fs->make<TH2D>("jetGenTimevTOFcorr", "GenTime v TOFcorr;GenTime;TOFcorr", 250, 0, 25, 250, 0, 25 );
    hist2d[110] = fs->make<TH2D>("jetEGenERatiovSCTime", "Jet E/GenE v SCTime;E/GenE;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[111] = fs->make<TH2D>("jetEGenERatiovDRTime", "Jet E/GenE v DRTime;E/GenE;DRTime", 80, 0, 2, 40, -15, 25 );
    hist2d[112] = fs->make<TH2D>("jetGenVarvSCTimeDiff", "Jet GenJet Var v SCJet GenJet TimeDif;Var;TimeDif", 270, -2, 25, 200, 0, 20.0 );
    hist2d[113] = fs->make<TH2D>("jetGenPurityvSCTimeDiff", "Jet GenJet Purity v SCJet GenJet TimeDif;Purity;TimeDif", 100, 0, 1, 200, 0, 20.0 );
    hist2d[114] = fs->make<TH2D>("jetGenPurityvGenJetVar", "Jet GenJet Purity v GenJet Var;Purity;Var", 100, 0, 1, 270, -2, 25 );
    hist2d[115] = fs->make<TH2D>("jetGenVarvGenJetNKids", "Jet GenJet Var v GenJet nKids;Var;nKids", 270, -2, 25, 100, 0, 100 );
    hist2d[116] = fs->make<TH2D>("jetGenPurityvGenJetNKids", "Jet GenJet Purity v GenJet nKids;Purity;nKids", 100, 0, 1, 100, 0, 100 );
    hist2d[117] = fs->make<TH2D>("genJetSCTimeDiffvDrMatchJet", "genJet SCTimeDiff v DrMatchJet;SCTimeDiff;DrMatchJet", 300, 0, 30, 320, 0, 3.2 );
    hist2d[118] = fs->make<TH2D>("jetGenTimevJetEMFrac", "GenTime v JetEMFrac;GenTime;JetEMFrac", 280, -15, 25, 150, 0, 1.5 );

    hist2d[119] = fs->make<TH2D>("jetSlopevRotSlope", "Jet Slope v rotated Slope;Slope;rotated Slope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    hist2d[120] = fs->make<TH2D>("jetSlopevDifSlope", "Jet Slope v dif w/ rotSlope;Slope;difSlope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    hist2d[121] = fs->make<TH2D>("jetPhivSlope", "Jet Phi v Slope;Phi;Slope ps/cm", 140, -3.5, 3.5, sldiv, slmin, slmax);

    hist2d[122] = fs->make<TH2D>("ebRhTimevEnergy", "ebRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );
    hist2d[123] = fs->make<TH2D>("eeRhTimevEnergy", "eeRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );

    hist2d[124] = fs->make<TH2D>("jetEmFracvGenjetDr", "Jet emFrac v GenJet Dr;emFrac;jetGenDr", 40, 0, 1, 200, 0, 0.5 );
    hist2d[125] = fs->make<TH2D>("jetEGenERatiovGenjetDr", "Jet E/GenE v GenJet Dr;E/GenE;jetGenDr", 80, 0, 2, 200, 0, 0.5 );

    hist2d[126] = fs->make<TH2D>("jetEtavSlope3D", "Jet Eta v Slope 3D;Eta;Slope ps/cm", 60, -1.5, 1.5, sl3div, sl3min, sl3max);
    hist2d[127] = fs->make<TH2D>("jetPhivSlope3D", "Jet Phi v Slope 3D;Phi;Slope ps/cm", 140, -3.5, 3.5, sl3div, sl3min, sl3max);



 ------------ method called once each job just after ending the event loop	------------
void LLPgammaAnalyzer_AOD::endJob(){

    normTH1D(hist1d[15]);
    normTH1D(hist1d[16]);
    normTH1D(hist1d[39]);
    normTH1D(hist1d[40]);

	normTH1D(hist1d[48]);
    normTH1D(hist1d[49]);
    normTH1D(hist1d[41]);
    normTH1D(hist1d[42]);

    normTH1D(hist1d[27]);
    normTH1D(hist1d[28]);

    normTH1D(hist1d[59]);
    normTH1D(hist1d[60]);

	hist2d[73]->Divide(hist2d[74]);
    hist2d[75]->Divide(hist2d[76]);
    hist2d[77]->Divide(hist2d[78]);

}>>>>void LLPgammaAnalyzer_AOD::endJob()
@
