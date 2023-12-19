#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack( hist_list, outname, date, layout, ptitle, y, x, l, t ):

    first = True
    #dofit = True
    dofit = False
    #sxtal = True
    sxtal = False
    paramn = []
    parnerror = []
    paramc = []
    parcerror = []
    params = []
    parserror = []
    thebinmid = []
    thebinerror = []
    f1 = []
    h1 = []
    n = 0

    setTDRStyle()
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
    if layout['logx'] : gStyle.SetOptLogx(1)
    if layout['logy'] : gStyle.SetOptLogy(1)
    c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
    c1.cd()
    #if layout['logx'] : c1.SetLogx()
    #if layout['logy'] : c1.SetLogy()
    c1.SetGridx(1)
    c1.SetGridy(1)
    c1.Update()
    c1.Draw()

    legend = TLegend(l[0],l[1],l[2],l[3]);
#    legend.SetHeader(layout['legtitle'], '')
    legend.SetName('legend')
    gStyle.SetLegendFont(42)
    gStyle.SetLegendTextSize(0.03)
    legend.SetBorderSize(0)
    #legend.SetLineColor(kBlack)

    lat = TLatex() 
    lat.SetNDC()

    hMax = y[0]
    hMin = y[1]

    if dofit :
        ns=str(n)
        hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',100,750,3)
        hfit.SetParName(0,'N')
        hfit.SetParameter(0,40.0)
        #hfit.SetParLimits(0,0,50)
        hfit.SetParLimits(0,0.0,100.0)
        #hfit.SetParameter(0,5)
        #hfit.SetParLimits(0,0,10)
        hfit.SetParName(1,'C')
        hfit.SetParameter(1,0.1)
        hfit.SetParLimits(1,0.0,1.0)
        #hfit.SetParLimits(1,0.01,10.0)
        #hfit.SetParLimits(1,0.02,1.0)
        #hfit.SetParameter(1,0.05)
        #hfit.SetParLimits(1,0.001,1.0)
        hfit.SetParName(2,'S')
        hfit.SetParameter(2,5.0)
        hfit.SetParLimits(2,0.0,25.0)

    mg = TMultiGraph();

    for histname, tree, infile, lego  in hist_list :
    
        f1.append(TFile.Open(infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
 
        orighist = f1[n].Get(hist)
        htitle = 'hist' + str(n)
        lenmybins = int(orighist.GetNbinsX())
        h1.append(TGraphErrors(lenmybins))
        norm = orighist.Integral()
        if norm == 0 : norm = 1

        for bn in range( 1,lenmybins-1):
            obinval = float(orighist.GetBinContent(bn))
            binval = obinval/norm
            binerr = float(orighist.GetBinError(bn))/norm
            binmid = float(orighist.GetBinCenter(bn))
            binwidth = float(orighist.GetBinWidth(bn))
            #binstart = float(orighist.GetBinLowEdge(bn))
            h1[n].SetPoint(bn,binmid,binval)
            if obinval == 0 : obinval = 1;
            widtherr = binwidth/(2*sqrt(obinval))
            h1[n].SetPointError(bn,widtherr,binerr)
            #print('Fill bin',bn,'at',binmid,'with',binval,'err',widtherr,'by',binerr,'for:',binstart,'to',binstart+binwidth,'width',binwidth) 
        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+1,kGray+2,kGray+3,kBlack]
        h1[n].SetLineColor(k[n])
        if dofit : hfit.SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        msz = 0.8
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        if( first ) :
                if dofit : h1[n].Fit(hfit,'RE')
                #if dofit : h1[n].Fit(hfit,'REQ')
                first = False
        else :
                if dofit : h1[n].Fit(hfit,'RE+')
                #if dofit : h1[n].Fit(hfit,'REQ+') 

        if dofit : 
                 paramn.append(str(abs(hfit.GetParameter(0))))
                 paramc.append(str(abs(hfit.GetParameter(1))))
                 params.append(str(abs(hfit.GetParameter(2))))
                 pne = hfit.GetParError(0)
                 pce = hfit.GetParError(1)
                 pse = hfit.GetParError(2)
                 #print('Fit info',paramn[n],pne,paramc[n],pce)
                 #print('Fit info',paramn[n],pne,paramc[n],pce,params[n],pse)
                 if pne < 0.01 : pne = 0.01
                 if pce < 0.0001 : pce = 0.0001
                 parnerror.append(str(pne))
                 parcerror.append(str(pce))
                 parserror.append(str(pse))
 
        legend.AddEntry(h1[n],lego,'epl');
        n += 1

        #End of loop
    
    for h in range(0,n):
        mg.Add(h1[h])

    mg.Draw('AP')

    #mg.UseCurrentStyle()
    #mg.SetMarkerStyle(n+25)
    #mg.SetMarkerStyle(6)
    mg.SetTitle(layout['title'])
#+';X axis '+layout['xtitle']+';Y Axis '+layout['ytitle']+';')
    mg.GetXaxis().CenterTitle(True)
    mg.GetXaxis().SetTitle(layout['xtitle'])
    mg.GetYaxis().CenterTitle(True)
    mg.GetYaxis().SetTitle(layout['ytitle'])
    mg.SetMinimum(y[0])
    mg.SetMaximum(y[1])
#   mg.GetXaxis().SetRangeUser(200.0,1100.0)
    mg.GetXaxis().SetRangeUser(x[0],x[1])
#    if layout['logx'] : mg.GetXaxis().SetMoreLogLabels()
#    if layout['logy'] : mg.GetYaxis().SetMoreLogLabels()

    if lego != 'none' : legend.Draw('same')  #   legend inclusion switch
    gPad.Modified()

    #lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    lat_cms = '#bf{CMS} #it{WkInPgrs}' + ptitle[0]
    #lat_title = ptitle[1]+' (13 TeV)'
    lat_title = ptitle[1]
    #lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + 2C^{2}'
    lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + #frac{S^{2}}{A_{eff}/#sigma_{n}} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (N/Eeff)^{2} + 2C^{2}'
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    lat.DrawLatex(0.15,0.9325,lat_cms);
    lat.DrawLatex((0.82-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.04);
    lat.DrawLatex(t[0],t[1],ptitle[2]);
    if dofit : 
        lat.SetTextSize(0.03);
        lat.DrawLatex(t[3],t[4]+.075,lat_form);
        for l in range(0,n):
            lat_param =	'#color['+str(k[l])+']{'
            lat_param = lat_param + 'N : '+paramn[l][0:6]+' #pm '+parnerror[l][0:6]+' [ns] '
            lat_param = lat_param + 'S : '+params[l][0:6]+' #pm '+parserror[l][0:6]+' [ns] '
            lat_param = lat_param + 'C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' [ns]}'
            lat.SetTextSize(0.03);
            lat.DrawLatex(t[3],t[4]-l*.035,lat_param);

    
    if layout['logx'] : c1.SetLogx()
    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.SaveAs( outname + '_' + date + '.C' )
    #c1.Show()
    c1.Close()

from overlay_hist_defs_v3 import *

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

legtitle = ''
#legtitle = 'KuStc'
#legtitle = 'KuNotStc'

rtitle = 'Run2 AODSIM'
#Ic_legtitle = ''
xtitle = ''
#xtitle = '#Delta_{Run}'
#xtitle = 'GeV'
#xtitle = '[ns]'
#ytitle = 'Ave Xtal Time [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#mu(t_{1}-t_{2}) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
#xtitle = 'HcalTowerSumEtBcConeDR04 [GeV]'
ytitle = 'a.u.'
htitle = ''
#islogx = True
islogx = False
islogy = True
#islogy = False

#---------------------------------------------------------------
#hl_mc_fms_loc = [
#     #['hist_name","tree_name",hist_file_location","legend_name"],
#     #["Data_sigma","",mc_full_loc+pcal+lstfr,"Full"],
#     #["Data_sigma","",mc_multi_loc+pcal+lstfr,"Multi"],
#     #["Data_sigma","",mc_single_loc+pcal+lstfr,"Single"],
#]

#rfname1 = "KUCMS_GMSB_L100t400_Skim_BaseHists.root"
#rfname2 = "KUCMS_GMSB_L100t400_gNino_Skim_BaseHists.root"
#rfname3 = "KUCMS_GMSB_L100t400_wzXino_Skim_BaseHists.root"
#rfname4 = "KUCMS_GMSB_L100t400_NotSUS_Skim_BaseHists.root"

#rfname1 = "KUCMS_GMSB_L100t400_Met150_Skim_BaseHists.root"
#rfname2 = "KUCMS_GMSB_L100t400_Met150_gNino_Skim_BaseHists.root"
#rfname3 = "KUCMS_GMSB_L100t400_Met150_wzXino_Skim_BaseHists.root"
#rfname4 = "KUCMS_GMSB_L100t400_Met150_NotSUS_Skim_BaseHists.root"

#rfname1 = "KUCMS_GMSB_L100t400_Met150_pho3_Skim_BaseHists.root"
#rfname2 = "KUCMS_GMSB_L100t400_Met150_pho3_gNino2_Skim_BaseHists.root"
#rfname3 = "KUCMS_GMSB_L100t400_Met150_pho3_wzXino_Skim_BaseHists.root"
#rfname4 = "KUCMS_GMSB_L100t400_Met150_pho3_NotSUS_Skim_BaseHists.root"

#rfname1 = "KUCMS_GMSB_L100t400_Met150_pho1m_gn_v6_Skim_BaseHists.root"
#rfname2 = "KUCMS_GMSB_L100t400_Met150_pho1m_not_v6_Skim_BaseHists.root"
#rfname3 = "KUCMS_GMSB_L100t400_Met150_pho1m_prt_v6_Skim_BaseHists.root"
#rfname4 = "KUCMS_GMSB_L100t400_Met150_pho1m_qgsr_v6_Skim_BaseHists.root"
#rfname5 = "KUCMS_GMSB_L100t400_Met150_pho1m_v6_Skim_BaseHists.root"
#rfname6 = "KUCMS_GMSB_L100t400_Met150_pho1m_wz_v6_Skim_BaseHists.root"
#rfname7 = "KUCMS_GMSB_L100t400_Met150_pho1m_xsr_v6_Skim_BaseHists.root"

#rfname1 = "KUCMS_GMSB_L100t400_Met150_jet1m_gn_v6_Skim_BaseHists.root"
#rfname2 = "KUCMS_GMSB_L100t400_Met150_jet1m_not_v6_Skim_BaseHists.root"
#rfname3 = "KUCMS_GMSB_L100t400_Met150_jet1m_prt_v6_Skim_BaseHists.root"
#rfname4 = "KUCMS_GMSB_L100t400_Met150_jet1m_qgsr_v6_Skim_BaseHists.root"
#rfname5 = "KUCMS_GMSB_L100t400_Met150_jet1m_v6_Skim_BaseHists.root"
#rfname6 = "KUCMS_GMSB_L100t400_Met150_jet1m_wz_v6_Skim_BaseHists.root"
#rfname7 = "KUCMS_GMSB_L100t400_Met150_jet1m_xsr_v6_Skim_BaseHists.root"

rfgmsb100 = "KUCMS_GMSB_L100_Met150_Other_v14_Skim_BaseHists.root"

rfgmsb1 = "hist_files/KUCMS_GMSB_L100t400_Met150_Signal_v10_Skim_BaseHists.root"
rfgmsb2 = "hist_files/KUCMS_GMSB_L100t400_Met150_XFSR_v10_Skim_BaseHists.root"
rfgmsb3 = "hist_files/KUCMS_GMSB_L100t400_Met150_Other_v10_Skim_BaseHists.root"
rfgmsb4 = "hist_files/KUCMS_GMSB_L100t400_Met150_UnMatched_v10_Skim_BaseHists.root"
#rfgmsb3 = "KUCMS_GMSB_L200_Met0_UnMatched_v9A_Skim_BaseHists.root"

rfgmsbsv1 = "hist_files/KUCMS_ShapeVar_Sig_SkimHists_v1.root"
rfgmsbsv2 = "hist_files/KUCMS_ShapeVar_GMSB_Bkg_SkimHists.root"
rfgmsbsv3 = "hist_files/KUCMS_ShapeVar_GJets_Bkg_SkimHists.root"

rfgjets1 = "hist_files/KUCMS_GJets_HT40tInf_Met150_Other_v10_Skim_BaseHists.root"
rfgjets2 = "hist_files/KUCMS_GJets_HT40tInf_Met150_UnMatched_v10_Skim_BaseHists.root"

rfgmsbroc1 = "hist_files/KUCMS_GMSB_L100t400_Met150_Signal_v11_Skim_BaseHists.root"
rfgmsbroc2 = "hist_files/KUCMS_GMSB_L100t400_Met150_notSig_v11a_Skim_BaseHists.root"

#y = [ 0.01, 1000000 ]
y = [ 0.00001, 1.0 ]
#y = [ 0.001, 0.1 ]
#y = [ 0.0, 0.035 ]
#y = [ 0.00001, 0.25 ]
#y = [ 0.0001, 0.4 ]
#y = [ 0.0, 0.01 ]
#y = [ 0.0, 0.005 ]
#x = [ 0.0, 0.001 ]
#x = [ 0.0, 0.05 ]
#x = [ 0.0, 0.1 ]
x = [ 0.0, 1.0 ]
#x = [ 0.0, 5.0 ]
#x = [ 0.0, 10.0 ]
#x = [ 0.0, 50.0 ]
#x = [ 0.0, 100.0 ]
#x = [ 0.0, 1500.0 ]
#x = [ -3.0, 3.0 ]
#x = [ -10.0, 25.0 ]
l = [ 0.7,0.7,0.925,0.9 ] # legend position top right
#l = [ 0.2,0.65,0.425,0.9 ] # legend position top left
t = [0.2,0.825,0.0,0.175,0.225] # titles position

#rhname = "evtMetPt"
#x = [ 0.0, 2500.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]

#rhname = "phoClstrRn"
#x = [ 0, 1.0 ]
#y = [ 0.0001, 0.15 ]
#l = [ 0.75,0.7,0.925,0.9 ]
#rhname = "phoEta"
#x = [ -3.0, 3.0 ]
#y = [ 0.00001, 0.0325 ]
#rhname = "phoPt"
#rhname = "phoEnergy" # same
#x = [ 0.0, 1500.0 ]
#x = [ 0.0, 250.0 ] 
#x = [ 0.0, 100.0 ] 
#y = [ 0.001, 0.1 ]
#rhname = "phoNrh"
#x = [ 0.0, 100.0 ]
#y = [ 0.00001, 0.25 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "phoSMin"
#x = [ 0, 1.0 ]
#y = [ 0.00001, 0.25 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "phoSieie"
#x = [ 0.0, 0.1 ]
#y = [ 0.00001, 0.25 ]
#l = [ 0.75,0.7,0.925,0.9 ]
#rhname = "phoTime"
#x = [ -10.0, 25.0 ]
#y = [ 0.00001, 0.25 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoHcalTowerSumEtBcConeDR04"
#rhname = "selPhoTrkSumPtHollowConeDR04" # same
#rhname = "selPhoTrkSumPtSolidConeDR04" # same
#x = [ 0.0, 5.0 ]
#x = [ 0.0, 25.0 ]
#y = [ 0.00001, 0.035 ]
#y = [ 0.0001, 0.5 ]
#y = [ 0.0001, 1 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#l = [ 0.7,0.25,0.925,0.5 ]
#rhname = "selPhoEcalRHSumEtConeDR04"
#x = [ 0.0, 25.0 ]
#y = [ 0.0001, 0.1 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoHadTowOverEM"
#x = [ 0.0, 0.05 ]
#x = [ 0.0, 0.25 ]
#y = [ 0.0001, 0.1 ]
#y = [ 0.0001, 0.5 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoSieip"
#x = [ 0.0, 0.001 ]
#l = [ 0.7,0.65,0.925,0.9 ]
##rhname = "selPhoSipip"
##x = [ 0.0, 0.1 ]

#rhname = "selPhoNrh_pt20"
#rhname = "selPhoNrh_pt30"
#rhname = "selPhoNrh_pt100"
#x = [ 0.0, 50.0 ]
#y = [ 0.0001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoEtaWidth_pt20"
#rhname = "selPhoEtaWidth_pt30"
#rhname = "selPhoEtaWidth_pt100"
#x = [ 0.0, 0.035 ]
#y = [ 0.0001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoPhiWidth_pt20"
#rhname = "selPhoPhiWidth_pt30"
#rhname = "selPhoPhiWidth_pt100"
#x = [ 0.0, 0.15 ]
#y = [ 0.0001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoR9_pt20"
#rhname = "selPhoR9_pt30"
#rhname = "selPhoR9_pt100"
#x = [ 0.2, 1.2 ]
#y = [ 0.0001, 1.0 ]
#l = [ 0.2,0.65,0.425,0.9 ] 
#rhname = "selPhoS4_pt20"
#rhname = "selPhoS4_pt30"
#rhname = "selPhoS4_pt100"
#x = [ 0.2, 1.2 ]
#y = [ 0.0001, 1.0 ]
#l = [ 0.2,0.65,0.425,0.9 ]
#rhname = "selPhoSAlp_pt20"
#rhname = "selPhoSAlp_pt30"
#rhname = "selPhoSAlp_pt100"
#x = [ -1.75, 1.75 ]
#y = [ 0.001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoSMaj_pt20"
#rhname = "selPhoSMaj_pt30"
#rhname = "selPhoSMaj_pt100"
#x = [ 0, 3.0 ]
#y = [ 0.001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoSMin_pt20"
#rhname = "selPhoSMin_pt30"
#rhname = "selPhoSMin_pt100"
#x = [ 0, 1.5 ]
#y = [ 0.0001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoSieie_pt20"
#rhname = "selPhoSieie_pt30"
rhname = "selPhoSieie_pt100"
#x = [ 0, 0.025 ]
#y = [ 0.0001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoSieip_pt20"
#rhname = "selPhoSieip_pt30"
#rhname = "selPhoSieip_pt100"
#x = [ -0.0005, 0.0005 ]
#y = [ 0.001, 1.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "selPhoSipip_pt20"
#rhname = "selPhoSipip_pt30"
rhname = "selPhoSipip_pt100"
x = [ 0, 0.025 ]
y = [ 0.0001, 1.0 ]
l = [ 0.7,0.65,0.925,0.9 ]

#y = [ 1, 1e8 ]
#rhname = "jetPt"
#rhname = "jetEnergy"
#rhname = "jetMass"
#x = [ 0.0, 500.0 ] 
#x = [ 0.0, 750.0 ]
#x = [ 0.0, 1500.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "jetEta"
#rhname = "jetPhi"
#x = [ -3.0, 3.0 ]
#rhname = "jetTime"
#x = [ -25.0, 25.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#rhname = "jetArea"
#rhname = "jetChEmEF"
#rhname = "jetChHM"
#rhname = "jetMuEF"
#rhname = "jetNeEmEF"
#rhname = "jetNeHEF"
#rhname = "jetNeHM"
#rhname = "jetChHEF"
#x = [ 0.0, 1500.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]

fhname = rhname
xtitle = fhname
outname = 'llpa_met150_v1_' + fhname
layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
ptitle=[' 2018 GMSB L100 ','','']

inhistlist = [
#    [ rhname, "", rfname5, "All Matched"],   
#    [ rhname, "", rfname1, "X0 -> Photon"],
#    [ rhname, "", rfname1, "SQ/Glino -> q->jet"],
#    [ rhname, "", rfname6, "X -> W/Z"],  
#    [ rhname, "", rfname4, "Gluino/Squark -->"],
#    [ rhname, "", rfname7, "X -->"], 
#    [ rhname, "", rfname3, "Not Susy"], 
#    [ rhname, "", rfname2, "Not Matched"],

#    [ rhname, "", rfgmsb1, "Signal"],
#    [ rhname, "", rfgmsb2, "XFSR"],
#    [ rhname, "", rfgmsb3, "Other"],
#    [ rhname, "", rfgmsb4, "UnMatched"],

#    [ rhname, "", rfgjets1, "Gjets"],
#    [ rhname, "", rfgjets2, "Gjets_UnMatched"],

#    [ rhname, "", rfgmsbroc1, "Signal"],
#    [ rhname, "", rfgmsbroc2, "NotSignal"],

#    [ rhname, "", rfgmsb100, "Other"],

    [ rhname, "", rfgmsbsv1, "Signal"],
    [ rhname, "", rfgmsbsv2, "notSignal"],
#    [ rhname, "", rfgmsbsv3, "GJets"],

]

dostack(inhistlist, outname, date, layout, ptitle,  y, x, l, t)

#ptitle=[' 2022 IOV5 359421-360089','','#splitline{EBEB}{CC Ave RH Time by Channel}'] #{GT 106X_dataRun2_v28}'
#y = [ 4.5, 0.5 ]
#x = [ 200.0, 700.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#t = [0.2,0.825,0.0,0.175,0.225]
#outname = 'downloads/tr_hl_r3_iov5cali_v7'
#dostack(hl_r3_iov5cali_v7, outname, date, Ic_layout, ptitle,  y, x, l, t)

#    legend = TLegend(0.25,0.20,0.52,0.525); # bottom left
#    legend = TLegend(0.4,0.205,0.6,0.525);   # bottom middle
#    legend = TLegend(0.4,0.60,0.6,0.90);   # top middle
#    legend = TLegend(0.645,0.50,0.825,0.9);   # top mid right
#    legend = TLegend(0.605,0.50,0.945,0.9);   # top right very wide
#    legend = TLegend(0.705,0.50,0.945,0.9);   # top right wide 
#    legend = TLegend(0.745,0.50,0.925,0.9);   # top right
#    legend = TLegend(0.745,0.40,0.925,0.9);   # top right tall
#    legend = TLegend(0.650,0.375,0.925,0.875);   # top mid right wide
#    legend = TLegend(0.62,0.60,0.8,0.9);   # top right
#    legend = TLegend(0.65,0.60,0.9,0.90);   # top right large

#      g_2->GetYaxis()->SetMoreLogLabels();
#      g_2->GetYaxis()->SetNoExponent();


