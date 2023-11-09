#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack( hist_list, outname, date, layout, ptitle, y, x, l, t, ry ):

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
    oh = []
    n = 0

    setTDRStyle()
    #gStyle.SetPadTickY(1)
    #gStyle.SetPadTickX(1)
    #if layout['logx'] : gStyle.SetOptLogx(1)
    #if layout['logy'] : gStyle.SetOptLogy(1)
    canvas = TCanvas( 'c1', 'canvas' , 0, 0, 800, 700 )
    canvas.cd()
    #c2 = TPad("c2","",0.05,0.05,0.95,0.35)
    c1 = TPad("c1","", 0.05, 0.4, 0.95, 0.90)
    #if layout['logx'] : c1.SetLogx()
    #if layout['logy'] : c1.SetLogy()
    c1.SetTopMargin(0.05)
    c1.SetBottomMargin(0.1)
    c1.SetGridx(1)
    c1.SetGridy(1)
    if layout['logx'] : c1.SetLogx()
    if layout['logy'] : c1.SetLogy()
    #c1.Update()
    c1.Draw()
    c1.cd()
    setTDRStyle()

    legend = TLegend(l[0],l[1],l[2],l[3]);
#    legend.SetHeader(layout['legtitle'], '')
    legend.SetName('legend')
    gStyle.SetLegendFont(42)
    gStyle.SetLegendTextSize(0.05)
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

    #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
    k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
    #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
    #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
    #k = [kBlack]
    #k = [kGray+1,kGray+2,kGray+3,kBlack]

    for histname, tree, infile, lego  in hist_list :
    
        f1.append(TFile.Open(infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
 
        orighist = f1[n].Get(hist)
        oh.append(orighist)
        htitle = 'hist' + str(n)
        lenmybins = int(orighist.GetNbinsX())
        h1.append(TGraphErrors(lenmybins))
        norm = orighist.Integral()
        if norm == 0 : norm = 1

        for bn in range( 0, lenmybins + 1 ):
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
        h1[n].SetLineColor(k[n])
        if dofit : hfit.SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        msz = 1.0
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

    #canvas.cd(0)
    mg.Draw('AP')

    #mg.UseCurrentStyle()
    #mg.SetMarkerStyle(n+25)
    #mg.SetMarkerStyle(6)
    mg.SetTitle(layout['title'])
#+';X axis '+layout['xtitle']+';Y Axis '+layout['ytitle']+';')
    #mg.GetXaxis().CenterTitle(True)
    #mg.GetXaxis().SetTitle(layout['xtitle'])
    mg.GetXaxis().SetLabelSize(0)
    mg.GetXaxis().SetLabelOffset(999)
    #mg.GetYaxis().SetLabelSize(20)
    mg.GetYaxis().ChangeLabel(1, -1, -1, -1, -1, -1, " ")
    mg.GetYaxis().SetLabelFont(43)
    mg.GetYaxis().SetLabelSize(20)#20
    mg.GetYaxis().CenterTitle(True)
    mg.GetYaxis().SetTitleFont(43)
    mg.GetYaxis().SetTitleSize(20)#32
    mg.GetYaxis().SetTitleOffset(2.5)
    mg.GetYaxis().SetTitle(layout['ytitle'])
    mg.SetMinimum(y[0])
    mg.SetMaximum(y[1])
#   mg.GetXaxis().SetRangeUser(200.0,1100.0)
    mg.GetXaxis().SetRangeUser(x[0],x[1])
#    if layout['logx'] : mg.GetXaxis().SetMoreLogLabels()
#    if layout['logy'] : mg.GetYaxis().SetMoreLogLabels()

    if lego != 'none' : legend.Draw('same')  #   legend inclusion switch
    if layout['logx'] : canvas.cd(0).SetLogx()
    if layout['logy'] : canvas.cd(0).SetLogy()
    canvas.Modified()
    canvas.Update()

    #c1 = TPad("c1","",0.05,0.4,0.95,0.95)
    c2 = TPad("c2","", 0.05, 0.05, 0.95, 0.45)
    c2.SetTopMargin(0.1)
    c2.SetBottomMargin(0.3)
    c2.SetGridx(1)
    c2.SetGridy(1)
    #c2.Update()
    #gStyle.SetOptLogy(0)
    if layout['logx'] : c2.SetLogx()
    c2.Draw()
    c2.cd()
    setTDRStyle()
    #canvas.cd(1).cd()

    oh[0].Sumw2()
    noh0 = oh[0].Integral()
    oh[0].Scale(1/noh0)

    first = True
    lk = 0
    ohn = []
    for ohs in oh[1:] :

        ohn.append(ohs)
        ohn[lk].Sumw2()
        noh1 = ohn[lk].Integral()
        ohn[lk].Scale(1/noh1)
        #rh = ohn[lk].Clone("rh")
        #rh.Sumw2()
        ohn[lk].Divide(oh[0])
        ohn[lk].SetLineColor(k[lk+1])
        ohn[lk].SetMarkerColor(k[lk+1])
        ohn[lk].SetMarkerStyle(lk+26)
        msz = 1.0
        if( lk == 2 ) : ohn[lk].SetMarkerSize(msz+0.3)
        elif( lk == 3 ) : ohn[lk].SetMarkerSize(msz+0.5)
        else : ohn[lk].SetMarkerSize(msz)
        ohn[lk].SetMarkerSize(msz)

        ohn[lk].SetMinimum(ry[0])#;  // Define Y ..
        ohn[lk].SetMaximum(ry[1])#; // .. range
        ohn[lk].GetXaxis().SetRangeUser(x[0],x[1])
        ohn[lk].SetStats(0)#;      // No statistics on lower plot

        ohn[lk].GetXaxis().ChangeLabel(1, -1, -1, -1, -1, -1, " ")
        ohn[lk].GetXaxis().SetLabelFont(43)
        ohn[lk].GetXaxis().SetLabelSize(20)#20
        ohn[lk].GetYaxis().SetNdivisions(505)
        ohn[lk].GetYaxis().ChangeLabel(1, -1, -1, -1, -1, -1, " ")
        ohn[lk].GetYaxis().SetLabelFont(43)
        ohn[lk].GetYaxis().SetLabelSize(25)#20
        ohn[lk].GetYaxis().SetNdivisions(505)

        ohn[lk].GetXaxis().CenterTitle(True)
        ohn[lk].GetXaxis().SetTitleFont(43)
        ohn[lk].GetXaxis().SetTitleSize(20)
        ohn[lk].GetXaxis().SetTitleOffset(1.4)
        ohn[lk].GetXaxis().SetTitle(layout['xtitle'])
        ohn[lk].GetYaxis().CenterTitle(True)
        ohn[lk].GetYaxis().SetTitleSize(20)#32
        ohn[lk].GetYaxis().SetTitleFont(43)
        ohn[lk].GetYaxis().SetTitleOffset(2.5)
        ohn[lk].GetYaxis().SetTitle("Ratio [CC/Rt]")

        if first : 
            ohn[lk].Draw("ep")#;       // Draw the ratio plot
            first = False
        else :
            ohn[lk].Draw("ep same")

        lk += 1 
        canvas.Modified()
        canvas.Update()

    canvas.cd()
    #lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    lat_cms = '#bf{CMS} #it{WkInPgrs}' + ptitle[0]
    #lat_title = ptitle[1]+' (13 TeV)'
    lat_title = ptitle[1]
    #lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + 2C^{2}'
    lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + #frac{S^{2}}{A_{eff}/#sigma_{n}} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (N/Eeff)^{2} + 2C^{2}'
    lat.SetTextSize(0.035);
    lat.SetTextFont(42);
    lat.DrawLatex(0.1,0.9325,lat_cms);
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

    canvas.Modified()
    canvas.Update()

    #canvas.Print( outname + '_' + date + '.pdf' )
    canvas.Print( outname + '_' + date + '.png' )
    #canvas.SaveAs( outname + '_' + date + '.root' )
    #canvas.SaveAs( outname + '_' + date + '.C' )
    #canvas.Show()
    canvas.Close()

#from overlay_hist_defs_v3 import *

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

legtitle = ''
#legtitle = 'KuStc'
#legtitle = 'KuNotStc'

rtitle = 'Run3 DQMReco' 
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
#ytitle = 'a.u.'
ytitle = 'N_{bin}/N_{tot}'
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

dqm_file1 = "DQM_V0001_R000369978__Global__CMSSW_X_Y_Z__RECO_EGamma0.root"
dqm_file2 = "DQM_V0001_R000369978__Global__CMSSW_X_Y_Z__RECO_EGamma0_newparam.root"
dqm_file3 = "DQM_V0001_R000369978__Global__CMSSW_X_Y_Z__RECO_EGamma0_ratio.root"
dqm_file4 = "DQM_V0001_R000369978__Global__CMSSW_X_Y_Z__RECO_EGamma0_ratiomod.root"

l = [ 0.7,0.7,0.925,0.9 ] # legend position top right
#l = [ 0.2,0.65,0.425,0.9 ] # legend position top left
#l = [0.25,0.20,0.52,0.525 ] # bottom left
#l =  [0.7,0.20,0.925,0.525 ] # bottom right

t = [0.2,0.825,0.0,0.175,0.225] # titles position

#rhtree = "DQMData/Run 369978/EcalBarrel/Run summary/EBClusterTask/"
#rhname = "EBCLT SC number"
#rhfilename = "EBCLT_SC_number"
#x = [ 0.0, 20.0 ]
#y = [ 0.0001, 1 ]
#ry = [ 0.5, 1.5 ]
#rhname = "EBCLT SC seed crystal energy"
#rhfilename = "EBCLT_SC_seed_crystal_energy"
#xtitle = "EBCLT SC seed crystal Energy [GeV]"
#x = [ 0.0, 250.0 ]
#y = [ 0.0001, 1 ]
#ry = [ 0.0, 1.5 ]
#rhtree = "DQMData/Run 369978/EcalBarrel/Run summary/EBRecoSummary/"
#rhname = "recHits_EB_recoFlag"
#rhfilename = "recHits_EB_recoFlag"
#xtitle = "recHits EB recoFlags"
#x = [ 0.0, 15.0 ]
#y = [ 0.000001, 1.1 ]
#ry = [ 0.0, 2.0 ]
#rhtree = "DQMData/Run 369978/Egamma/Run summary/Electrons/Ele5_TagAndProbe/"
#rhname = "ele201_mee_os"
#rhfilename = "ele5_tp_ele201_mee_os"
#xtitle = "Ele5 Tag&Probe : ele201_mee_os [GeV]"
#x = [ 50.0, 130.0 ]
#y = [ 0.001, 1.0 ]
#ry = [ 0.5, 2.0 ]
#rhtree = "DQMData/Run 369978/Egamma/Run summary/Electrons/Ele2_All/"
#rhname = "ele201_mee_os"
#rhfilename = "ele2_all_ele201_mee_os"
#xtitle = "Ele2 All : ele201_mee_os [GeV]" 
#x = [ 0.0, 150.0 ]
#y = [ 0.0001, 1.0 ]
#ry = [ 0.4, 1.6 ]
#rhtree = "DQMData/Run 369978/Egamma/Run summary/gedPhotonAnalyzer/AllPhotons/Et above 20 GeV/"
#rhname = "h_01_phoEBarrel"
#rhfilename = "gedPho_all_h_01_phoEBarrel_energy"
#xtitle = "gedPhoton All : h_01_phoEBarrel Energy [GeV]"
#x = [ 0.0, 500.0 ]
#y = [ 0.00001, 1.0 ]
#ry = [ 0.25, 1.75 ]
#rhname = "h_10_r9Barrel"
#rhfilename = "gedPho_all_h_10_r9Barrel"
#xtitle = "gedPhoton All : h_10_r9Barrel R9"
#x = [ 0.0, 1.2 ]
#y = [ 0.00000001, 1.0 ]
#ry = [ 0.0, 6.0 ]
#rhtree = "DQMData/Run 369978/Egamma/Run summary/PiZeroAnalyzer/"
#rhname = "Pt1Pi0EB"
#rhfilename = "egamma_pizero_Pt1Pi0EB"
#xtitle = "PiZero Analyzer : 1st photon Pt [GeV]"
#x = [ 0.0, 20 ]
#y = [ 0.0001, 1.0 ]
#ry = [ 0.75, 1.25 ]
#rhtree = "DQMData/Run 369978/Egamma/Run summary/gedPhotonAnalyzer/GoodCandidatePhotons/Et above 20 GeV/"
#rhname = "h_05_nPhoBarrel"
#rhfilename = "gedPhoton_goodcandidate_h_05_nPhoBarrel"
#xtitle = "gedPhoton Good Canadidte : # Photon EB"
#x = [ 0.0, 6 ]
#y = [ 0.0000001, 1.5 ]
#ry = [ -0.5, 2.0 ]
#rhtree = "DQMData/Run 369978/EcalBarrel/Run summary/EBTimingTask/"
#rhname = "EBTMT timing 1D summary"
rhtree = "DQMData/Run 369978/EcalBarrel/Run summary/EBSummaryClient/"
#rhname = "EBTMT timing mean 1D summary"
#rhfilename = "EBTMT_timing_mean_1D_summary"
#xtitle = "EBSummaryClient : EBTMT timing mean 1D summary"
#x = [ -20.0, 20.0 ]
#y = [ 0.000001, 1.5 ]
#ry = [ -1, 5 ]
rhname = "EBTMT timing rms 1D summary"
rhfilename = "EBTMT_timing_rms_1D_summary"
xtitle = "EBSummaryClient : EBTMT timing rms 1D summary"
x = [ 0.0, 4 ]
y = [ 0.00001, 1.5 ]
ry = [ -0.5, 100 ]

outname = 'dqm_v1_' + rhfilename
layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
ptitle=[' 2023D EGamma0 ','','']

inhistlist = [

    [ rhname, rhtree, dqm_file3, "Rt"],
    [ rhname, rhtree, dqm_file1, "Orig CC"],
    [ rhname, rhtree, dqm_file2, "New Params CC"],

]

dostack(inhistlist, outname, date, layout, ptitle,  y, x, l, t, ry)

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


