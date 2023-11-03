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
    f1 = []
    f2 = []
    f1b = []
    f2b = []
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

    mg = TMultiGraph();
    first = True;

    for histname, tree, infileb1, infileb2, infile1, infile2, lego, cutside in hist_list :
    
        f1b.append(TFile.Open(infileb1))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname

        f2b.append(TFile.Open(infileb2))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
		
        f1.append(TFile.Open(infile1))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
 
        f2.append(TFile.Open(infile2))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname

        orighist1 = f1[n].Get(hist)
        orighist2 = f2[n].Get(hist)
        orighist1b = f1b[n].Get(hist)
        orighist2b = f2b[n].Get(hist)
        htitle = 'hist' + str(n)
        lenmybins = int(orighist1.GetNbinsX())
        h1.append(TGraphErrors(lenmybins))
        norm1 = orighist1b.Integral()
        norm2 = orighist2b.Integral()
        if norm1 == 0 : norm1 = 1
        if norm2 == 0 : norm2 = 1
        print( "Norms : ", norm1, norm2 )

        for bn in range( 2, lenmybins):
            sigpass = 0
            bkgfail = 0
            rint1 = orighist1.Integral(bn,lenmybins)
            lint1 = orighist1.Integral(1,bn-1)
            rval1 = rint1/norm1
            lval1 = lint1/norm1
            rint2 = orighist2.Integral(bn,lenmybins)
            lint2 = orighist2.Integral(1,bn-1)
            rval2 = rint2/norm2
            lval2 = lint2/norm2
            #print( "`Intergrals : ", lint1, lint2 )
            if( cutside == "<" ):
                sigpass = lval1
                bkgpass = lval2
            else :
                sigpass = rval1
                bkgpass = rval2
            h1[n].SetPoint(bn,sigpass,bkgpass)
            #print( " set point : ", bn, sigpass, bkgpass,  )
        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+1,kGray+2,kGray+3,kBlack]
        h1[n].SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        msz = 0.8
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        legend.AddEntry(h1[n],lego+' '+cutside,'epl');
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
xtitle = 'Signal Eff.'
#xtitle = '#Delta_{Run}'
#xtitle = 'GeV'
#xtitle = '[ns]'
#ytitle = 'Ave Xtal Time [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#mu(t_{1}-t_{2}) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
#xtitle = 'HcalTowerSumEtBcConeDR04 [GeV]'
#ytitle = 'a.u.'
ytitle = 'Bkg Eff.'
htitle = ''
#islogx = True
islogx = False
#islogy = True
islogy = False

#---------------------------------------------------------------
#hl_mc_fms_loc = [
#     #['hist_name","tree_name",hist_file_location","legend_name"],
#     #["Data_sigma","",mc_full_loc+pcal+lstfr,"Full"],
#     #["Data_sigma","",mc_multi_loc+pcal+lstfr,"Multi"],
#     #["Data_sigma","",mc_single_loc+pcal+lstfr,"Single"],
#]

#l = [ 0.7,0.7,0.925,0.9 ] # legend position top right
l = [ 0.2,0.55,0.425,0.75 ] # legend position top left
#l = [ 0.25,0.20,0.52,0.525 ] # legend position bottom left
t = [0.2,0.84,0.0,0.175,0.225] # titles position

#rhname1 = "selPhoHcalTowerSumEtBcConeDR04"
#rhname1 = "pho30ptHcalTowerSumEtBcConeDR04"
#rhname1 = "pho100ptHTSumEtBcConeDR04"
rhname1l = [ "selPhoHcalTowerSumEtBcConeDR04", "pho30ptHcalTowerSumEtBcConeDR04", "pho100ptHTSumEtBcConeDR04" ]
cut1 = '<'
##rhname = "selPhoTrkSumPtHollowConeDR04" # same
#rhname2 = "selPhoTrkSumPtSolidConeDR04" # same
#rhname2 = "pho30ptTrkSumPtSolidConeDR04" # same
#rhname2 = "pho100ptTrkSumPtSolidConeDR04" # same
rhname2l = [ "selPhoTrkSumPtSolidConeDR04", "pho30ptTrkSumPtSolidConeDR04", "pho100ptTrkSumPtSolidConeDR04" ]
cut2 = '<'
#rhname3 = "selPhoEcalRHSumEtConeDR04"
#rhname3 = "pho30ptEcalRHSumEtConeDR04"
#rhname3 = "pho100ptEcalRHSumEtConeDR04"
rhname3l = [ "selPhoEcalRHSumEtConeDR04", "pho30ptEcalRHSumEtConeDR04", "pho100ptEcalRHSumEtConeDR04" ]
cut3 = '<'
#rhname4 = "selPhoHadTowOverEM"
#rhname4 = "pho30ptHadTowOverEM"
#rhname4 = "pho100ptHadTowOverEM"
rhname4l = [ "selPhoHadTowOverEM", "pho30ptHadTowOverEM", "pho100ptHadTowOverEM" ]
cut4 = '<'

rfgmsbrocj0 = "KUCMS_JetHt_18D_Met150_notSig_v15_iso0_Skim_BaseHists.root"
rfgmsbrocj1P = "KUCMS_JetHt_18D_Met150_notSig_v15_iso1_Skim_BaseHists.root"
rfgmsbrocj1a = "KUCMS_JetHt_18D_Met150_notSig_v15_iso2_Skim_BaseHists.root"
rfgmsbrocj1b = "KUCMS_JetHt_18D_Met150_notSig_v15_iso3_Skim_BaseHists.root"
rfgmsbrocj1c = "KUCMS_JetHt_18D_Met150_notSig_v15_iso4_Skim_BaseHists.root"

rfgmsbrocs0 = "KUCMS_GMSB_L100_Met150_Signal_v15_iso0_Skim_BaseHists.root"
rfgmsbrocs1P = "KUCMS_GMSB_L100_Met150_Signal_v15_iso1_Skim_BaseHists.root"
rfgmsbrocs1a = "KUCMS_GMSB_L100_Met150_Signal_v15_iso2_Skim_BaseHists.root"
rfgmsbrocs1b = "KUCMS_GMSB_L100_Met150_Signal_v15_iso3_Skim_BaseHists.root"
rfgmsbrocs1c = "KUCMS_GMSB_L100_Met150_Signal_v15_iso4_Skim_BaseHists.root"

rfgmsbrocg0 = "KUCMS_GMSB_L100_Met150_notSig_v15_iso0_Skim_BaseHists.root"
rfgmsbrocg1P = "KUCMS_GMSB_L100_Met150_notSig_v15_iso1_Skim_BaseHists.root"
rfgmsbrocg1a = "KUCMS_GMSB_L100_Met150_notSig_v15_iso2_Skim_BaseHists.root"
rfgmsbrocg1b = "KUCMS_GMSB_L100_Met150_notSig_v15_iso3_Skim_BaseHists.root"
rfgmsbrocg1c = "KUCMS_GMSB_L100_Met150_notSig_v15_iso4_Skim_BaseHists.root"

x = [ 0.0, 1.1 ]
y = [ 0.0, 1.1 ]

rfgmsbroc1  = rfgmsbrocs0
rfgmsbroc2g = rfgmsbrocg0
rfgmsbroc2j  = rfgmsbrocj0

title = ' 2018 GMSB L100 '
ptcut = 'Pt < '
#isocut = 'hadTowOverEM < 0.13'
isocut = ''
ptlist = [ '20','30','100']
outfileend = '_pt'

for rhname1, rhname2, rhname3, rhname4, pt in zip( rhname1l, rhname2l, rhname3l, rhname4l, ptlist ) :

  ptitle=[ title, ptcut+pt, isocut ]
  rfgmsbroc2  = rfgmsbroc2g
  bkgrd = 'GMSB '
  outname = 'llpa_met150_multiroc_gmsb'+outfileend+pt
  layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  
  inhistlist = [
  
  #    [ rhname1, "", base1-sig, base2-bkgrd, rfgmsbroc1-sig, rfgmsbroc2-bkgrd, rhname1, cut1, usebase1 ],
      [ rhname1, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname1, cut1 ],
      [ rhname2, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname2, cut2 ],
      [ rhname3, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname3, cut3 ],
      [ rhname4, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname4, cut4 ]
  
  ]
  
  dostack(inhistlist, outname, date, layout, ptitle,  y, x, l, t)
  
  ptitle=[ title, ptcut+pt, isocut ]
  rfgmsbroc2  = rfgmsbroc2j
  outname = 'llpa_met150_multiroc_jetht_'+outfileend+pt
  bkgrd = 'JetHt '
  layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  
  inhistlist = [
  
  #    [ rhname1, "", base1-sig, base2-bkgrd, rfgmsbroc1-sig, rfgmsbroc2-bkgrd, rhname1, cut1, usebase1 ],
      [ rhname1, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname1, cut1 ],
      [ rhname2, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname2, cut2 ],
      [ rhname3, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname3, cut3 ],
      [ rhname4, "", rfgmsbroc1, rfgmsbroc2, rfgmsbroc1, rfgmsbroc2, bkgrd+rhname4, cut4 ]
  
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


