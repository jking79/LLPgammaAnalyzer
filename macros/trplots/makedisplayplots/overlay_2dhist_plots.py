#  jack w king 3

from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack( hist_list, outname, date, layout, ptitle, z, y, x, l, t ):

    first = True
    #dofit = True
    dofit = False
    f1 = []
    h1 = []
    n = 0

    setTDRStyle()
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
    gStyle.SetPadRightMargin(0.15)
    gStyle.SetTitleYOffset(0.8)
    #if layout['logx'] : gStyle.SetOptLogx(1)
    #if layout['logy'] : gStyle.SetOptLogy(1)
    c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
    #c1 = TCanvas( 'c1', 'canvas' , 200, 10, 1000, 1000 )
    c1.cd()
    if layout['logx'] : c1.SetLogx()
    if layout['logy'] : c1.SetLogy()
    if layout['logz'] : c1.SetLogz()
    c1.SetGridx(1)
    c1.SetGridy(1)
    c1.Update()
    c1.Draw()

#    legend = TLegend(l[0],l[1],l[2],l[3]);
#    legend.SetHeader(layout['legtitle'], '')
#    legend.SetName('legend')
#    gStyle.SetLegendFont(42)
#    gStyle.SetLegendTextSize(0.04)
#    legend.SetBorderSize(0)
    #legend.SetLineColor(kBlack)

    lat = TLatex() 
    lat.SetNDC()

    hMax = y[0]
    hMin = y[1]

    if dofit :
        ns=str(n)
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',75,2250,2)
        hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',6,100,2)
        hfit.SetParName(0,'N')
        #hfit.SetParameter(0,20)
        #hfit.SetParLimits(0,10,50)
        hfit.SetParameter(0,5)
        hfit.SetParLimits(0,0,10)
        hfit.SetParName(1,'C')
        #hfit.SetParameter(1,0.1)
        #hfit.SetParLimits(1,0.02,0.3)
        hfit.SetParameter(1,0.05)
        hfit.SetParLimits(1,0.001,1.0)

    for histname, tree, infile, lego  in hist_list :
    
        f1.append(TFile.Open(infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
        h1.append(f1[n].Get(hist))
        print( '____________________________________________________________________' )
        print( 'Infile Name : ', infile ) 
        print( 'Address : ', f1[n] )
        print( 'Hist Name : ', hist )
        print( 'Address : ', h1[n] )
        print( '____________________________________________________________________' )

#       print( f1 )
#       print( h1 )

        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #h1[n].SetMarkerStyle(6)
        h1[n].SetTitle(layout['title'])
        h1[n].GetXaxis().CenterTitle(True)
        h1[n].GetXaxis().SetTitle(layout['xtitle'])
        h1[n].GetYaxis().CenterTitle(True)
        h1[n].GetYaxis().SetTitle(layout['ytitle'])
#       k = [kMagenta+2,kBlue+1,kAzure+4,kBlack,kYellow+1,kViolet+7,kOrange+7,kRed+2,kGreen+3, kGray]
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        #k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kBlack]
#       k = [kMagenta+2,kBlue+2,kGreen+2]
        h1[n].SetLineColor(k[n])
        if dofit : hfit.SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        #msz = 0.2
        msz = 0.4
        #msz = 0.8
        #msz = 1.2
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        h1[n].SetMinimum(z[0])
        h1[n].SetMaximum(z[1])
#       h1[n].GetXaxis().SetRangeUser(200.0,1100.0)
        h1[n].GetXaxis().SetRangeUser(x[0],x[1])
        h1[n].GetYaxis().SetRangeUser(y[0],y[1])
        if layout['logx'] : h1[n].GetXaxis().SetMoreLogLabels()
        if layout['logy'] : h1[n].GetYaxis().SetMoreLogLabels()
        if( first ) :
                h1[n].Draw('colz')
                #h1[n].Draw('hist')
                if dofit : h1[n].Fit(hfit,'REQ')
                first = False
        else :
                h1[n].Draw('epSAME')
                if dofit : h1[n].Fit(hfit,'REQ+') 
        #c1.Update()
        #c1.Draw()

        if dofit : 
                 paramn = str(hfit.GetParameter(0))
                 paramc = str(hfit.GetParameter(1))
                 pne = hfit.GetParError(0)
                 pce = hfit.GetParError(1)
                 print(paramn,pne,paramc,pce)
                 if pne < 0.01 : pne = 0.01
                 if pce < 0.0001 : pce = 0.0001
                 parnerror = str(pne)
                 parcerror = str(pce)
                 #print( 'C: ' + param + ' +/- ' + parerror )
                 lat_param = '#color['+str(k[n])+']{N : ' + paramn[0:4] + ' #pm ' + parnerror[0:4] + ' [ns]  C : ' + paramc[0:6] + ' #pm ' + parcerror[0:6] + ' [ns]}'
                 lat.SetTextSize(0.03);
                 lat.DrawLatex(t[3],t[4]-n*.035,lat_param);    
        
#        legend.AddEntry(h1[n],lego,'epl');
        c1.Update()
        c1.Draw()
        n += 1

        #End of loop

    #legend.Draw('same')
    #c1.SetGridx(1)
    #c1.SetGridy(1)
    c1.cd()
    c1.Update()
    c1.Draw()

    lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    #lat_cms = '#bf{CMS} #it{Work in Progress}' + ptitle[0]
    #lat_title = 'Run2018D 3206730-320824' #   7 fb^{-1} (#sqrt{s} = 13 TeV)'
    #lat_title = 'Run2018D 1Tier miniAOD'
    lat_title = ptitle[1]+' (13 TeV)'
    #lat_form = '#sigma^{2}_{i}=(#frac{N}{A_{eff}/#sigma_{n}})^{2}+2C^{2}'
    lat_form = '#sigma^{2}_{i}=(N/ETeff)^{2}+2C^{2}'
    #lat.SetTextSize(0.045);
    #lat.SetTextFont(132);
    #lat.DrawLatex(0.48+max(0.,0.32-lat.GetXsize()),0.947,lat_string);
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    #lat.DrawLatex(0.12,0.9325,lat_cms);
    lat.DrawLatex(0.15,0.9325,lat_cms);
    #lat.DrawLatex(0.2,0.9325,lat_cms);
    #lat.SetTextSize(0.04);
    #lat.SetTextFont(42);
    #lat.DrawLatex(0.58,0.9325,lat_title);
    lat.DrawLatex((0.72-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.04);
    lat.DrawLatex(t[0],t[1],ptitle[2]);
    if dofit : 
        lat.SetTextSize(0.03);
        lat.DrawLatex(t[3],t[4]+.05,lat_form);
        #lat.SetTextAlign(12)
        #lat.DrawLatex(.74,0.0325,layout['xtitle']);

    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.Show()
    c1.Close()

from overlay_hist_defs_v3 import *

#dispho_2t_eg2018B_preplot_sp

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

#legtitle = 'Run2016 SP'
#legtitle = 'Run2016 DEG'
#legtitle = 'Run2016G Raw 94 '
#legtitle = 'Run2017 SP'
#legtitle = 'Run2017B Raw 94X'
#legtitle = 'Run2017 DEG'
#legtitle = 'Run2018 EG' 
#legtitle = 'Run2018B EG'
#legtitle = 'Run2018C EG'
#legtitle = 'Run2018D EG'
#legtitle = 'Run2018A Raw '
#legtitle = 'Run2018A Mini '
#legtitle = 'Full dR '
#legtitle = 'dR < 2.5 '
#legtitle = 'dR > 2.5 '
legtitle = ''
#legtitle = 'KuStc'
#legtitle = 'KuNotStc'

rtitle = 'Run2 (13TeV)'

gi_legtitle = 'Global Inclusive'
li_legtitle = 'Local Inclusive'
ls_legtitle = 'Local Same'
Ic_legtitle = ''
loc2_legtitle = ' LS'
#xtitle = 'A_{eff}/#sigma_{n}'
#xtitle = 'ETeff [GeV]'
#xtitle = 'GeV'
#xtitle = 'pCalo time Wt RMS [ns]'
#xtitle = 'multiplcity'
#ytitle = 'RecHit Energy [GeV]'
#xtitle = 'Adjusted pCalo Wt Ave time [ns]'
#xtitle = 'Adjusted pCalo time [ns]'
#ytitle = 'Adjusted pCalo time [ns]'
ytitle = 'CC rechit time [ns]'
#xtitle = 'ADC units'
#xtitle = '#delta(amplitude) [ADC units]'
#xtitle = '# of hits'
#ytitle = '#sigma(t_{1}-t_{2}) [ns]'
#ytitle = 'multiplcity'
#ytitle = 'sum pCalo energy [GeV]'
#ytitle = 'pCalo energy [GeV]'
#xtitle = 'pCalo energy [GeV]'
#xtitle = 'pCalo rechit time [ns]'
xtitle = 'SwissCross'
#ytitle = 'pCalo energy / sum pCalo energy [GeV]'
#xtitle = 'sum pCalo energy [GeV]'
#htitle = 'A_{eff}/#sigma_{n} vs #sigma_{#delta_{t}}'
htitle = ''
#islogx = True
islogx = False
#islogy = True
islogy = False
islogz = True
#islogz = False

gll_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
loc_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle + li_legtitle }
#s2legtitle = '#splitline{'+legtitle+loc2_legtitle+'}{Cali w/ TOF}'
#s2legtitle = '#splitline{'+legtitle+loc2_legtitle+'}{Cali w/out TOF}'
#s2legtitle = '#splitline{'+legtitle+loc2_legtitle+'}{Cali TOF Comp}'
s2_li_legtitle = '#splitline{'+legtitle+'}{'+li_legtitle+'}'
li_legtitle = legtitle+' '+li_legtitle
s2_gi_legtitle = '#splitline{'+legtitle+'}{'+gi_legtitle+'}'
gi_legtitle = legtitle+' '+gi_legtitle
Ic_legtitle = legtitle+' '+Ic_legtitle

loc2_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : li_legtitle }
glo2_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : gi_legtitle }
Ic_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'logz' : islogz, 'legtitle' : Ic_legtitle }
glo_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle + gi_legtitle }


hl_swisscross = [
     ["SwcrVRtTimeOOT","","../../egammares_diag_23D_dvm_Partial_v11_EB.root",""],
]

#ptitle=[' DYJetsToLL','','']
ptitle=['','','#splitline{Energy < 10 GeV}{#splitline{Run 369999}{EB Topo + kOOT Cuts}}']
#ptitle=[' HTo2LongLivedTo4b','','#splitline{10 GeV < rechit Energy < 120 GeV}{With OOT Amplitude Correction}']
#ptitle=[' HTo2LongLivedTo4b','','#splitline{10 GeV < rechit Energy < 120 GeV}{Without OOT Amplitude Correction}']
#ptitle=[' HTo2LongLivedTo4b','','#splitline{10 GeV < rechit Energy < 120 GeV}{Adjusted OOT Amplitude Correction}']
#ptitle=[' DYJetsToLL','','#splitline{#splitline{ECAL Barrel}{#splitline{Calibrated pCalo time}{pCalo energy > 1.5 GeV}}}{4 Crystal}']
#ptitle=[' DYJetsToLL','','#splitline{#splitline{ECAL Barrel}{pCalo time Wt RMS}}{4 Crystal}']
#ptitle=[' DYJetsToLL','','#splitline{#splitline{ECAL Barrel}{pCalo Wt Time}}{iEta 66, iPhi 316}']
#ptitle=[' DYJetsToLL','','#splitline{#splitline{ECAL Barrel}{pCalo Wt Time}}{iEta -4, iPhi 184}']
#ptitle=[' DYJetsToLL','','#splitline{#splitline{ECAL Barrel}{pCalo Wt Time}}{iEta 26, iPhi 357}']
#ptitle=[' DYJetsToLL','','#splitline{#splitline{ECAL Barrel}{pCalo Wt Time}}{iEta -48, iPhi 61}']
#z = [10, 1000000]
z = [1, 10000]
#y = [ 60, 0 ]
#y = [ 280, 0.8 ]
#x = [ -0.01, 0.5 ]
#x = [ -0.5, 20.5 ]
#x = [ -1.5, 1.5 ]
#x = [ 0, 25 ]
#y = [0, 1]
#y = [-1.75, 1.75]
#x = [ -0.05, 0.4 ]
#x = [ -2.25, 19 ]
#y = [-2.25, 14]
#x = [ -1.25, 8.25 ]
#y = [-1.25, 6 ]
x = [ 0.1, 1.1 ]
y = [-15, 15 ]
#x = [ -1.25, 3.5 ]
#y = [-1.25, 2.375 ]
#l = [ 0.7,0.65,0.925,0.9 ]
l = [ 0.7,0.66,0.9,0.9]
#t = [0.7,0.775,0.0]
#t = [0.225,0.85,0.0]
t = [0.2,0.225,0.0]
#outname = 'downloads/tr_pcaloMtime_hist'
#outname = 'downloads/tr_pcaloMulti3_hist'
#outname = 'downloads/tr_pcaloWtRMS3_hist'
#outname = 'downloads/tr_pcalo_adj_tve3_hist'
#outname = 'downloads/tr_2d_pcwtt_v_pce'
#outname = 'downloads/tr_2d_pctvne'
#outname = 'downloads/tr_2d_pcrmsvmulti'
#outname = 'downloads/tr_2d_pcrmsvspce'
#outname = 'downloads/tr_2d_spcevmulti'
#outname = 'downloads/tr_2d_apctvnpce'
#outname = 'downloads/tr_2d_delay_orig_zoom_pcvcc'
outname = 'tr_2d_scvttopo_r3_eg0'

dostack(hl_swisscross, outname, date, Ic_layout, ptitle, z, y, x, l, t)
#dostack(hl_sp_plots, outname, date, Ic_layout, ptitle, z, y, x, l, t)

