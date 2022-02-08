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
    if layout['logy'] : c1.SetLogy()
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
        hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',75,1275,2)
        #hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',75,1275,3)
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',6,100,2)
        hfit.SetParName(0,'N')
        hfit.SetParameter(0,20)
        hfit.SetParLimits(0,10,50)
        #hfit.SetParameter(0,5)
        #hfit.SetParLimits(0,0,10)
        hfit.SetParName(1,'C')
        hfit.SetParameter(1,0.1)
        hfit.SetParLimits(1,0.02,0.3)
        #hfit.SetParameter(1,0.05)
        #hfit.SetParLimits(1,0.001,1.0)
        #hfit.SetParName(2,'S')
        #hfit.SetParameter(2,100.0)
        #hfit.SetParLimits(2,0.00001,100000.0)

    mg = TMultiGraph();

    thebinmid.append([  88.01, 111.9, 136.4, 161.4, 196.6, 247.8, 298.5, 348.9, 419.7, 528.9, 656.9, 801.5,   990, 1487.5,  1975])
    thebinerror.append([ 7.09, 7.156, 7.154, 7.172, 14.28, 14.39, 14.43,  14.4, 28.53, 35.37, 40.58, 43.16, 38.06,      0,     0])
    thebinmid.append([  90.43, 112.6, 137.1, 161.8, 196.8, 247.1, 297.6, 347.9, 417.6, 528.2,   663, 825.1,  1022,   1320,  1975])
    thebinerror.append([6.521, 7.156, 7.192, 7.192, 14.22, 14.29, 14.33, 14.34,  28.4, 35.39, 42.34, 53.72, 63.78,  51.57,     0])
    thebinmid.append([  94.04,   115, 137.7, 162.2, 198.1, 247.7, 297.7, 347.9,   418, 528.4, 664.2, 833.9,  1064,   1387,  1779])
    thebinerror.append([5.056, 6.691, 7.158, 7.201, 14.32, 14.31, 14.33, 14.36, 28.43, 35.41, 42.68, 56.53,  84.3,  81.16, 65.90])

    for histname, tree, infile, timefitfile, tdhistfile, lego in hist_list :
     
        f1.append(TFile.Open(infile+timefitfile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
 
        xbinsfile = TFile.Open(infile+tdhistfile)
        xbinmean = []
        xbinerror = []
        xbhist = []
        for xbn in range(0,17):
          binfilename = 'bin'+str(xbn)+'hist'
          print( 'pulling :' + binfilename )
          xbhist.append(xbinsfile.Get(binfilename))
          xbinmean.append(xbhist[xbn].GetMean(1))
          xbinerror.append(xbhist[xbn].GetStdDev(1))

        orighist = f1[n].Get(hist)
        nbins = 15 #orighist.GetNbinsX()
        low = 0.0 #orighist.GetMinimumBin()
        high = 100.0 #orighist.GetMaximumBin()
        bscale = (high-low)/(10*nbins)
        #print('Find Bins',nbins,low,high,bscale) 
        htitle = 'hist' + str(n)
        #mybins = array([0,2,4,6,8,10,12,14,16,20,24,30,38,52,68,100],dtype='float64')
        #mybins = [0,2,4,6,8,10,12,14,16,20,24,30,38,52,68,100]
        #binwidths = [     1,1,1,1, 1, 1, 1, 1, 2, 2, 3, 4, 7, 8, 16]
        #numExcludeBins = 0
        mybins = [75,  100, 125, 150, 175, 225, 275, 325, 375, 475, 600, 750,  950, 1275, 1700, 2250]
        binwidths = [   12.5,12.5,12.5,12.5,25.0,25.0,25.0,25.0,50.0,62.5,75.0,100.0,162.5,212.5,275.0]
        numExcludeBins = 0 

        #binwidths = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]        

        lenmybins = len(mybins) - numExcludeBins
        #myhist = TH1D(htitle,"",150,low,high)
        #h1.append(TGraphErrors(lenmybins))
        h1.append(TGraphErrors(lenmybins))
        #h1.append(TH2D(htitle,"",(10*nbins),low,high,240,0,2.4))
        #h1.append(TH1F(htitle,"",len(mybins)-1,mybins))

        for bn in range( 1,lenmybins):
            binval = float(orighist.GetBinContent(bn)) #*binwidths[bn-1]*2
            #binval = float(orighist.GetBinContent(bn))
            #binerr = 0.0
            binerr = float(orighist.GetBinError(bn)) #*binwidths[bn-1]*2
            #binmid = float(orighist.GetBinCenter(bn))
            binmid = float(xbinmean[bn])
            #binmid = thebinmid[n][bn-1]  #float(orighist.GetBinCenter(bn))
            if binval > 0.0 :     
                #print('Fill bin',bn,'at',binmid,'with',binval,binerr)
                h1[n].SetPoint(bn,binmid,binval)
                #ebin = int(binmid/bscale)+1
                #widtherr = binwidths[bn-1]/sqrt(3)
                widtherr = float(xbinerror[bn])
                #widtherr = thebinerror[n][bn-1] #binwidths[bn-1]/sqrt(3)
                h1[n].SetPointError(bn,widtherr,binerr)
                print('Fill bin',bn,'at',binmid,'with',binval,'for',widtherr,'by',binerr) 
                #print('Fill bin',bn,'at',binmid,'with',binval,'for',ebin,'with',binerr)  

#        h1.append(f1[n].Get(hist))

#       print( f1 )
#       print( h1 )

        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #h1[n].SetMarkerStyle(6)
        #h1[n].SetTitle(layout['title'])
        #h1[n].GetXaxis().CenterTitle(True)
        #h1[n].GetXaxis().SetTitle(layout['xtitle'])
        #h1[n].GetYaxis().CenterTitle(True)
        #h1[n].GetYaxis().SetTitle(layout['ytitle'])
#       k = [kMagenta+2,kBlue+1,kAzure+4,kBlack,kYellow+1,kViolet+7,kOrange+7,kRed+2,kGreen+3, kGray]
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlack]
#       k = [kMagenta+2,kBlue+2,kGreen+2]
        h1[n].SetLineColor(k[n])
        if dofit : hfit.SetLineColor(k[n])
        #if dofit : hfit.SetLineColor(kRed+2)
        h1[n].SetMarkerColor(k[n])
        #msz = 0.4
        msz = 0.8
        #msz = 1.2
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        #h1[n].SetMinimum(y[1])
        #h1[n].SetMaximum(y[0])
#       h1[n].GetXaxis().SetRangeUser(200.0,1100.0)
        #h1[n].GetXaxis().SetRangeUser(x[0],x[1])
        #if layout['logx'] : h1[n].GetXaxis().SetMoreLogLabels()
        #if layout['logy'] : h1[n].GetYaxis().SetMoreLogLabels()
        #c1.cd()
        #c1.Update()
        if( first ) :
                #h1[n].Draw('ep')
                #h1[n].Draw('AP')
                if dofit : h1[n].Fit(hfit,'REQ')
                first = False
        else :
                #h1[n].Draw('epSAME')
                #h1[n].Draw('SAME AP')
                if dofit : h1[n].Fit(hfit,'REQ+') 
        #c1.Update()
        #c1.Draw()

        if dofit : 
                 paramn.append(str(hfit.GetParameter(0)))
                 paramc.append(str(hfit.GetParameter(1)))
                 #params.append(str(hfit.GetParameter(2)))
                 pne = hfit.GetParError(0)
                 pce = hfit.GetParError(1)
                 #pse = hfit.GetParError(2)
                 print('Fit info',paramn[n],pne,paramc[n],pce)
                 #print('Fit info',paramn[n],pne,paramc[n],pce,params[n],pse)
                 if pne < 0.01 : pne = 0.01
                 if pce < 0.0001 : pce = 0.0001
                 parnerror.append(str(pne))
                 parcerror.append(str(pce))
                 #parserror.append(str(pse))
                 #print( 'C: ' + param + ' +/- ' + parerror )
                 #lat_param = '#color['+str(k[n])+']{N : ' + paramn[0:4] + ' #pm ' + parnerror[0:4] + ' [ns]  C : ' + paramc[0:6] + ' #pm ' + parcerror[0:6] + ' [ns]}'
                 #lat.SetTextSize(0.03);
                 #lat.DrawLatex(t[3],t[4]-n*.035,lat_param);    
                 #c1.Modified()
                 #c1.Update()       
 
        legend.AddEntry(h1[n],lego,'epl');
        #c1.Update()
        #c1.Draw()
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
    mg.SetMinimum(y[1])
    mg.SetMaximum(y[0])
#   mg.GetXaxis().SetRangeUser(200.0,1100.0)
    mg.GetXaxis().SetRangeUser(x[0],x[1])
    if layout['logx'] : mg.GetXaxis().SetMoreLogLabels()
    if layout['logy'] : mg.GetYaxis().SetMoreLogLabels()

    #mg.Draw('AP')
    legend.Draw('same')  #   legend inclusion switch

    gPad.Modified()
    #c1.SetGridx(1)
    #c1.SetGridy(1)
    #c1.cd()
    #c1.Update()
    #mg.Draw('AP')

    lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    #lat_cms = '#bf{CMS} #it{Work in Progress}' + ptitle[0]
    #lat_title = 'Run2018D 3206730-320824' #   7 fb^{-1} (#sqrt{s} = 13 TeV)'
    #lat_title = 'Run2018D 1Tier miniAOD'
    lat_title = ptitle[1]+' (13 TeV)'
    lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + #frac{S^{2}}{A_{eff}/#sigma_{n}} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (N/Eeff)^{2} + 2C^{2}'
    #lat.SetTextSize(0.045);
    #lat.SetTextFont(132);
    #lat.DrawLatex(0.48+max(0.,0.32-lat.GetXsize()),0.947,lat_string);
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    #lat.DrawLatex(0.12,0.9325,lat_cms);
    lat.DrawLatex(0.15,0.9325,lat_cms);
    #lat.DrawLatex(0.25,0.9325,lat_cms);
    #lat.SetTextSize(0.04);
    #lat.SetTextFont(42);
    #lat.DrawLatex(0.58,0.9325,lat_title);
    lat.DrawLatex((0.82-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.04);
    lat.DrawLatex(t[0],t[1],ptitle[2]);
    if dofit : 
        lat.SetTextSize(0.03);
        lat.DrawLatex(t[3],t[4]+.075,lat_form);
        #lat.SetTextAlign(12)
        #lat.DrawLatex(.74,0.0325,layout['xtitle']);
        for l in range(0,n):
            lat_param ='#color['+str(k[l])+']{N : '+paramn[l][0:4]+' #pm '+parnerror[l][0:4]+' [ns]  C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' }'
            #lat_param ='#color['+str(k[l])+']{N : '+paramn[l][0:4]+' #pm '+parnerror[l][0:4]+' [ns]  C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' [ns] S : '+params[l][0:6]+' #pm '+parserror[l][0:6]+' }'
            lat.SetTextSize(0.03);
            lat.DrawLatex(t[3],t[4]-l*.035,lat_param);

    
    if layout['logx'] : c1.SetLogx()
    ####c1.BuildLegend()
    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.SaveAs( outname + '_' + date + '.C' )
    #c1.Show()
    c1.Close()

from overlay_hist_defs_v2 import *

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
xtitle = 'A_{eff}/#sigma_{n}'
#xtitle = 'Eeff [GeV]'
#xtitle = 'GeV'
#xtitle = 'ns'
#ytitle = '#sigma(t_{1}-t_{2}) [ns]'
ytitle = '#occupancy(t_{1}-t_{2}) '
#ytitle = ''
#htitle = 'A_{eff}/#sigma_{n} vs #sigma_{#delta_{t}}'
htitle = ''
islogx = True
#islogx = False
islogy = True
#islogy = False

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
Ic_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : Ic_legtitle }
glo_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle + gi_legtitle }



ptitle=[' 2017','','#splitline{EBEB}{Same RO unit}']
y = [ 10000000, 0.01 ]
x = [ 75.0, 2250.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.82,0.0,0.175,0.28]
outname = 'tr_mini_occ_78_plot'
dostack(hl_mini_87_occ, outname, date, Ic_layout, ptitle, y, x, l, t)

#--------------------------------------------------------------------------------------------
ptitle=[' 2018','','#splitline{#splitline{EB+ iEta 66-85}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 10000000, 1 ]
x = [ 75.0, 2250.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.82,0.0,0.175,0.28]
outname = 'downloads/tr_mini_occ_dzeta_6690_plot'
#dostack(hl_mini8_occ_dzeta6690, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB+ iEta 01-25}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_0025_plot'
#dostack(hl_mini8_dzeta0025, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB+ iEta 26-45}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_2645_plot'
#dostack(hl_mini8_dzeta2645, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB+ iEta 46-65}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_4665_plot'
#dostack(hl_mini8_dzeta4665, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB+ iEta 66-85}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_6690_plot'
#dostack(hl_mini8_dzeta6690, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB- iEta 01-25}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_m2500_plot'
#dostack(hl_mini8_dzetam2500, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB- iEta 26-45}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_m45m26_plot'
#dostack(hl_mini8_dzetam45m26, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB- iEta 46-65}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_m65m46_plot'
#dostack(hl_mini8_dzetam65m46, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB- iEta 66-85}{Runs 318100 - 322800}}{Z -> ee}']
y = [ 0.4, 0.18 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_dzeta_m90m66_plot'
#dostack(hl_mini8_dzetam90m66, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2017','','ECAL Barrel']
y = [ 1.0, 0.1 ]
x = [ 75.0, 1275.0 ]
l = [ 0.75,0.65,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.275]
outname = 'downloads/tr_mini_7se_comp_locd_gt'
#dostack(hl_7se_comp_ld, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB Abs}{Runs 318100 - 322800}}{Same RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_eta_8a4a_locs_plot'
#dostack(hl_mini_eta_8a_4a_ls, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB Minus}{Runs 318100 - 322800}}{Same RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_eta_8a4m_locs_plot'
#dostack(hl_mini_eta_8a_4m_ls, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB Plus}{Runs 318100 - 322800}}{Same RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_eta_8a4p_locs_plot'
#dostack(hl_mini_eta_8a_4_ls, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB Abs}{Runs 318100 - 322800}}{Different RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_eta_8a4a_locd_plot'
#dostack(hl_mini_eta_8a_4a_ld, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB Minus}{Runs 318100 - 322800}}{Different RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_eta_8a4m_locd_plot'
#dostack(hl_mini_eta_8a_4m_ld, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{EB Plus}{Runs 318100 - 322800}}{Different RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.815,0.0,0.175,0.225]
outname = 'downloads/tr_mini_eta_8a4p_locd_plot'
#dostack(hl_mini_eta_8a_4_ld, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{#splitline{ECAL Barrel}{Runs 315000 - 320200}}{Different RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0,0.175,0.225]
outname = 'downloads/tr_mini_eta_8abc_locd_plot'
#dostack(hl_mini_eta_8abc, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{ECAL Barrel}{Z->ee}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_mini_er_8_glo_plot'
#dostack(hl_mini_er_8_glo, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{ECAL Barrel}{Same RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_mini_er_8_loc_plot'
#dostack(hl_mini_er_8_loc, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{ECAL Barrel}{Different RO unit}']
y = [ 0.5, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_mini_er_8_loc_dif_plot'
#dostack(hl_mini_er_8_loc_dif, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018A','','#splitline{ECAL Barrel}{Same RO unit}']
y = [ 0.35, 0.1 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0,0.175,0.28]
outname = 'downloads/tr_mini_er_8a_loc_plot'
#dostack(hl_mini_er_8_loc, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{ECAL Barrel}{Z->ee}']
y = [ 1000000000, 1 ]
x = [ 75.0, 2250.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_mini_oer_8_glo_plot'
#dostack(hl_mini_oer_8_glo, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{ECAL Barrel}{Different RO unit}']
y = [ 1000000000, 1 ]
x = [ 75.0, 2250.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_mini_oer_8_loc_dif_plot'
#dostack(hl_mini_oer_8_loc_dif, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018','','#splitline{ECAL Barrel}{Same RO unit}']
y = [ 1000000000, 1 ]
x = [ 75.0, 2250.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0,0.175,0.28]
outname = 'downloads/tr_mini_oer_8a_loc_plot'
#dostack(hl_mini_oer_8_loc, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018A','','#splitline{ECAL Barrel}{Different RO unit}']
y = [ 0.35, 0.1 ]
x = [ 100.0, 750.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0,0.175,0.28]
outname = 'downloads/tr_mini_er_8a_loc_dif_plot'
#dostack(hl_mini_er_8a_loc_dif, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018D','','#splitline{ECAL Barrel}{#splitline{Z -> ee}{Run 320807-320824}}']
y = [ 0.50, 0.0 ]
x = [ 75.0, 475.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_kurh_s2_glocomp_nology_plot'
#dostack(hl_kurhs_v21_8d_s2_glocomp, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018D','','#splitline{ECAL Barrel}{#splitline{Z -> ee}{Run 320807-320824}}']
y = [ 5.0, 0.01 ]
x = [ 75.0, 475.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_kurh_s2_glocomp_plot'
#dostack(hl_kurhs_v21_8d_s2_glocomp, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018D','','#splitline{ECAL Barrel}{#splitline{Z -> ee}{Run 320807-320824}}']
y = [ 5.0, 0.01 ]
x = [ 75.0, 475.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_kurh_s2_glo_plot'
#dostack(hl_kurhs_v21_8d_s2_glo, outname, date, Ic_layout, ptitle, y, x, l, t)

ptitle=[' 2018D','','#splitline{ECAL Barrel}{#splitline{Same RO unit}{Run 320807-320824}}']
y = [ 2.0, 0.00001 ]
x = [ 75.0, 950.0 ]
l = [ 0.7,0.66,0.9,0.9]
t = [0.325,0.825,0.0]
outname = 'downloads/tr_kurh_s2_plot'
#dostack(hl_kurhs_v21_8d_s2, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2018D','','#splitline{ECAL Barrel}{#splitline{Same RO unit}{Run 320807-320824}}']
#y = [ 0.4, 0.08 ]
#x = [ 75.0, 950.0 ]
#l = [ 0.7,0.72,0.9,0.9 ]
#t = [0.325,0.825,0.0]
#outname = 'downloads/tr_kurh_s1_plot'
#dostack(hl_kurhs_v21_8d_s1, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2018D','','#splitline{ECAL Barrel}{#splitline{Same RO unit}{Run 320807-320824}}']
#y = [ 0.4, 0.12 ]
#x = [ 75.0, 950.0 ]
#l = [ 0.7,0.54,0.9,0.9 ]
#t = [0.325,0.825,0.0]
#outname = 'downloads/tr_calicomp_mini_Run2018D_ica'
#dostack(hl_comp_mini_ica_80_789, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2018C','','#splitline{ECAL Barrel}{#splitline{Z -> ee}{Run 319337-320065}}']
#y = [ 0.5, 0.18  ]
#x = [ 75.0, 600.0 ]
#l = [ 0.7,0.62,0.9,0.9 ]
#t = [0.325,0.825,0.0]
#outname = 'downloads/tr_8c_mini_ele_match_plot'
#dostack(hl_glo_ele, outname, date, Ic_layout, ptitle, y, x, l, t)



#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}'] #'ECAL Barrel Local Resolution'
#y = [ 0.65, 0.06 ]
#x = [ 75.0, 1275.0 ]
#l = [ 0.785,0.55,0.925,0.9 ]
#t = [0.35,0.825,0.09,0.175,0.32]
#outname = 'downloads/tr_mini_6all_loc_gt'
#dostack(hl_6all_loc, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}'] #'ECAL Barrel Local Resolution'
#y = [ 0.65, 0.06 ]
#x = [ 75.0, 1275.0 ]
#l = [ 0.785,0.55,0.925,0.9 ]
#t = [0.35,0.825,0.09,0.175,0.28]
#outname = 'downloads/tr_mini_7all_loc_gt'
#dostack(hl_7all_loc, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}'] #'ECAL Barrel Local Resolution'
#y = [ 0.65, 0.06 ]
#x = [ 75.0, 1275.0 ]
#l = [ 0.785,0.55,0.925,0.9 ]
#t = [0.35,0.825,0.09,0.175,0.28]
#outname = 'downloads/tr_mini_8all_loc_gt'
#dostack(hl_8all_loc, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Z -> ee}']
y = [ 1.0, 0.1 ]
x = [ 6.0, 100.0 ]
l = [ 0.775,0.55,0.925,0.9 ]
t = [0.325,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_2_Eeff_gt'
#dostack(hl_2_e_otg, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Z -> ee}']
y = [ 1.0, 0.1 ]
x = [ 75.0, 1275.0 ]
l = [ 0.775,0.55,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_2_glo_gt'
#dostack(hl_2_otg, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Different RO unit}']
y = [ 1.0, 0.1 ]
x = [ 75.0, 1275.0 ]
l = [ 0.775,0.55,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_2_locd_gt'
#dostack(hl_2_otld, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}']
y = [ 1.0, 0.1 ]
x = [ 75.0, 1275.0 ]
l = [ 0.775,0.55,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_2_loc_gt'
#dostack(hl_2_otl, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Different RO unit}']
y = [ 1.0, 0.1 ]
x = [ 75.0, 1275.0 ]
l = [ 0.775,0.65,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_678_locd_gt'
#dostack(hl_678_otld, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}']
y = [ 1.0, 0.1 ]
x = [ 6.0, 100.0 ]
l = [ 0.775,0.55,0.925,0.9 ]
t = [0.325,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_678_Eeff_otl'
#dostack(hl_678_e_otl, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Z -> ee}']
y = [ 1.0, 0.1 ]
x = [ 6.0, 100.0 ]
l = [ 0.775,0.65,0.925,0.9 ]
t = [0.325,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_678_Eeff_otg'
#dostack(hl_678_e_otg, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}']
y = [ 1.0, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.775,0.65,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_678_loc_gt'
#dostack(hl_678_otl, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}']
y = [ 1.0, 0.1 ]
x = [ 75.0, 950.0 ]
l = [ 0.775,0.65,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_7ul_comp_loc_gt'
#dostack(hl_7ul_comp, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Z -> ee}']
y = [ 1.0, 0.1 ]
x = [ 75.0, 1275.0 ]
l = [ 0.775,0.65,0.925,0.9 ]
t = [0.35,0.825,0.09,0.175,0.225]
outname = 'downloads/tr_mini_678_glo_gt'
#dostack(hl_678_otg, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=[' 2018D','','#splitline{ECAL Barrel}{#splitline{Local}{Run 320807-320824}}']
#y = [ 0.4, 0.009 ]
#x = [ 75.0, 950.0 ]
#l = [ 0.7,0.60,0.85,0.85 ]
#t = [0.325,0.77,0.0]
#outname = 'downloads/tr_kurh_v19_18d_plot'
#dostack(kurh_8d_v19, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2018D','','#splitline{#splitline{ECAL Barrel}{E5 Calibration}}{#splitline{Local}{Run 320807-320824}}']
#y = [ 0.4, 0.009 ]
#x = [ 75.0, 950.0 ]
#l = [ 0.7,0.60,0.85,0.85 ]
#t = [0.325,0.77,0.0]
#outname = 'downloads/tr_8d_v20v19_comp_plot'
#dostack(kurh_8d_vc, outname, date, Ic_layout, ptitle, y, x, l, t)

#eletd6b='downloads/dispho_ot_mini_Run2016B_ele.root'
#eletd6c='downloads/dispho_ot_mini_Run2016C_ele.root'
#eletd6d='downloads/dispho_ot_mini_Run2016D_ele.root'
#eletd6e='downloads/dispho_ot_mini_Run2016E_ele.root'
#eletd6f='downloads/dispho_ot_mini_Run2016F_ele.root'
#eletd6g='downloads/dispho_ot_mini_Run2016G_ele.root'
#eletd6h='downloads/dispho_ot_mini_Run2016H_ele.root'
#
#eletd7b='downloads/dispho_ot_mini_Run2017B_ele.root'
#eletd7c='downloads/dispho_ot_mini_Run2017C_ele.root'
#eletd7d='downloads/dispho_ot_mini_Run2017D_ele.root'
#eletd7e='downloads/dispho_ot_mini_Run2017E_ele.root'
#eletd7f='downloads/dispho_ot_mini_Run2017F_ele.root'
#
#eletd8a='downloads/dispho_ot_mini_Run2018A_ele.root'
#eletd8b='downloads/dispho_ot_mini_Run2018B_ele.root'
#eletd8c='downloads/dispho_ot_mini_Run2018C_ele.root'
#eletd8d='downloads/dispho_ot_mini_Run2018D_ele.root'
#
##eletd=[ eletd6b, eletd6c, eletd6d, eletd6e, eletd6f, eletd6g, eletd6h, eletd7b, eletd7c, eletd7d, eletd7e, eletd7f, eletd8a, eletd8b, eletd8c, eletd8d ]
#eletd=[eletd6b]
#
#for elefile in eletd :
#
#    hl_glo_eletd = [
#         ["noEleMatchHist","",elefile, "default"],
#         ["eleMatchTrueHist","",elefile, "e match true"],
#         ["eleMatchFalseHist","",elefile, "e match false"],
#    ]
#
#    ptitle=[' '+elefile[25:-9],'','#splitline{ECAL Barrel}{Global}']
#    y = [ 10000000.0, 0.9 ]
#    x = [ -10.0, 10.0 ]
#    l = [ 0.7,0.60,0.85,0.85 ]
#    t = [0.20,0.77,0.0]
#    outname = 'downloads/tr_'+elefile[25:-5]+'_match_tds_plot'
#    dostack(hl_glo_eletd, outname, date, Ic_layout, ptitle, y, x, l, t)


#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Global}']
#y = [ 0.58, 0.12 ]
#x = [ 0.0, 950.0 ]
#l = [ 0.675,0.55,0.925,0.9 ]
#t = [0.35,0.79,0.08,0.225,0.27]
#outname = 'downloads/tr_678_mini_ele_match_plot'
#dostack(hl_678_glo_ele, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=[' 2018C','','#splitline{ECAL Barrel}{#splitline{Global}{Run 319337-320065}}']
#y = [ 0.9, 1000000.0 ]
#x = [ -10.0, 10.0 ]
#l = [ 0.7,0.60,0.85,0.85 ]
#t = [0.20,0.77,0.0]
#outname = 'downloads/tr_8c_mini_ele_match_td_plot'
#dostack(hl_glo_eletd, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2018C','','#splitline{ECAL Barrel}{#splitline{Global}{Run 319337-320065}}']
#y = [ 0.001, 10000.0 ]
#x = [ -5.0, 5.0 ]
#l = [ 0.7,0.60,0.85,0.85 ]
#t = [0.20,0.77,0.0]
#outname = 'downloads/tr_8c_mini_ele_match_bins_plot'
#dostack(hl_glo_ele_bins, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2017D','','#splitline{ECAL Barrel}{Comp}']      ##splitline{Global}{Run 319337-320065}}']
#y = [ 0.45, 0.1  ]
#x = [ 75.0, 950.0 ]
#l = [ 0.7,0.60,0.85,0.85 ]
#t = [0.325,0.77,0.0]
#outname = 'downloads/tr_7d_mini_ele_match_plot'
#dostack(hl_7d_ele, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2018C','','#splitline{ECAL Barrel}{#splitline{Global}{Run 319337-320065}}']
#y = [ 0.35, 0.17  ]
#x = [ 75.0, 600.0 ]
#l = [ 0.7,0.60,0.85,0.85 ]
#t = [0.325,0.77,0.0]
#outname = 'downloads/tr_8c_mini_ele_match_plot'
#dostack(hl_glo_ele, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{#splitline{Global}{Z->ee mass}}']
#y = [ 1000000, 500  ]
#x = [ 0.0, 160.0 ]
#l = [ 0.8,0.60,0.90,0.85 ]
#t = [0.225,0.78,0.1]
#outname = 'downloads/tr_678_mini_global_plots'
#dostack(hl_gp678, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle='Run2018D miniAOD '
#y = [ 0.22, 0.1 ]
#x = [ 100.0, 1500.0 ]
#l = [ 0.745,0.40,0.925,0.9 ]
#t = 0.68
#outname = 'downloads/tr_calicomp_mini_Run2018D_ica'
#dostack(hl_comp_mini_ica_80_789, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=['','','ECAL Barrel Run2016']
#y = [ 0.56, 0.1 ]
#x = [ 1.0, 1275.0 ]
#l = [ 0.58,0.50,0.75,0.9 ]
#t = [0.225,0.15,0.0]
#outname = 'downloads/tr_6all_mini_glsd'
#####outname = 'downloads/tr_6all_log_mini_glsd'
#dostack(hl_6_glsd, outname, date, Ic_layout, ptitle, y, x, l, t)
### <<<<<<<<<<<<<<<<<<<<<<<<<<
#ptitle=[' 2017','','ECAL Barrel']
#y = [ 0.56, 0.0 ]
#######y = [ 0.5, 0.09 ]
#x = [ 0, 1275.0 ]
#l = [ 0.64,0.60,0.925,0.9 ]
#t = [0.275,0.83,0,0.225,0.225]
#outname = 'downloads/tr_7all_mini_glsd'
#####outname = 'downloads/tr_7all_mini_log_glsd'
#dostack(hl_7_glsd, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2017','','ECAL Barrel']
#y = [ 0.56, 0.0 ]
#######y = [ 0.5, 0.09 ]
#x = [ 0, 1275.0 ]
#l = [ 0.64,0.60,0.925,0.9 ]
#t = [0.275,0.83,0,0.225,0.225]
#outname = 'downloads/tr_7all2_mini_glsd'
#####outname = 'downloads/tr_7all_mini_log_glsd'
#dostack(hl_7_glsd2, outname, date, Ic_layout, ptitle, y, x, l, t)

####
#ptitle=['','','ECAL Barrel Run2018']
#####y = [ 0.42, 0.12 ]
#y = [ 0.56, 0.1 ]
#x = [ 1, 1275.0 ]
#l = [ 0.58,0.50,0.75,0.9 ]
#t = [0.225,0.15,0]
#outname = 'downloads/tr_8all_mini_glsd'
#####outname = 'downloads/tr_8all_log_mini_glsd'
#dostack(hl_8_glsd, outname, date, Ic_layout, ptitle, y, x, l, t)
#
#ptitle=[' 2018','','#splitline{ECAL Barrel}{#splitline{Same RO unit}{Run 319337-320065}}']   #'Run2018C 319337-320065' KuReco Local '
#y = [ 0.33, 0.0 ]
#x = [ 0.0, 325.0 ]
#l = [ 0.655,0.60,0.925,0.9 ]
#t = [0.325,0.775,0,0.225,0.245]
#outname = 'downloads/tr_kurh_comp_glo'
#dostack(hl_kurh_8c_comp_glo, outname, date, Ic_layout, ptitle, y, x, l, t)
#
#ptitle=[' 2018','','#splitline{ECAL Barrel}{#splitline{Same RO unit}{Run 319337-320065}}']   #'Run2018C 319337-320065' KuReco Local '
#y = [ 0.33, 0.0 ]
#x = [ 0.0, 950.0 ]
#l = [ 0.655,0.60,0.925,0.9 ]
#t = [0.325,0.775,0,0.225,0.245]
#outname = 'downloads/tr_kurh_comp_loc'
#dostack(hl_kurh_8c_comp_loc, outname, date, Ic_layout, ptitle, y, x, l, t)
#
#ptitle=[' 2018','','ECAL Barrel']
#y = [ 0.33, 0.1 ]
#x = [ 0.0, 950.0 ]
#l = [ 0.755,0.60,0.925,0.9 ]
#t = [0.325,0.775,0,0.225,0.257]
#outname = 'downloads/tr_kurh_calicomp_loc'
#dostack(hl_kurh_8c_comp_cali, outname, date, Ic_layout, ptitle, y, x, l, t)
#dostack(hl_kurh_8c_comp_cali_glo, outname, date, Ic_layout, ptitle, y, x, l, t)

#ptitle=[' 2018D','','ECAL Barrel']
#y = [ 0.4, 0.1 ]
#x = [ 0.0, 1200.0 ]
#l = [ 0.745,0.60,0.925,0.9 ]
#t = [0.225,0.85,0]
#outname = 'downloads/tr_ot_8d_calicomp_loc'
#dostack(hl_mini_ica_809l, outname, date, Ic_layout, ptitle, y, x, l, t)
#outname = 'downloads/tr_ot_8d_calicomp_glo'
#dostack(hl_mini_ica_809g, outname, date, Ic_layout, ptitle, y, x, l, t)
#
#ptitle='ECAL Barrel Global Resolution'
#y = [ 0.4, 0.09 ]
#x = [ 10.0, 2000.0 ]
#l = [ 0.745,0.3,0.925,0.9 ]
#t = 0.68
#outname = 'downloads/tr_mini_678all_loc_gt'
#dostack(hl_678all_loc, outname, date, Ic_layout, ptitle,  y, x, l, t)
#
#ptitle='ECAL Barrel Global Resolution'
#y = [ 0.5, 0.09 ]
#x = [ 10.0, 2000.0 ]
#l = [ 0.745,0.3,0.925,0.9 ]
#t = [0.225,0.15]
#outname = 'downloads/tr_mini_678all_glo_gt'
#dostack(hl_678all_glo, outname, date, Ic_layout, ptitle,  y, x, l, t)
#
#ptitle='ECAL Barrel Local Resolution'
#y = [ 0.45, 0.09 ]
#x = [ 10.0, 1275.0 ]
#l = [ 0.695,0.45,0.835,0.9 ]
#t = [0.225,0.15]
#outname = 'downloads/tr_mini_6all_loc_gt'
#dostack(hl_6all_loc, outname, date, Ic_layout, ptitle,  y, x, l, t)
##
#ptitle='ECAL Barrel Global Resolution'
#y = [ 0.45, 0.09 ]
#x = [ 10.0, 1275.0 ]
#l = [ 0.185,0.15,0.325,0.55 ]
#t = [0.525,0.15]
#outname = 'downloads/tr_mini_6all_glo_gt'
#dostack(hl_6all_glo, outname, date, Ic_layout, ptitle,  y, x, l, t)
##
#ptitle='ECAL Barrel Local Resolution'
#y = [ 0.45, 0.09 ]
#x = [ 10.0, 1275.0 ]
#l = [ 0.685,0.45,0.825,0.9 ]
#t = [0.225,0.15]
#outname = 'downloads/tr_mini_7all_loc_gt'
#dostack(hl_7all_loc, outname, date, Ic_layout, ptitle,  y, x, l, t)
##
#ptitle='ECAL Barrel Global Resolution'
#y = [ 0.45, 0.09 ]
#x = [ 10.0, 1275.0 ]
#l = [ 0.785,0.45,0.925,0.9 ]
#t = [0.225,0.15]
#outname = 'downloads/tr_mini_7all_glo_gt'
#dostack(hl_7all_glo, outname, date, Ic_layout, ptitle,  y, x, l, t)
##
#ptitle='ECAL Barrel Local Resolution'
#y = [ 0.45, 0.09 ]
#x = [ 10.0, 1275.0 ]
#l = [ 0.785,0.45,0.925,0.9 ]
#t = [0.225,0.15]
#outname = 'downloads/tr_mini_8all_loc_gt'
#dostack(hl_8all_loc, outname, date, Ic_layout, ptitle,  y, x, l, t)
##
#ptitle='ECAL Barrel Global Resolution'
#y = [ 0.45, 0.09 ]
#x = [ 10.0, 1275.0 ]
#l = [ 0.785,0.45,0.925,0.9 ]
#t = [0.225,0.15]
#outname = 'downloads/tr_mini_8all_glo_gt'
#dostack(hl_8all_glo, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Global}']
#y = [ 0.62, 0.0 ]
#x = [ 0.0, 950.0 ]
#l = [ 0.775,0.55,0.925,0.9 ]
#t = [0.35,0.79,0.08,0.225,0.27]
#outname = 'downloads/tr_mini_678_glo_gt'
#dostack(hl_678_glo, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Global}']
#y = [ 0.62, 0.0 ]
#x = [ 0.0, 950.0 ]
#l = [ 0.775,0.55,0.925,0.9 ]
#t = [0.35,0.79,0.08,0.225,0.27]
#outname = 'downloads/tr_mini_678_glo2_gt'
#dostack(hl_678_glo2, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}']
#y = [ 0.4, 0.1 ]
#x = [ 75.0, 1700.0 ]
#l = [ 0.775,0.55,0.925,0.9 ]
#t = [0.35,0.79,0.08,0.225,0.62]
#outname = 'downloads/tr_mini_678_lg_gt'
#dostack(hl_678_lg, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=['','Run 2 ','#splitline{ECAL Barrel}{Same RO unit}']
#y = [ 0.38, 0.0 ]
#x = [ 0.0, 950.0 ]
#l = [ 0.775,0.55,0.925,0.9 ]
#t = [0.35,0.79,0.08,0.225,0.27]
#outname = 'downloads/tr_mini_678_lg2_gt'
#dostack(hl_678_lg2, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle=[' 2017','','#splitline{ECAL Barrel}{Same RO unit}']
#y = [ 0.38, 0.0 ]
#x = [ 0.0, 950.0 ]
#l = [ 0.725,0.55,0.925,0.9 ]
#t = [0.35,0.79,0,0.225,0.275]
#outname = 'downloads/tr_mini_678_compcali_gt'
#dostack(hl_678_compcali, outname, date, Ic_layout, ptitle,  y, x, l, t)

#ptitle='Run2016C 4.3%'
#y = [ 0.42, 0.09 ]
#y = [ 0.42, 0.24 ] # global
#y = [ 12.0, 0.09 ]
#outname = 'downloads/tr_kurh_loc_Run2016c_ica'
#outname = 'downloads/tr_kurh_loc_zm_Run2016c_ica'
#outname = 'downloads/tr_kurh_glo_Run2016c_ica'
#outname = 'downloads/tr_kurh_glo_zm_Run2016c_ica'
#dostack(hl_kurh_6c_loc, outname, date, Ic_layout, ptitle,  y)
#dostack(hl_kurh_6c_glo, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018C 319337-320065'
#y = [ 0.38, 0.12 ]
#y = [ 1.4, 0.12 ]
#outname = 'downloads/tr_kurh_loc_Run2018c_ica'
#outname = 'downloads/tr_kurh_loc_zm_Run2018c_ica'
#outname = 'downloads/tr_kurh_glo_Run2018c_ica'
#outname = 'downloads/tr_kurh_glo_zm_Run2018c_ica'
#dostack(hl_kurh_8c_loc, outname, date, Ic_layout, ptitle,  y)
#dostack(hl_kurh_8c_glo, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018A miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_gt102_Run2018A_ica'
#dostack(hl_gt102_mini_ica_802, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018C 319337-320065 KuStc'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_isoot_Run2018C_ica'
#dostack(hl_isoot_mini_ica_806, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2017F miniAOD'
#y = [ 0.8, 0.25 ]
#outname = 'downloads/tr_eb_epem_Run2017F_ica'
#dostack(hl_ebpm_7f, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018A miniAOD'
#y = [ 0.9, 0.05 ]
#outname = 'downloads/tr_eb_epem_Run2018A_ica'
#dostack(hl_ebpm_8a, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2018B miniAOD'
#y = [ 0.9, 0.05 ]
#outname = 'downloads/tr_eb_epem_Run2018B_ica'
#dostack(hl_ebpm_8b, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2018C miniAOD'
#y = [ 0.9, 0.05 ]
#outname = 'downloads/tr_eb_epem_Run2018C_ica'
#dostack(hl_ebpm_8c, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2018D miniAOD'
#y = [ 0.9, 0.05 ]
#outname = 'downloads/tr_eb_epem_Run2018D_ica'
#dostack(hl_ebpm_8d, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2017D miniAOD'
#y = [ 0.625, 0.25 ]
#outname = 'downloads/tr_eb_epem_Run2017D_ica'
#dostack(hl_ebpm_7d, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2017B miniAOD'
#y = [ 0.625, 0.25 ]
#outname = 'downloads/tr_eb_epem_Run2017B_ica'
#dostack(hl_ebpm_7b, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2016B miniAOD'
#y = [ 0.625, 0.25 ]
#outname = 'downloads/tr_eb_epem_Run2016B_ica'
#dostack(hl_ebpm_6b, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2016E miniAOD'
#y = [ 0.625, 0.25 ]
#outname = 'downloads/tr_eb_epem_Run2016E_ica'
#dostack(hl_ebpm_6e, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2017F miniAOD'
#y = [ 0.625, 0.25 ]
#outname = 'downloads/tr_eb_plus_minus_mini_Run2017F_ica'
#dostack(hl_ebpm_7f, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018D miniAOD'
#y = [ 0.4, 0.15 ]
#outname = 'downloads/tr_eb_plus_minus_mini_Run2018D_ica'
#dostack(hl_ebpm_8d, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018D miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_tofswap_mini_Run2018D_ica'
#dostack(hl_swapTOF_8d, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2017F miniAOD'
#y = [ 0.625, 0.25 ]
#outname = 'downloads/tr_tofswap_mini_Run2017F_ica'
#dostack(hl_swapTOF_7f, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2016G miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2016G_ica'
#dostack(hl_mini_ica_605l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2016F miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2016F_ica'
#dostack(hl_mini_ica_604l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2016E miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2016E_ica'
#dostack(hl_mini_ica_603l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2016D miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2016D_ica'
#dostack(hl_mini_ica_602l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2016C miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2016C_ica'
#dostack(hl_mini_ica_601l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2016B miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2016B_ica'
#dostack(hl_mini_ica_600l, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2017F miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2017F_ica'
#dostack(hl_mini_ica_705l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2017E miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2017E_ica'
#dostack(hl_mini_ica_704l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2017D miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2017D_ica'
#dostack(hl_mini_ica_703l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2017C miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2017C_ica'
#dostack(hl_mini_ica_702l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2017B miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2017B_ica'
#dostack(hl_mini_ica_701l, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018D miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2018D_ica'
#dostack(hl_mini_ica_809l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2018C miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2018C_ica'
#dostack(hl_mini_ica_806l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2018B miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2018B_ica'
#dostack(hl_mini_ica_803l, outname, date, Ic_layout, ptitle,  y)
#
#ptitle='Run2018A miniAOD'
#y = [ 0.35, 0.1 ]
#outname = 'downloads/tr_calicomp_mini_Run2018A_ica'
#dostack(hl_mini_ica_802l, outname, date, Ic_layout, ptitle,  y)

#ptitle='Time Photon Delta 0,1'
#y = [ 450000, 0 ]
#outname = 'downloads/tr_glop_timed'
#dostack(hl_glop_timed, outname, date, Ic_layout, ptitle,  y)

#ptitle='Time Photon 1'
#y = [ 425000, 0 ]
#outname = 'downloads/tr_glop_time1'
#dostack(hl_glop_time1, outname, date, Ic_layout, ptitle,  y)

#ptitle='Time Photon 0'
#y = [ 425000, 0 ]
#outname = 'downloads/tr_glop_time0'
#dostack(hl_glop_time0, outname, date, Ic_layout, ptitle,  y)

#ptitle='TOF Photon Delta 0,1'
#y = [ 1000000.0, 0.01 ]
#outname = 'downloads/tr_glop_tofd'
#dostack(hl_glop_tofd, outname, date, Ic_layout, ptitle,  y)

#ptitle='TOF Photon 1'
#y = [ 1000000.0, 0.01 ]
#outname = 'downloads/tr_glop_tof1'
#dostack(hl_glop_tof1, outname, date, Ic_layout, ptitle,  y)

#ptitle='TOF Photon 0'
#y = [ 1000000.0, 0.01 ]
#outname = 'downloads/tr_glop_tof0'
#dostack(hl_glop_tof0, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018C 1Tier miniAOD'
#y = [ 0.4, 0.1 ]
#outname = 'downloads/tr_calicomp_pho_mini_Run2018C_ica'
#dostack(hl_comp_mini_ica_8020a, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2016C SP 1Tier miniAOD'
#y = [ 0.4, 0.1 ]
#outname = 'downloads/tr_comp_mini_Run2016C_SP_cali_DE_ica'
#dostack(hl_comp_mini_ica_608, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018C 1Tier miniAOD'
#y = [ 0.4, 0.1 ]
#outname = 'downloads/tr_comp_mini_Run2018C_caliD_ica'
#dostack(hl_comp_mini_ica_8019, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018A 1Tier miniAOD'
#y = [ 0.4, 0.1 ]
#outname = 'downloads/tr_comp_mini_Run2018A_caliall_ica'
#dostack(hl_comp_mini_ica_8017, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018B 1Tier miniAOD'
#y = [ 0.4, 0.1 ]
#outname = 'downloads/tr_comp_mini_Run2018B_calliall_ica'
#dostack(hl_comp_mini_ica_8018, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018C 1Tier miniAOD'
#y = [ 0.4, 0.1 ]
#outname = 'downloads/tr_comp_mini_Run2018C_calliall_ica'
#dostack(hl_comp_mini_ica_8020, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018D 1Tier miniAOD'
#y = [ 0.4, 0.1 ]
#outname = 'downloads/tr_comp_mini_Run2018D_caliall_ica'
#dostack(hl_comp_mini_ica_8023, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2018D miniAOD'
#y = [ 0.19, 0.125 ]
#outname = 'downloads/tr_calicomp_mini_Run2018D_Stats_ica'
#dostack(hl_comp_mini_ica_809_stats, outname, date, Ic_layout, ptitle,  y)

#ptitle='Run2016C SP 1Tier miniAOD'
#outname = 'downloads/tr_comp_mini_Run2016C_SP_ica'
#dostack(hl_comp_mini_ica_606, outname, date, Ic_layout, ptitle)

#outname = 'downloads/tr_comp_mini_Run2017F_ica'
#dostack(hl_comp_mini_ica_705, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2017E_ica'
#dostack(hl_comp_mini_ica_704, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2017D_ica'
#dostack(hl_comp_mini_ica_703, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2017C_ica'
#dostack(hl_comp_mini_ica_702, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2017B_ica'
#dostack(hl_comp_mini_ica_701, outname, date, Ic_layout)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#ptitle='Run2018D miniAOD EBEB'
#y = [ 0.22, 0.1 ]
#x = [ 100.0, 1500.0 ]
#l = [ 0.745,0.40,0.925,0.9 ]
#t = 0.58
#outname = 'downloads/tr_calicomp_mini_Run2018D_ica'
#dostack(hl_comp_mini_ica_80_789, outname, date, Ic_layout, ptitle, y, x, l, t)

#outname = 'downloads/tr_calicomp_mini_Run2018C_ica'
#dostack(hl_comp_mini_ica_80_456, outname, date, Ic_layout)

#outname = 'downloads/tr_calicomp_mini_Run2018A_ica'
#dostack(hl_comp_mini_ica_80_012, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2016G_ica'
#dostack(hl_comp_mini_ica_605, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2016F_ica'
#dostack(hl_comp_mini_ica_604, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2016E_ica'
#dostack(hl_comp_mini_ica_603, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2016D_ica'
#dostack(hl_comp_mini_ica_602, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2016C_ica'
#dostack(hl_comp_mini_ica_601, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2016B_ica'
#dostack(hl_comp_mini_ica_600, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2018D_ica'
#dostack(hl_comp_mini_ica_809, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2018C_ica'
#dostack(hl_comp_mini_ica_806, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2018B_ica'
#dostack(hl_comp_mini_ica_803, outname, date, Ic_layout)

#outname = 'downloads/tr_comp_mini_Run2018A_ica'
#dostack(hl_comp_mini_ica_802, outname, date, Ic_layout)

#outname = 'downloads/tr_icg_comp_mini'
#dostack(hl_loc_mini_icg, outname, date, Ic_layout)

#outname = 'downloads/tr_glo_ica_comp_mini'
#dostack(hl_glo_mini_ica, outname, date, glo2_layout)

#outname = 'downloads/tr_loc_ica_comp_mini'
#dostack(hl_loc_mini_ica, outname, date, loc2_layout)

#outname = 'downloads/tr_glo_ica_comp_kustc'
#dostack(hl_glo_kustc_ica, outname, date, glo2_layout)

#outname = 'downloads/tr_glo_ica_comp_kunotstc'
#dostack(hl_glo_kunotstc_ica, outname, date, glo2_layout)

#outname = 'downloads/tr_loc_ica_comp_kunotstc'
#dostack(hl_loc_kunotstc_ica, outname, date, loc2_layout)

#outname = 'downloads/tr_loc_ica_comp_kustc'
#dostack(hl_loc_kustc_ica, outname, date, loc2_layout)

#outname = 'downloads/tr_glo_icg_comp_kunotstc'
#dostack(hl_glo_kunotstc_icg, outname, date, glo2_layout)

#outname = 'downloads/tr_loc_icg_comp_kunotstc'
#dostack(hl_loc_kunotstc_icg, outname, date, loc2_layout)

#outname = 'downloads/tr_icgv_comp'
#dostack(hist_list_preplot_kustc_icgv, outname, date, loc2_layout)

#outname = 'downloads/dispho_tt_eg2018B_local_ic_preplot/tr_dispho_preplot_local_ic_tt_eg2018B'
#dostack(hist_list_preplot_local_comp_tt_2018B_ic, outname, date, loc2_layout)

#outname = 'downloads/dispho_tt_eg2018B_global_ic_preplot/tr_dispho_preplot_global_ic_tt_eg2018B'
#dostack(hist_list_preplot_global_comp_tt_2018B_ic, outname, date, glob2_layout)

#outname = 'downloads/dispho_2t_eg2018B_global_preplot/tr_dispho_preplot_global_comp_2t_eg2018B'
#dostack(hist_list_preplot_global_comp_2t_2018B, outname, date, glob2_layout)

#outname = 'downloads/dispho_tt_2018B_preplot_localta/tr_dispho_preplot_local_comp_2t_eg2018B'
#dostack(hist_list_preplot_local_comp_tt_2018Bp, outname, date, loc2_layout)

#outname = 'downloads/dispho_tt_2018C_preplot_localta/tr_dispho_preplot_local_comp_2t_eg2018C'
#dostack(hist_list_preplot_local_comp_tt_2018Cp, outname, date, loc2_layout)

#outname = 'downloads/dispho_tt_2018D_preplot_localta/tr_dispho_preplot_local_comp_2t_eg2018D'
#dostack(hist_list_preplot_local_comp_tt_2018Dp, outname, date, loc2_layout)

#outname = 'downloads/dispho_tt_2018B_preplot_globalta/tr_dispho_preplot_global_comp_2t_eg2018B'
#dostack(hist_list_preplot_global_comp_tt_2018Bp, outname, date, glob2_layout)

#outname = 'downloads/dispho_tt_2018C_preplot_globalta/tr_dispho_preplot_global_comp_2t_eg2018C'
#dostack(hist_list_preplot_global_comp_tt_2018Cp, outname, date, glob2_layout)

#outname = 'downloads/dispho_2t_eg2018D_global_icv2_i25_e5e3_preplot/tr_dispho_preplot_global_comp_2t_eg2018D'
#dostack(hist_list_preplot_global_comp_tt_2018Dp, outname, date, glob2_layout)

#outname = 'downloads/dispho_tt_2018B_preplot_localta/tr_dispho_preplot_local_comp_2t_eg2018B'
#dostack(hist_list_preplot_local_comp_tt_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_tt_2018B_preplot_localta/tr_dispho_preplot_local_comp_2t_eg2018B'
#dostack(hist_list_preplot_local_comp_tt_2018B_proposal, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_global_preplot/tr_dispho_preplot_global_all_2t_eg2018B'
#dostack(hist_list_preplot_global_all_2t_2018B, outname, date, glob2_layout)

#outname = 'downloads/dispho_2t_eg2018B_local_icv2_i25_e5e5_preplot/tr_dispho_preplot_local_all_2t_eg2018B'
#dostack(hist_list_preplot_local_all_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_preplot_sp/tr_dispho_preplot_local_all_2t_eg2018B'
#dostack(hist_list_preplot_local_all_2t_2018B_sp, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_global_preplot/tr_dispho_preplot_global_stats_2t_eg2018B'
#dostack(hist_list_preplot_global_stats_2t_2018B, outname, date, glob2_layout)

#outname = 'downloads/dispho_2t_eg2018B_local_preplot/tr_dispho_preplot_local_stats_2t_eg2018B'
#dostack(hist_list_preplot_local_stats_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_notof_kuStc_preplot/tr_dispho_preplot_kustc_2t_eg2018B'
#dostack(hist_list_preplot_kustc_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_notof_preplot/tr_dispho_preplot_comp_2t_eg2018B'
#dostack(hist_list_preplot_comp2_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_tof_preplot/tr_dispho_preplot_tofcomp_2t_eg2018B'
#dostack(hist_list_preplot_tofcomp_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_notof_preplot/tr_dispho_preplot_notof_2t_eg2018B'
#dostack(hist_list_preplot_notof_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_tof_preplot/tr_dispho_preplot_tof_2t_eg2018B'
#dostack(hist_list_preplot_tof_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_preplot/tr_dispho_preplotNot_2t_eg2018B'
#dostack(hist_list_preplotNot_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B_preplot/tr_dispho_preplotStc_2t_eg2018B'
#dostack(hist_list_preplotStc_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B/tr_dispho_ls_2t_eg2018B'
#dostack(hist_list_same_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B/tr_dispho_2tStc_eg2018B'
#dostack(hist_list_2tStc_2018B, outname, date, loc2_layout)

#outname = 'downloads/dispho_2t_eg2018B/tr_dispho_2t_eg2018B'
#dostack(hist_list_2t_2018B, outname, date, loc2_layout)

#outname = 'downloads/rmax_cut_plots_v3/tr_R2017Bpart_raw94_rmaxcuts'
#dostack(hist_list_rmax_v3, outname, date, loc2_layout)

#-----------------------------
#outname = 'downloads/Run2016degtof/tr_R2016degtof_gll_Emin10_FulldR'
#dostack(hist_list16deg, outname, date, gll_layout)

#outname = 'downloads/Run2016deg/tr_mu_R2016deg_gll_Emin10_FulldR'
#dostack(hist_list16deg_mu, outname, date, gll_layout)

#outname = 'downloads/Run2016raw94tof/tr_R2016rawtof_gll_Emin10_FulldR'
#dostack(hist_list16r, outname, date, gll_layout)

#outname = 'downloads/Run2016raw94/tr_mu_R2016raw_gll_Emin10_FulldR'
#dostack(hist_list16r_mu, outname, date, gll_layout)

#------------------------------
#outname = 'downloads/Run2017raw94tof/tr_R2017rawtof_gll_Emin10_FulldR'
#dostack(hist_list17r, outname, date, gll_layout)

#outname = 'downloads/Run2017degtof/tr_R2016degtof_gll_Emin10_FulldR'
#dostack(hist_list17deg, outname, date, gll_layout)

#outname = 'downloads/Run2017deg/tr_mu_R2017deg_gll_Emin10_FulldR'
#dostack(hist_list17deg_mu, outname, date, gll_layout)

#outname = 'downloads/Run2017raw94/tr_R2017raw_gll_Emin10_FulldR'
#dostack(hist_list17r, outname, date, gll_layout)

#outname = 'downloads/Run2017raw94/tr_mu_R2017raw_gll_Emin10_FulldR'
#dostack(hist_list17r_mu, outname, date, gll_layout)

#------------------------------
#outname = 'downloads/Run2018Arawtof/tr_R2018rawtof_gll_Emin10_FulldR'
#dostack(hist_list18r, outname, date, gll_layout)

#outname = 'downloads/Run2018Aminitof/tr_R2018egminitof_gll_Emin10_FulldR'
#dostack(hist_list18mini, outname, date, gll_layout)

#outname = 'downloads/Run2018raw/tr_mu_R2018raw_gll_Emin10_FulldR'
#dostack(hist_list18r_mu, outname, date, gll_layout)

#outname = 'downloads/Run2018egmini/tr_mu_R2018egmini_gll_Emin10_FulldR'
#dostack(hist_list18mini_mu, outname, date, gll_layout)

#----------------------------
#outname = 'downloads/tr_eg_tof_notracker_global_Emin10_FulldR'
#dostack(hist_list_glob, outname, date, glob_layout)

#outname = 'downloads/tr_eg_tof_local_Emin10_FulldR'
#dostack(hist_list_loc, outname, date, loc_layout)

#outname = 'downloads/tr_global_Emin10_ltdR'
#dostack(hist_list_glob_lt, outname, date, glob_layout)

#outname = 'downloads/tr_global_Emin10_gtdR'
#dostack(hist_list_glob_gt, outname, date, glob_layout)


#--------------------------------
#outname = 'downloads/03April19_Run2016/tr_R2016_gll_Emin10'
#dostack(hist_list16, outname, date, gll_layout)

#outname = 'downloads/03April19_Run2016/tr_R2016_lv_Emin10'
#dostack(hist_list16lv, outname, date, loc_layout)

#outname = 'downloads/03April19_Run2016/tr_R2016_gv_Emin10'
#dostack(hist_list16gv, outname, date, glob_layout)

#outname = 'downloads/06Feb19_Run2017/tr_R2017_gll_Emin10'
#dostack(hist_list17, outname, date, gll_layout)

#outname = 'downloads/06Feb19_Run2017/tr_R2017_lv_Emin10'
#dostack(hist_list17lv, outname, date, layout)

#outname = 'downloads/06Feb19_Run2017/tr_R2017_gv_Emin10'
#dostack(hist_list17gv, outname, date, layout)

#outname = 'downloads/08April19_Run2018/tr_R2018_gll_Emin10'
#dostack(hist_list18, outname, date, gll_layout)

#outname = 'downloads/08April19_Run2018/tr_R2018_lv_Emin10'
#dostack(hist_list18lv, outname, date, loc_layout)

#outname = 'downloads/08April19_Run2018/tr_R2018_gv_Emin10'
#dostack(hist_list18gv, outname, date, glob_layout)

#outname = 'downloads/08April19_Run2018/tr_R2018_lt_Emin10'
#dostack(hist_list18lt, outname, date, loc_layout)

#outname = 'downloads/08April19_Run2018/tr_R2018_gt_Emin10'
#dostack(hist_list18gt, outname, date, glob_layout)


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


