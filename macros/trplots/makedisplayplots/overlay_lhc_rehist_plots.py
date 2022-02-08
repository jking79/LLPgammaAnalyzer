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
    dofit = False
    #dofit = False
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
    gStyle.SetTitleYOffset(1.0)
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
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',75,950,2)
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',125,1275,2)
        hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',75,2250,3)
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',6,100,2)
        hfit.SetParName(0,'N')
        hfit.SetParameter(0,50)
        #hfit.SetParLimits(0,10,50)
        hfit.SetParLimits(0,0,100)
        #hfit.SetParameter(0,5)
        #hfit.SetParLimits(0,0,10)
        hfit.SetParName(1,'C')
        hfit.SetParameter(1,0.5)
        hfit.SetParLimits(1,0,1)
        #hfit.SetParLimits(1,0.02,0.3)
        #hfit.SetParLimits(1,0.02,1.0)
        #hfit.SetParameter(1,0.05)
        #hfit.SetParLimits(1,0.001,1.0)
        hfit.SetParName(2,'S')
        hfit.SetParameter(2,0)
        hfit.SetParLimits(2,-20,20)

    mg = TMultiGraph();

    #thebinmid.append([  88.01, 111.9, 136.4, 161.4, 196.6, 247.8, 298.5, 348.9, 419.7, 528.9, 656.9, 801.5,   990, 1487.5,  1975])
    #thebinerror.append([ 7.09, 7.156, 7.154, 7.172, 14.28, 14.39, 14.43,  14.4, 28.53, 35.37, 40.58, 43.16, 38.06,      0,     0])
    #thebinmid.append([  90.43, 112.6, 137.1, 161.8, 196.8, 247.1, 297.6, 347.9, 417.6, 528.2,   663, 825.1,  1022,   1320,  1975])
    #thebinerror.append([6.521, 7.156, 7.192, 7.192, 14.22, 14.29, 14.33, 14.34,  28.4, 35.39, 42.34, 53.72, 63.78,  51.57,     0])
    #thebinmid.append([  94.04,   115, 137.7, 162.2, 198.1, 247.7, 297.7, 347.9,   418, 528.4, 664.2, 833.9,  1064,   1387,  1779])
    #thebinerror.append([5.056, 6.691, 7.158, 7.201, 14.32, 14.31, 14.33, 14.36, 28.43, 35.41, 42.68, 56.53,  84.3,  81.16, 65.90])



    for histname, tree, infile, lego  in hist_list :
   
        f1.append(TFile.Open(infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
        print('Opening', infile, 'at', f1[n]) 

        rt5sum = 0;
        rt5cnt = 0;
        rt5calihist = f1[n].Get('ratio_e5')
        #rt5calihist = f1[n].Get('Hist2')
        for bn in range(0,3564) :
        #for bn in range(3,51) :
            thisval = rt5calihist.GetBinContent(bn)
            if thisval != 0.0 :
                rt5sum = rt5sum + thisval
                rt5cnt = rt5cnt + 1
        rt5cali = rt5sum/rt5cnt

        cc5sum = 0;
        cc5cnt = 0;
        cc5calihist = f1[n].Get('kucc_e5')
        #cc5calihist = f1[n].Get('rfcHist2')
        for bn in range(0,3564) :
        #for bn in range(3,51) :
            thisval = cc5calihist.GetBinContent(bn)
            if thisval != 0.0 :
                cc5sum = cc5sum + thisval
                cc5cnt = cc5cnt + 1
        cc5cali = cc5sum/cc5cnt

        orighist = f1[n].Get(hist)
        print('Getting', hist, 'at', orighist)
        nbins = 3564 # 3564 orighist.GetNbinsX()
        bstep = 1
        #nbins = 60
        #bstep = 3564/60
        #startbin = 1450
        startbin = 0
        #endbin = 1650 
        #endbin = 1520
        endbin = nbins
        #low = 0 #orighist.GetMinimumBin()
        #high = nbins #orighist.GetMaximumBin()
        bscale = 1 #(high-low)/(nbins)
        print('Find Bins',endbin-startbin,startbin,endbin,bscale) 
        htitle = 'hist' + str(n)
        #mybins = array([0,2,4,6,8,10,12,14,16,20,24,30,38,52,68,100],dtype='float64')
        #mybins = [0,2,4,6,8,10,12,14,16,20,24,30,38,52,68,100]
        #binwidths = [1,1,1,1, 1, 1, 1, 1, 2, 2, 3, 4, 7, 8, 16]
        #numExcludeBins = 0
        ##mybins = array([75,  100, 125, 150, 175, 225, 275, 325, 375, 475, 600, 750,  950, 1275, 1700, 2250],dtype='float64')
        #mybins = [4, 6, 7, 8, 9, 12, 20]
        #binwidths = [1,0.5,0.5,0.5,1.5,4]
        #mybins = [75,  100, 125, 150, 175, 225, 275, 325, 375, 475, 600, 750,  950, 1275, 1700, 2250]
        #mybins = [75,  100, 125, 150, 175, 225, 275, 325, 375, 475, 600, 950, 2250]
        #mybins = [0, 75,  100, 125, 150, 175, 225, 275, 325, 375, 475, 600, 950, 2250]
        #binwidths = [   12.5,12.5,12.5,12.5,25.0,25.0,25.0,25.0,50.0,62.5,75.0,100.0,162.5,212.5,275.0]
        #binwidths = [ 12.5,12.5,12.5,12.5,25.0,25.0,25.0,25.0,50.0,62.5,175.0,650.0]
        #binwidths = [ 37.5,12.5,12.5,12.5,12.5,25.0,25.0,25.0,25.0,50.0,62.5,175.0,650.0]
        #mybins = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
        #mybins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
        #binwidths = [0.5,0.5,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25]
        #numExcludeBins = 0

        #binwidths = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]        
        lenmybins = endbin - startbin
        #lenmybins = len(mybins) - numExcludeBins
        #myhist = TH1D(htitle,"",150,low,high)
        #h1.append(TGraphErrors(lenmybins))
        h1.append(TGraphErrors(lenmybins))
        #h1.append(TH2D(htitle,"",(10*nbins),low,high,240,0,2.4))
        #h1.append(TH1F(htitle,"",len(mybins)-1,mybins))

        bincnt = 0
        sbinval = 0
        sbinerr = 0
        sbinmid = 0
        for bn in range(1,lenmybins):
            doplot = False
            #binval = float(orighist.GetBinContent(bn))*binwidths[bn-1]*2
            binval = float(orighist.GetBinContent(bn+startbin))
            #binerr = 0.0
            binerr = float(orighist.GetBinError(bn+startbin))
            binmid = float(orighist.GetBinCenter(bn+startbin))
            #binmid = thebinmid[n][bn-1]  #float(orighist.GetBinCenter(bn))
            #if binval > 0.0 : 
            bincnt = bincnt + 1    
            sbinval = sbinval + binval
            sbinerr = sbinerr + binerr*binerr
            sbinmid = sbinmid + binmid
            #if( bincnt == 50 ) : 
            if( False ) :
                bincnt = 0
                doplot = True
                binval = sbinval/50
                binerr = sqrt(sbinerr)/50
                binmid = sbinmid/50
                sbinval = 0
                sbinerr = 0
                sbinmid = 0
            if abs(binerr) < 0.1 and binval != 0.0 :
            #if abs(binerr) < 0.1 and binval != 0.0 and doplot == True :
            #if abs(binerr) < 0.1 and binval != 0.0 and binmid > 2:
                if( sxtal ):  
                    print( '--- adjusting bin value' )
                    binval = binval/sqrt(2)
                doplot = False
                if 'ratio' in hist : binval = binval - rt5cali
                else : binval = binval - cc5cali
                #if 'rfc' in hist : binval = binval - cc5cali
                #else : binval = binval - rt5cali
                binmid = binmid*bstep
                h1[n].SetPoint(bn,binmid,binval)
                #ebin = int(binmid/bscale)+1
                widtherr = 0;
                #widtherr = binwidths[bn-1]/sqrt(3)
                #widtherr = (0.25)/sqrt(3)
                #widtherr = thebinerror[n][bn-1] #binwidths[bn-1]/sqrt(3)
                #widtherr = bscale/sqrt(3)
                h1[n].SetPointError(bn,widtherr,binerr)
                #print('Fill bin',bn,'at',binmid,'with',binval,'err',widtherr,'by',binerr,'for:',mybins[bn-1],'to',mybins[bn],'width',binwidths[bn-1]) 
                print('Fill bin',bn,'at',binmid,'with',binval,'error',binerr,'by',widtherr)  

#        h1.append(f1[n].Get(hist))

#       print( f1 )
#       print( h1 )
       
        print('Setting markers')
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
        #k = [kGreen+2,kYellow+1,kBlue+2,kBlack,kYellow+7,kGreen+3,kBlue+4,kGreen+2,kAzure+4,kGreen+2]
        #k = [kBlue+4,kBlue+1,kCyan+2,kGreen+3,kGreen-3,kYellow+3,kYellow+1]
        k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3]
        #k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kGreen-7,kGreen+2,kBlue-7,kBlue+2,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kAzure-3,kAzure-6,kSpring-6,kSpring-7]
        #k = [kBlack]
#       # = [kMagenta+2,kBlue+2,kGreen+2]
        h1[n].SetLineColor(k[n])
        if dofit : hfit.SetLineColor(k[n])
        #if dofit : hfit.SetLineColor(kRed+2)
        h1[n].SetMarkerColor(k[n])

        #msz = 0.1
        msz = 0.5
        #msz = 0.8
        #msz = 1.0
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.8)
        elif( n == 3 ) : h1[n].SetMarkerSize(msz+0.8)
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
        print('Checking for fits')
        if( first ) :
                #h1[n].Draw('ep')
                #h1[n].Draw('AP')
                if dofit : h1[n].Fit(hfit,'RE')
                #if dofit : h1[n].Fit(hfit,'REQ')
                first = False
        else :
                #h1[n].Draw('epSAME')
                #h1[n].Draw('SAME AP')
                if dofit : h1[n].Fit(hfit,'RE+')
                #if dofit : h1[n].Fit(hfit,'REQ+') 
        #c1.Update()
        #c1.Draw()

        if dofit : 
                 paramn.append(str(hfit.GetParameter(0)))
                 paramc.append(str(hfit.GetParameter(1)))
                 params.append(str(abs(hfit.GetParameter(2))))
                 pne = hfit.GetParError(0)
                 pce = hfit.GetParError(1)
                 pse = hfit.GetParError(2)
                 #print('Fit info',paramn[n],pne,paramc[n],pce)
                 print('Fit info',paramn[n],pne,paramc[n],pce,params[n],pse)
                 if pne < 0.01 : pne = 0.01
                 if pce < 0.0001 : pce = 0.0001
                 parnerror.append(str(pne))
                 parcerror.append(str(pce))
                 parserror.append(str(pse))
                 #print( 'C: ' + param + ' +/- ' + parerror )
                 #lat_param = '#color['+str(k[n])+']{N : ' + paramn[0:4] + ' #pm ' + parnerror[0:4] + ' [ns]  C : ' + paramc[0:6] + ' #pm ' + parcerror[0:6] + ' [ns]}'
                 #lat.SetTextSize(0.03);
                 #lat.DrawLatex(t[3],t[4]-n*.035,lat_param);    
                 #c1.Modified()
                 #c1.Update()       
 
        print('Adding legend')
        legend.AddEntry(h1[n],lego,'epl');
        #c1.Update()
        #c1.Draw()
        n += 1

        #End of loop
    print('Building Multi Hist')
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

    print('Drawing on Multi Hist')
    #mg.Draw('AP')
    if lego != 'none' : legend.Draw('same')  #   legend inclusion switch

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
    #lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + 2C^{2}'
    lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + #frac{S^{2}}{A_{eff}/#sigma_{n}} + 2C^{2}'
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
            #lat_param ='#color['+str(k[l])+']{N : '+paramn[l][0:4]+' #pm '+parnerror[l][0:4]+' [ns]  C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' }'
            lat_param =	'#color['+str(k[l])+']{'
            lat_param = lat_param + 'S : '+params[l][0:4]+' #pm '+parserror[l][0:4]+' [ns] '
            lat_param = lat_param + 'N : '+paramn[l][0:4]+' #pm '+parnerror[l][0:4]+' [ns] '
            lat_param = lat_param + 'C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' [ns]}'
            lat.SetTextSize(0.03);
            lat.DrawLatex(t[3],t[4]-l*.035,lat_param);

    
    if layout['logx'] : c1.SetLogx()
    ####c1.BuildLegend()
    c1.Modified()
    c1.Update()

    print('Saving output')
    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.SaveAs( outname + '_' + date + '.C' )
    #c1.Show()
    c1.Close()

#from overlay_hist_defs_v2 import *
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
#xtitle = ''
#xtitle = 'A_{eff}/#sigma_{n}'
xtitle = 'BX'
#xtitle = 'E_{eff} [GeV]'
#xtitle = 'pCaloHit E_{eff} [GeV]'
#xtitle = 'pCaloHit energy [GeV]'
#xtitle = 'GeV'
ytitle = 'Adj. Mean BX Time [ns]'
#ytitle = 'Mean BX Time [ns]'
#ytitle = '#sigma(t_{1}-t_{2}) [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
#ytitle = ''
#htitle = 'A_{eff}/#sigma_{n} vs #sigma_{#delta_{t}}'
htitle = ''
#islogx = True
islogx = False
#islogy = True
islogy = False

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
#---------------------------------------------------------------

ptitle=[' 320804-320809','','#splitline{EPEM}{}'] #{GT 106X_dataRun2_v28}']
y = [ 12.0, -2.0 ]
#y = [ 2.75, -1.25 ]
x = [ 1450, 1650 ]
l = [ 0.7,0.65,0.925,0.9 ]
t = [0.2,0.825,0.0,0.175,0.23]
#t = [0.325,0.82,0.0,0.175,0.28]
outname = 'downloads/tr_lhc_18A_v25_talk'
dostack(hl_lhc_v25_ratio, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=[' 320804-320809','','#splitline{EPEM}{E > 2 GeV}'] #{GT 106X_dataRun2_v28}']
#y = [ 12.0, -2.0 ]
y = [ 12.0, -2.0 ]
x = [ 1450, 1650 ]
l = [ 0.7,0.65,0.925,0.9 ]
t = [0.2,0.825,0.0,0.175,0.23]
#t = [0.325,0.82,0.0,0.175,0.28]
outname = 'downloads/tr_lhc_18A_ratio2_note'
#dostack(hl_lhc_v25_ratio, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=[' 320804-320809','','#splitline{EPEM}{}'] #{GT 106X_dataRun2_v28}']
#y = [ 12.0, -2.0 ]
y = [ 12.0, -2.0 ]
x = [ 1460, 1520 ]
l = [ 0.7,0.65,0.925,0.9 ]
t = [0.2,0.825,0.0,0.175,0.23]
#t = [0.325,0.82,0.0,0.175,0.28]
outname = 'downloads/tr_lhc_18A_ratioall_note'
#dostack(hl_lhc_v25_ratio, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=[' 320804-320809','','#splitline{EPEM}{}'] #{GT 106X_dataRun2_v28}']
#y = [ 12.0, -2.0 ]
y = [ 13.0, -2.0 ]
x = [ 0, 3600 ]
l = [ 0.7,0.65,0.925,0.9 ]
t = [0.2,0.825,0.0,0.175,0.23]
#t = [0.325,0.82,0.0,0.175,0.28]
outname = 'downloads/tr_lhc_18A_ratioall_full50bin_note'
#dostack(hl_lhc_v25_ratio, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=[' 320804-320809','','#splitline{EPEM}{Train Position}'] #{GT 106X_dataRun2_v28}']
#y = [ 12.0, -2.0 ]
y = [ 0.475, 0.34 ]
x = [ 20, 3350 ]
l = [ 0.7,0.65,0.925,0.9 ]
t = [0.2,0.825,0.0,0.175,0.23]
#t = [0.325,0.82,0.0,0.175,0.28]
outname = 'downloads/tr_lhc_18A_ratioall_stp_note'
#dostack(hl_lhc_v25_st_ratio, outname, date, Ic_layout, ptitle,  y, x, l, t)

ptitle=[' 320688','','#splitline{EPEM}{Train Position}'] #{GT 106X_dataRun2_v28}']
#y = [ 12.0, -2.0 ]
y = [ 0.06, -0.04 ]
x = [ 20, 3350 ]
l = [ 0.7,0.65,0.925,0.9 ]
t = [0.2,0.825,0.0,0.175,0.23]
#t = [0.325,0.82,0.0,0.175,0.28]
outname = 'downloads/tr_lhc_18A_rfc_stp_note'
#dostack(hl_lhc_rfc_st, outname, date, Ic_layout, ptitle,  y, x, l, t)




