from ROOT import *
from math import *
import time
from jwk_tdr_style_py import *

setTDRStyle()
gStyle.SetPadTickY(1)
gStyle.SetPadTickX(1)
gStyle.SetOptStat()
#if layout['logx'] : gStyle.SetOptLogx(1)
#if layout['logy'] : gStyle.SetOptLogy(1)
c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
c1.cd()
#if layout['logx'] : c1.SetLogx()
#if layout['logy'] : 
#c1.SetLogy()
c1.SetLogz()
c1.SetGridx(1)
c1.SetGridy(1)
#c1.Update()
#c1.Draw()
#legend = TLegend(l[0],l[1],l[2],l[3]);
#legend.SetName('legend')
gStyle.SetLegendFont(42)
gStyle.SetLegendTextSize(0.03)
#legend.SetBorderSize(0)

#lat = TLatex()
#lat.SetNDC()

#k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
k = [kMagenta+2,kBlue+1,kAzure+4,kBlack,kYellow+1,kViolet+7,kOrange+7,kRed+2,kGreen+3, kGray]

x_range = [-26, 26 ]
x_bin = 520
#x_label = 'ncjitter time [ns]'
#h1 = TH1F('hpx', 'adjTime', 200, -10, 10)

est = 0
#edv = 15
#eed = 600
#edv = 10
#eed = 200
edv = 10
eed = 100
##############
#est = 400
#edv = 16
#eed = 2000
############3
he = TH1F('hpe', 'RH Energy; Energy [GeV];a.u.', edv, est, eed )
hr1 = TH1F('hr1', 'kOOT % CC;Energy [GeV];kOOT %', edv, est, eed )
hs = TH1F('hps', 'kOOT Same', edv, est, eed )
hr2 = TH1F('hr2', 'kOOT % Rt;Energy [GeV];kOOT %', edv, est, eed )
#hd = TH1F('hpd', 'kOOT Denom', edv, est, eed )

#hn = TH1F('hpn', 'kOOT Ratio;Energy [GeV];kOOT %', 40, 0, 200)
#hd = TH1F('hpd', 'kOOT Denom', 40, 0, 200)
#hn = TH1F('hpn', 'kOOT Ratio;Energy [GeV];kOOT %', 40, 0, 20)
#hd = TH1F('hpd', 'kOOT Denom', 40, 0, 20)

#h2 = TH2F('h2px', 'time v energy', 520, -26, 26, 375, 0, 750 );
#h2 = TH2F('h2px', 'time v energy', 200, -10, 10, 50, 0, 25 );
#h2 = TH2F('h2px', 'time v amp', 300, -15, 15, 300, 0, 600 );
#h2 = TH2F('h2px', 'time v energy', 300, -15, 15, 50, 0, 25 );
#h2 = TH2F('h2px', 'energy v adcToGeV0', 50, 0, 500, 120, 0, 12.0 );
h2 = TH2F('h2px', 'CC v RT', 53, -26.5, 26.5, 53, -26.5, 26.5 );

#h1 = TH1F('hpx1', 'adjTimeCC', 520, -26, 26)
h1 = TH1F('hpx1', 'amplitude', 100, 0, 200)
#h1.GetXaxis().SetTitle("jitter + itc [ns]") 
#h1.GetXaxis().SetTitle("itime constant [ns]")
#h1.GetXaxis().SetTitle(" kOOT False corTime [ns]")
#h1.GetXaxis().SetTitle("corTime (threshold) [ns]")
h1.GetXaxis().SetTitle("Amplitude")
#h2 = TH1F('hpx2', 'adjTimeRT', 520, -26, 26)
#h2.GetXaxis().SetTitle("corTime (threshold) [ns]")
#h2.GetYaxis().SetTitle("Energy [GeV]")
h2.GetXaxis().SetTitle("CC time [ns]")
h2.GetYaxis().SetTitle("Rt time [ns]")
#outname = 'jitterTime'
#h1.SetTitle(outname)


#inputfile = 'koot_encoding_log_full.txt'
#inputfile = 'koot_eb_encoding_log.txt'
#inputfile = 'koot_eb_366850_log.txt'
#inputfile = 'koot_eb_ratio_366850_log.txt'

#inputfile = 'koot_eb_cc_amp_366850_log.txt'
#inputfile = 'koot_10k_eb_uncc_amp_366850_log.txt'
#inputfile = 'koot_10k_eb_crcc_amp_366850_log.txt'
#inputfile = 'koot_10k_eb_rt_amp_366850_log_v3.txt'

#inputfile1 = 'koot_10k_eb_crcc_energy_366850_log.txt'

#inputfile1 = 'koot_10k_eb_uncc_energy_366850_log.txt'
#inputfile1 = 'koot_10k_eb_uncc_t3_energy_366850_log.txt'
#inputfile1 = 'koot_10k_eb_uncc_t4_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_energy_366850_log.txt'

#inputfile1 = 'koot_47k_eb_uncc_t25_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t35_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t40_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t45_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c6n225_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c6n255_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c7n255_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c8n255_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c7n225_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c8n225_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c75n255_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c9n255_energy_366850_log.txt'
#inputfile1 = 'koot_47k_eb_uncc_t30c85n255_energy_366850_log.txt'
inputfile1 = 'koot_47k_eb_uncc_t25c6n285_energy_366850_log.txt'

#inputfile2 = 'koot_10k_eb_rt_energy_366850_log.txt'
inputfile2 = 'koot_47k_eb_rt_t50_energy_366850_log.txt'

rhcnt = 0
tcnt = 0

with open( inputfile1, 'r' ) as file1, open( inputfile2, 'r' ) as file2 :
  for line1, line2 in zip( file1, file2 ):

    info1 = line1.split()	
    jitter1 = float(info1[0])
    nocorjitter1 = float(info1[1])
    cortime1 = float(info1[2])   
    koot1 = int(info1[3])
    sigthr1 = float(info1[4]) 
    sigmat1 = float(info1[5])
    thres1 = float(info1[6])
    nterm1 = float(info1[7])
    cterm1 = float(info1[8])
    energy1 = float(info1[9]) #/float(info1[13])
    lasercalib1 = float(info1[10])
    interCalib1 = float(info1[11])
    adcToGeV01 = float(info1[12])
    pedrms121 = float(info1[13])
    itimeconst1 = float(info1[14])
    offsetTime1 = float(info1[15])
    amplitude1 = float(info1[16])
    rhid1 = float(info1[17])


    info2 = line2.split()
    jitter2 = float(info2[0])
    nocorjitter2 = float(info2[1])
    cortime2 = float(info2[2])
    koot2 = int(info2[3])
    sigthr2 = float(info2[4])
    sigmat2 = float(info2[5])
    thres2 = float(info2[6])
    nterm2 = float(info2[7])
    cterm2 = float(info2[8])
    energy2 = float(info2[9]) # /float(info2[13])
    lasercalib2 = float(info2[10])
    interCalib2 = float(info2[11])
    adcToGeV02 = float(info2[12])
    pedrms122 = float(info2[13])

    rhcnt += 1
    #if koot1 and energy1 > 10 and energy1 > 120: 
    if energy1 == energy2 :
        tcnt += 1

    #if energy2 < 20 and cortime2 < -2 and cortime2 > -4 : he.Fill(energy2);
    #if energy2 < 5 and cortime2 < -2 and cortime2 > -4 : print( info2 )
    #if jitter1 != 0 : print( info1 )

    if jitter2 != 0.0 and energy1 == energy2 : 
        energy = energy1/pedrms121
        he.Fill(energy)
        if( koot1 == 1 ) : hr1.Fill(energy)
        if( koot2 == 1 ) : hr2.Fill(energy)
        if( koot1 == koot2 ) : hs.Fill(energy)

    #h2.Fill( energy1, pedrms122 )

    #testtime = nocorjitter
    #testtime = jitter1 + itimeconst1
    #testtime = itimeconst1
    #testtime = cortime1
    #if testtime < -2 and energy > 12 and energy < 18  and koot == 1 : print( info )
    #print( jitter, ncorjitter, co2ime, koot, itc, osc )
    #if koot1 == 0 : #and testtime > -20 :
    if jitter2 == 0 and amplitude1 > 100 : # and energy1 > 2.0 : # and koot1 == 0 :
    ##if amp > thres :
    #if abs(testtime) < 50 :
    #    #print( testtime, thres, amp )
    #     h1.Fill( testtime )
    #     h1.Fill( pedrms121 )
    #     h1.Fill( amplitude1 )
          h2.Fill( cortime1, cortime2 )
    #     hg.Fill( testtime, energy1 )


print( "precnt rh koot true E > 120 : ", float(tcnt)/rhcnt )
#print( "precnt rh koot true 10 > E < 120 : ", float(tcnt)/rhcnt )

#he.UseCurrentStyle()
#he.SetMarkerStyle(1+25)
#he.SetLineColor(k[1])
#he.SetMarkerColor(k[1])
#legend.AddEntry(h1[n],lego,'epl')
#legend.Draw('same') 

#h1.Draw('ep')
h2.Draw('colz')

#hr1.SetMinimum(0)
#hr1.SetMaximum(1.2)
#mg.GetXaxis().SetRangeUser(x[0],x[1])

'''
te1 = TEfficiency( hr1, he )
te1.UseCurrentStyle()
te1.SetTitle("kOOT Ratio;Energy [GeV];kOOT %")
te1.SetMarkerStyle(1+25)
te1.SetLineColor(k[1])
te1.SetMarkerColor(k[1])
te1.Draw()
te2 = TEfficiency( hr2, he )
te2.UseCurrentStyle()
te2.SetTitle("kOOT Ratio;Energy [GeV];kOOT %")
te2.SetMarkerStyle(2+25)
te2.SetLineColor(k[0])
te2.SetMarkerColor(k[0])
te2.Draw('same')
te3 = TEfficiency( hs, he )
te3.UseCurrentStyle()
te3.SetTitle("kOOT Ratio;Energy [GeV];kOOT %")
te3.SetMarkerStyle(3+25)
te3.SetLineColor(k[3])
te3.SetMarkerColor(k[3])
te3.Draw('same')
'''

outname = 'koot_c85n255_t30'

#he.UseCurrentStyle()
#he.SetMarkerStyle(1+25)
#he.SetLineColor(k[1])
#he.SetMarkerColor(k[1])
#he.Draw('')
#outname = 'rhenergy'

#h2.Draw('colz')
#he.Draw( 'EP')

c1.Modified()
c1.Update()

'''
teg = te1.GetPaintedGraph();
teg.SetMinimum(0)
teg.SetMaximum(1.2)
c1.Update()
'''

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')
#c1.Print( outname + '_' + date + '.pdf' )
c1.Print( outname + '_' + date + '.png' )
#c1.SaveAs( outname + '_' + date + '.root' )
#c1.SaveAs( outname + '_' + date + '.C' )
#c1.Show()
c1.Close()



