from ROOT import *

can_x_width = 700
can_y_width = 600

def fixOverlay():
        gPad.RedrawAxis()


def setTDRStyle():
        tdrStyle = TStyle("tdrStyle","Style for P-TDR")

##// For the canvas:
        tdrStyle.SetCanvasBorderMode(0)
        tdrStyle.SetCanvasColor(kWhite)
        tdrStyle.SetCanvasDefH(can_x_width)  #//Height of canvas
        tdrStyle.SetCanvasDefW(can_y_width)  #//Width of canvas
        tdrStyle.SetCanvasDefX(0)    #//POsition on screen
        tdrStyle.SetCanvasDefY(0)

#// For the Pad:
        tdrStyle.SetPadBorderMode(0)
        #// tdrStyle.SetPadBorderSize(Width_t size = 1) 
        tdrStyle.SetPadColor(kWhite)
        tdrStyle.SetPadGridX(False)
        tdrStyle.SetPadGridY(False)
        tdrStyle.SetGridColor(0)
        tdrStyle.SetGridStyle(3)
        tdrStyle.SetGridWidth(1)

#// For the frame:
        tdrStyle.SetFrameBorderMode(0)
        tdrStyle.SetFrameBorderSize(1)
        tdrStyle.SetFrameFillColor(0)
        tdrStyle.SetFrameFillStyle(0)
        tdrStyle.SetFrameLineColor(1)
        tdrStyle.SetFrameLineStyle(1)
        tdrStyle.SetFrameLineWidth(1)


#// set the paper & margin sizes
        tdrStyle.SetPaperSize(20,26);
        tdrStyle.SetPadTopMargin(0.09);
        tdrStyle.SetPadRightMargin(0.15);
        tdrStyle.SetPadBottomMargin(0.12);
        tdrStyle.SetPadLeftMargin(0.15);

#// For the histo:
        #tdrStyle.SetHistFillColor(1) 
        #tdrStyle.SetHistFillStyle(0) 
        tdrStyle.SetHistLineColor(1)
        tdrStyle.SetHistLineStyle(0)
        tdrStyle.SetHistLineWidth(1)
        #tdrStyle.SetLegoInnerR(Float_t rad = 0.5) 
        #tdrStyle.SetNumberContours(Int_t number = 20) 
        tdrStyle.SetMarkerStyle(8)
        tdrStyle.SetMarkerSize(0.8)
        tdrStyle.SetEndErrorSize(2)
        #// tdrStyle.SetErrorMarker(20) 
        tdrStyle.SetErrorX(0.5)

#//For the fit/function:
        tdrStyle.SetOptFit(0)
        tdrStyle.SetFitFormat("5.4g")
        tdrStyle.SetFuncColor(2)
        tdrStyle.SetFuncStyle(1)
        tdrStyle.SetFuncWidth(1)

#//For the date:
        tdrStyle.SetOptDate(0)
        #// tdrStyle.SetDateX(Float_t x = 0.01) 
        #// tdrStyle.SetDateY(Float_t y = 0.01) 

#// For the statistics box:
        tdrStyle.SetOptFile(0)
        tdrStyle.SetOptStat(0)  #// To display the mean and RMS:   SetOptStat("mr");
        tdrStyle.SetStatColor(kWhite)
        tdrStyle.SetStatFont(42)
        tdrStyle.SetStatFontSize(0.025)
        tdrStyle.SetStatTextColor(1)
        tdrStyle.SetStatFormat("6.4g")
        tdrStyle.SetStatBorderSize(1)
        tdrStyle.SetStatH(0.1)
        tdrStyle.SetStatW(0.15)
        #// tdrStyle.SetStatStyle(Style_t style = 1001) 
        #tdrStyle.SetStatX( 0.825 ) 
        #tdrStyle.SetStatY( 0.85 ) 


#// For the Global title:
        tdrStyle.SetOptTitle(0)
        tdrStyle.SetTitleFont(42)
        tdrStyle.SetTitleColor(1)
        tdrStyle.SetTitleTextColor(1)
        tdrStyle.SetTitleFillColor(10)
        tdrStyle.SetTitleFontSize(0.065)
        #// tdrStyle.SetTitleH(0)  // Set the height of the title box
        #// tdrStyle.SetTitleW(0)  // Set the width of the title box
        #tdrStyle.SetTitleX(0.5)  #// Set the position of the title box
        #tdrStyle.SetTitleAlign(23) 
        #tdrStyle.SetTitleY(0.985)  #// Set the position of the title box
        #//  tdrStyle.SetTitleStyle(1001) 
        #tdrStyle.SetTitleBorderSize(0) 

#// For the axis titles:

        tdrStyle.SetTitleColor(1, "XYZ")
        tdrStyle.SetTitleFont(42, "XYZ")
        tdrStyle.SetTitleSize(0.05, "XYZ")
        #tdrStyle.SetTitleXSize(0.06)  #// Another way to set the size?  0.02
        #tdrStyle.SetTitleYSize(0.06) 
        #tdrStyle.SetTitleXOffset(0.9)
        tdrStyle.SetTitleXOffset(1.0)
        #tdrStyle.SetTitleYOffset(1.1)
        tdrStyle.SetTitleYOffset(1.4)
        #//  tdrStyle.SetTitleOffset(1.1, "Y")  // Another way to set the Offset

#// For the axis labels:

        tdrStyle.SetLabelColor(1, "XYZ")
        tdrStyle.SetLabelFont(42, "XYZ")
        tdrStyle.SetLabelOffset(0.007, "XYZ")
        #tdrStyle.SetLabelSize(0.05, "XYZ")
        tdrStyle.SetLabelSize(0.03, "XYZ")

#// For the axis:

        tdrStyle.SetAxisColor(1, "XYZ")
        tdrStyle.SetStripDecimals(kTRUE)
        tdrStyle.SetTickLength(0.03, "XYZ")
        tdrStyle.SetNdivisions(505, "X")
        tdrStyle.SetPadTickX(1)   #// To get tick marks on the opposite side of the frame
        tdrStyle.SetPadTickY(1)

#// Change for log plots:
        #tdrStyle.SetOptLogx(0)
        #tdrStyle.SetOptLogy(0)
        #tdrStyle.SetOptLogz(0)

#// Postscript options:
        #//#tdrStyle.SetPaperSize(20.,20.) 
        #// tdrStyle.SetLineScalePS(Float_t scale = 3) 
        #// tdrStyle.SetLineStyleString(Int_t i, const char* text) 
        #// tdrStyle.SetHeaderPS(const char* header) 
        #// tdrStyle.SetTitlePS(const char* pstitle) 

        #// tdrStyle.SetBarOffset(Float_t baroff = 0.5) 
        #// tdrStyle.SetBarWidth(Float_t barwidth = 0.5) 
        #// tdrStyle.SetPaintTextFormat(const char* format = "g") 
        #// tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0) 
        #// tdrStyle.SetTimeOffset(Double_t toffset) 
        #// tdrStyle.SetHistMinimumZero(kTRUE) 

        #//tdrStyle.SetHatchesLineWidth(5) 
        #//tdrStyle.SetHatchesSpacing(0.05)

        #// color palette
        tdrStyle.SetPalette(kBird)
        tdrStyle.SetNumberContours(255)

        tdrStyle.cd()
        gROOT.ForceStyle()
