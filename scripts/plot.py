# Implementation of the plotting step of the analysis
#
# The plotting combines the histograms to plots which allow us to study the
# inital dataset based on observables motivated through physics.


import ROOT
from array import array
from math import sqrt
import argparse
from numpy import log as ln

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument('input', help='histogram to be processed')
parser.add_argument('--fit', help='fit to be processed')
parser.add_argument('--output_path', help='path to output folder')
parser.add_argument('--ddb_region', help='path to output folder')
args = parser.parse_args()

labels = {
    # "njets": "Number of leading jets",
    "jpt": "p_{T} (GeV)",
    "jeta": "#eta",
    "jphi": "#phi",
    "jm": "mass (GeV)",
    "jmsoftdrop": "m_{SD} (GeV)",
    "jbtag": "b-tag",
    "jrho": "rho",
}

pt_i = {
        '350': '0',
        '400': '1',
        '450': '2',
        '500': '3',
        '550': '4',
        '600': '5',
        '675': '6',
        '800': '7',
        '1200': '8'
        }

pt_bins = ['350', '400', '450', '500', '550', '600', '675', '800', '1200']

_tag = args.ddb_region

colours = {
    "Hbb": ROOT.kPink+7, #ROOT.TColor.GetColor("#BF2229"),
    "W": ROOT.kGreen+3, #ROOT.TColor.GetColor("#00A88F"),
    "Z": ROOT.kRed, #ROOT.TColor.GetColor(248, 206, 104),
    "tt": ROOT.kMagenta+3, #ROOT.TColor.GetColor(222, 90, 106),
    "QCD": ROOT.kGray+2 ,#ROOT.TColor.GetColor(250, 202, 255),
}

def getHistogram(tfile, name, variable, pt_bin1, tag=_tag):
    name = "{}_{}_pt{}".format(name, tag, pt_i[pt_bin1])
    h = tfile.Get(name)
    if not h:
        raise Exception("Failed to load histogram {}.".format(name))
    return h

def getFitHistogram(tfile, name, pt_bin1):
    histdirname = "shapes_fit_s" + "/ptbin{0}{1}".format(pt_i[pt_bin1], args.ddb_region) + "/" + name
    h = tfile.Get(histdirname)
    if not h:
        raise Exception("Failed to load histogram {}.".format(histdirname))
    h.Scale(7.0)
    return h


# Main function of the plotting step
#
# The major part of the code below is dedicated to define a nice-looking layout.
# The interesting part is the combination of the histograms to the QCD estimation.
# There, we take the data histogram from the control region and subtract all known
# processes defined in simulation and define the remaining part as QCD. Then,
# this shape is extrapolated into the signal region with a scale factor.
def main(variable, pt_bin1, pt_bin2):
    tfile = ROOT.TFile(args.input, "READ")
    tfile_fit = ROOT.TFile(args.fit, "READ")

    # Styles
    ROOT.gStyle.SetOptStat(0)

    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(1000)
    ROOT.gStyle.SetCanvasDefW(600)
    ROOT.gStyle.SetCanvasDefX(0)
    ROOT.gStyle.SetCanvasDefY(0)

    # ROOT.gStyle.SetPadTopMargin(0.02)
    # ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    ROOT.gStyle.SetEndErrorSize(2)
    ROOT.gStyle.SetMarkerStyle(20)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleXOffset(1.00)
    ROOT.gStyle.SetTitleYOffset(1.40)

    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")

    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetPaperSize(20.0, 40.0)
    ROOT.gStyle.SetHatchesLineWidth(5)
    ROOT.gStyle.SetHatchesSpacing(0.05)

    ROOT.gStyle.SetLegendFont(42)

    ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "Y")

    ##
    # Get histogram for each sample
    ##

    Hbb = getFitHistogram(tfile_fit, "ggF", pt_bin1)
    for name in [
        "ZH",
        "WH",
        "ttH",
        "VBF",
    ]:
        Hbb.Add(getFitHistogram(tfile_fit, name, pt_bin1))

    QCD = getFitHistogram(tfile_fit, "qcd", pt_bin1)

    W = getFitHistogram(tfile_fit, "Wjets", pt_bin1)

    Z = getFitHistogram(tfile_fit, "Zjets", pt_bin1)
    for name in [
        "VV",
    ]:
        Z.Add(getFitHistogram(tfile_fit, name, pt_bin1))

    tt = getFitHistogram(tfile_fit, "ttbar", pt_bin1)
    for name in [
        "singlet",
    ]:
        tt.Add(getFitHistogram(tfile_fit, name, pt_bin1))

    data = getHistogram(tfile, "data", variable, pt_bin1)
    if (_tag == "pass"):
      for i in range(11, 15):
        data.SetBinContent(i, 0)
        data.SetBinError(i,0)

    total_bkg = getFitHistogram(tfile_fit, "total_background", pt_bin1)

    ##
    # Compute signal to background ratio
    ##

    N_sig = Hbb.Integral(11,15)
    N_bkg = total_bkg.Integral(11,15)

    print("N_sig: ", N_sig)
    print("N_bkg: ", N_bkg)
    print("N_sig/sqrt(N_bkg): ", N_sig/sqrt(N_bkg))
    print("sqrt(2*(N_sig+N_bkg)*ln(1+N_sig/N_bkg)-2*N_sig): ", sqrt(2*(N_sig+N_bkg)*ln(1+N_sig/N_bkg)-2*N_sig))

    ##
    # Format histograms
    ##

    for x, l in [(Hbb, "Hbb")]:
        x.SetLineWidth(0)
        x.SetFillColor(colours[l])

    for x, l in [(QCD, "QCD"), (W, "W"), (Z, "Z"), (tt, "tt")]:
        x.SetLineWidth(2)
        x.SetLineColor(colours[l])
        if l == "Z":
            x.SetLineStyle(3)
        else:
            x.SetLineStyle(2)

    data.SetMarkerStyle(20)
    data.SetLineColor(ROOT.kBlack)

    total_bkg.SetLineWidth(2)
    total_bkg.SetLineColor(ROOT.kBlue)

    ##
    # Create a stack out of total_bkg background and signal
    ##

    stack_bkg_signal = ROOT.THStack("", "")
    for x in [total_bkg, Hbb]:
       stack_bkg_signal.Add(x)  

    N_data = data.GetEntries()
    
    ratio = data.Clone("ratio")
    for bin in range(data.GetNbinsX()+1):
        data_v = data.GetBinContent(bin)
        bkg_v = total_bkg.GetBinContent(bin)
        sqrt_v = sqrt(N_data)
        rat_v = (data_v - bkg_v) / sqrt_v
        rat_err = data.GetBinError(bin) / sqrt_v
        print("{0} {1} {2} {3} {4} {5}".format(bin, data_v, bkg_v, sqrt_v, rat_v, rat_err))
        ratio.SetBinContent(bin, rat_v)
        ratio.SetBinError(bin, rat_err)
        if (_tag == "fail"):
           ratio.SetBinError(bin, 1.0)
        if (_tag == "pass"):
            if bin >= 11 and bin < 15:
                ratio.SetBinContent(bin, 0)

    # ratio.Add(total_bkg, -1)
    # ratio.Scale(1/N_data)

    ##
    # Draw histograms, stack and data
    ##

    c = ROOT.TCanvas("", "", 600, 700)

    pad1 = ROOT.TPad ('hist', '', 0., 0.3, 1.0, 1.0)
    pad1.SetFillColor(0)
    pad1.SetFillStyle(0)
    pad1.SetBottomMargin(0.03)
    pad1.Draw()

    c.cd()
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1.0, 0.3)
    pad2.SetFillColor(0)
    pad2.SetFillStyle(0)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()

    pad1.cd()

    stack_bkg_signal.Draw("hist")

    name = variable
    if name in labels:
        title = labels[name]
    else:
        title = name
    ratio.GetXaxis().SetTitle(title)
    ratio.GetYaxis().SetTitle("#frac{Data - Bkg}{#sigma_{Data}}")
    ratio.GetXaxis().SetLabelSize(0.1)
    ratio.GetYaxis().SetLabelSize(0.1)
    ratio.GetXaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetTitleOffset(0.75)
    stack_bkg_signal.GetXaxis().SetLabelSize(0)
    stack_bkg_signal.GetYaxis().SetTitle("Events / 7 GeV")
    stack_bkg_signal.SetMaximum(max(stack_bkg_signal.GetMaximum(), data.GetMaximum()) * 1.4)
    stack_bkg_signal.SetMinimum(1.0)

    for x in [QCD, W, Z, tt]:
        x.Draw("HIST SAME")
    data.Draw("E1P SAME")

    # Add legend
    legend = ROOT.TLegend(0.65, 0.55, 0.90, 0.86)
    legend.SetNColumns(1)
    legend.AddEntry(W, "W", "l")
    legend.AddEntry(Z, "Z", "l")
    legend.AddEntry(tt, "t#bar{t}", "l")
    legend.AddEntry(QCD, "Multijet", "l")
    legend.AddEntry(total_bkg, "Total background", "f")
    legend.AddEntry(Hbb, "H(b#bar{b})", "f")
    legend.AddEntry(data, "Data", "lep")
    legend.SetBorderSize(0)
    legend.Draw()

    # Add title
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatex(0.6, 0.935, "2017 (13 TeV)")
    latex.DrawLatex(0.16, 0.935, "#bf{CMS}")
    latex.DrawLatex(0.25, 0.80, pt_bin1 + " < p_{T} < " + pt_bin2 + " GeV")
    if _tag == "fail":
        latex.DrawLatex(0.25, 0.75, "Failing region")
    elif _tag == "pass":
        latex.DrawLatex(0.25, 0.75, "Passing region")

    pad2.cd()
    
    ratio.Draw("E1P")
    Hbb.Scale(sqrt(1/N_data))
    Hbb.Draw("HIST SAME")
    ratio.SetMaximum(max(ratio.GetMaximum(), Hbb.GetMaximum()) * 1.1)
    if (_tag == "fail"):
       ratio.SetMaximum(1.1)

    # Save
    c.SaveAs(args.output_path + "/{}_2017_ptbins_pt{}to{}_{}.png".format(variable,pt_bin1,pt_bin2,_tag))


# Loop over all variable names and make a plot for each
if __name__ == "__main__":
    for i in range(len(pt_bins)-1):
      main("jmsoftdrop", pt_bins[i], pt_bins[i+1])
    # for variable in labels.keys():
    #     main(variable)

