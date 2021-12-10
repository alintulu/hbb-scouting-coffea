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

_tag = args.ddb_region

colours = {
    "Hbb": ROOT.kPink+7, #ROOT.TColor.GetColor("#BF2229"),
    "W": ROOT.kGreen+3, #ROOT.TColor.GetColor("#00A88F"),
    "Z": ROOT.kRed, #ROOT.TColor.GetColor(248, 206, 104),
    "tt": ROOT.kMagenta+3, #ROOT.TColor.GetColor(222, 90, 106),
    "QCD": ROOT.kGray+2 ,#ROOT.TColor.GetColor(250, 202, 255),
}

def getHistogram(tfile, name, variable, tag=_tag):
    name = "{}_{}".format(name, tag)
    h = tfile.Get(name)
    if not h:
        raise Exception("Failed to load histogram {}.".format(name))
    return h

def getFitHistogram(tfile, name):
    h = tfile.Get(name)
    if not h:
        raise Exception("Failed to load histogram {}.".format(name))
    return h


# Main function of the plotting step
#
# The major part of the code below is dedicated to define a nice-looking layout.
# The interesting part is the combination of the histograms to the QCD estimation.
# There, we take the data histogram from the control region and subtract all known
# processes defined in simulation and define the remaining part as QCD. Then,
# this shape is extrapolated into the signal region with a scale factor.
def main(variable):
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

    Hbb = getHistogram(tfile, "ggF", variable)
    for name in [
        "ZH",
        "WH",
        "ttH",
        "VBF",
    ]:
        Hbb.Add(getHistogram(tfile, name, variable))

    histdirname = "shapes_fit_s"
    if _tag == "pass":
        histdirname += "/ptbin0pass"
    elif _tag == "fail":
        histdirname += "/ptbin0fail"
    QCD = getFitHistogram(tfile_fit, histdirname + "/qcd")
    QCD.Scale(7.0)

    W = getHistogram(tfile, "Wjets", variable)

    Z = getHistogram(tfile, "Zjets", variable)
    for name in [
        "VV",
    ]:
        Z.Add(getHistogram(tfile, name, variable))

    tt = getHistogram(tfile, "ttbar", variable)
    for name in [
        "singlet",
    ]:
        tt.Add(getHistogram(tfile, name, variable))

    data = getHistogram(tfile, "data", variable)
    if (_tag == "pass"):
      for i in range(11, 14):
        data.SetBinContent(i, 0)

    total_bkg = QCD.Clone("hnew")
    total_bkg.Add(W)
    total_bkg.Add(Z)
    total_bkg.Add(tt)

    ##
    # Compute signal to background ratio
    ##

    N_sig = 0
    N_bkg = 0
    for bin in range(11, 14):
      N_sig += Hbb.GetBinContent(bin)
      N_bkg += total_bkg.GetBinContent(bin)

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

    print("Variable: {}".format(variable))
    stack_sum = stack_bkg_signal.GetStack().Last()
    bin_center = array('f', [])
    mc_data_diff = array('f', [])
    data_number_events = 0
    for bin in range(stack_sum.GetXaxis().GetNbins()):
        print("Bin {} has {} events".format(bin, data.GetBinContent(bin)))
        data_number_events += data.GetBinContent(bin)
    for bin in range(stack_sum.GetXaxis().GetNbins()):
        bin_center.append(stack_sum.GetXaxis().GetBinCenter(bin))
        if bin >= 11 and bin < 14:
            mc_data_diff.append(0)
        else:
            mc_data_diff.append((data.GetBinContent(bin) - stack_sum.GetBinContent(bin)) / sqrt(data_number_events))

    mc_data_diff_graph = ROOT.TGraph(len(bin_center), bin_center, mc_data_diff)

    ##
    # Draw histograms, stack and data
    ##

    c = ROOT.TCanvas("", "", 600, 700)

    pad1 = ROOT.TPad ('hist', '', 0., 0.3, 1.0, 1.0)
    pad1.SetBottomMargin(0.03)
    pad1.Draw()

    c.cd()
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1.0, 0.3)
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
    mc_data_diff_graph.GetXaxis().SetTitle(title)
    mc_data_diff_graph.GetYaxis().SetTitle("#frac{Data - Bkg}{#sigma_{Data}}")
    mc_data_diff_graph.GetXaxis().SetLabelSize(0.1)
    mc_data_diff_graph.GetYaxis().SetLabelSize(0.1)
    mc_data_diff_graph.GetXaxis().SetTitleSize(0.1)
    mc_data_diff_graph.GetYaxis().SetTitleSize(0.1)
    mc_data_diff_graph.GetYaxis().SetTitleOffset(0.5)
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
    latex.DrawLatex(0.25, 0.80, "350 < p_{T} < 1200 GeV")
    if _tag == "fail":
        latex.DrawLatex(0.25, 0.75, "Failing region")
    elif _tag == "pass":
        latex.DrawLatex(0.25, 0.75, "Passing region")

    pad2.cd()

    mc_data_diff_graph.Draw("AP")
    Hbb.Scale(sqrt(1/data_number_events))
    Hbb.Draw("HIST SAME")
    mc_data_diff_graph.Draw("E1P SAME")
    mc_data_diff_graph.GetYaxis().SetTitleOffset(0.75)
    mc_data_diff_graph.SetMaximum(max(max(mc_data_diff), Hbb.GetMaximum()) * 1.1)

    # Save
    c.SaveAs(args.output_path + "/{}_2017_pt350_{}.png".format(variable,_tag))


# Loop over all variable names and make a plot for each
if __name__ == "__main__":
    main("jmsoftdrop")
    # for variable in labels.keys():
    #     main(variable)

