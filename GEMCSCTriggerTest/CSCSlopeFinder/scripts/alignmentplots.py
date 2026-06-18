import ROOT
import os
# fname = "/eos/user/p/pflanaga/singlemu2_50_eta_restricted/output.root"

# f = ROOT.TFile(fname)

# #event = f.Get("GEMCSCBendingAngleTester/AllLCTs") #ToC needs All LCTs
# #event = f.Get("GEMCSCBendingAngleTester/HighestEMTFpT")
# event = f.Get("GEMCSCBendingAngleTester/OnlyRecos") #All LCTs but only ones matched to a reco track (should be faster for eff study)

# files = ["/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_16may2025noalignment/260516_192445/"+folder+"/"+file for folder in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_16may2025noalignment/260516_192445/") for file in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_16may2025noalignment/260516_192445/"+folder) if file.find('log')==-1]
# files = ["/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_18may2025alignment/260519_001457/"+folder+"/"+file for folder in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_18may2025alignment/260519_001457/") for file in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_18may2025alignment/260519_001457/"+folder) if file.find('log')==-1]
files = ["/eos/user/p/pflanaga/singlemu_rerun16june2026/"+folder+"/output.root" for folder in os.listdir("/eos/user/p/pflanaga/singlemu_rerun16june2026/")]

ROOT.gROOT.SetBatch(1)

ROOT.EnableImplicitMT(8)
rdf = ROOT.RDataFrame("GEMCSCBendingAngleTester/OnlyRecos",files)
ROOT.RDF.Experimental.AddProgressBar(rdf)


if not os.path.exists("plots/"):
  os.makedirs("plots/")


plotdir = "plots/alignment_both_unaligned_simulation/"
if not os.path.exists(plotdir):
  os.makedirs(plotdir)

mode_cut = f"((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
base_cut = "has_emtf_track_match"# & (LCT_emu_quality > 3) & (LCT_emu_gemLayer == 1)"
cut = f"{base_cut}  & {mode_cut}"
cut = cut + " & (l1muon_match_pt>30) & (l1muon_match_pt<80)"
# cut = cut + " & (l1muon_match_charge==-1)"
# cut = cut + " & (LCT_CSC_ME1b == 1)"
cut = cut + " & (LCT_emu_eighthstrip < 512)"
# cut = cut + " & (LCT_CSC_endcap==2)"
rdf_tmp = rdf.Filter(cut)
rdf_tmp = rdf_tmp.Define("chargedpt","l1muon_match_pt*l1muon_match_charge")

hist2dendcap1 = ROOT.TH2D("alignmentendcap1", "alignment endcap1;chamber;bending angle",36,0,36,61,-30.5,30.5)
ptchargeendcap1 = ROOT.TH2D("ptchargeendcap1","muon pt*charge in endcap1;chamber;pt*charge",36,0,36,100,-100,100)
hist2dendcap1.SetDirectory(0)
hist2d_modelendcap1 = ROOT.RDF.TH2DModel(hist2dendcap1)
ptcharge_modelendcap1 = ROOT.RDF.TH2DModel(ptchargeendcap1)
etahistendcap1 = ROOT.TH2D("etaendcap1","eta even chambers endcap1;bendingangle;eta",61,-30.5,30.5,20,1.5,2.2)
etahist_modelendcap1 = ROOT.RDF.TH2DModel(etahistendcap1)
etahistoddendcap1 = ROOT.TH2D("etaoddendcap1","eta odd chambers endcap1;bendingangle;eta",61,-30.5,30.5,20,1.5,2.2)
etahist_modeloddendcap1 = ROOT.RDF.TH2DModel(etahistoddendcap1)

alignmentendcap1 = rdf_tmp.Filter("(LCT_CSC_endcap==1)").Define("LCT_emu_charged_slope","LCT_emu_slope*(LCT_emu_bend ? 1 : -1)").Histo2D(hist2d_modelendcap1,"LCT_CSC_chamber","LCT_emu_charged_slope")
ptchargeplotendcap1 = rdf_tmp.Filter("(LCT_CSC_endcap==1)").Histo2D(ptcharge_modelendcap1,"LCT_CSC_chamber","chargedpt")
etaplotevenendcap1 = rdf_tmp.Filter("(LCT_CSC_endcap==1)").Define("LCT_emu_charged_slope","LCT_emu_slope*(LCT_emu_bend ? 1 : -1)").Define("abseta","abs(l1muon_match_eta)").Filter("LCT_CSC_chamber%2==0").Histo2D(etahist_modelendcap1,"LCT_emu_charged_slope","abseta")
etaplotoddendcap1 = rdf_tmp.Filter("(LCT_CSC_endcap==1)").Define("LCT_emu_charged_slope","LCT_emu_slope*(LCT_emu_bend ? 1 : -1)").Define("abseta","abs(l1muon_match_eta)").Filter("LCT_CSC_chamber%2==1").Histo2D(etahist_modeloddendcap1,"LCT_emu_charged_slope","abseta")
ptplotevenendcap1 = rdf_tmp.Filter("(LCT_CSC_endcap==1)").Filter("LCT_CSC_chamber%2==0").Histo1D(("muonptevenendcap1","muonpt even chambers endcap1;pt;count",50,0,100),"l1muon_match_pt")
ptplotoddendcap1 = rdf_tmp.Filter("(LCT_CSC_endcap==1)").Filter("LCT_CSC_chamber%2==1").Histo1D(("muonptoddendcap1","muonpt odd chambers endcap1;pt;count",50,0,100),"l1muon_match_pt")

hist2dendcap2 = ROOT.TH2D("alignmentendcap2", "alignment endcap2;chamber;bending angle",36,0,36,61,-30.5,30.5)
ptchargeendcap2 = ROOT.TH2D("ptchargeendcap2","muon pt*charge in endcap1;chamber;pt*charge",36,0,36,100,-100,100)

hist2dendcap2.SetDirectory(0)
hist2d_modelendcap2 = ROOT.RDF.TH2DModel(hist2dendcap2)
ptcharge_modelendcap2 = ROOT.RDF.TH2DModel(ptchargeendcap2)

etahistendcap2 = ROOT.TH2D("etaendcap2","eta even chambers endcap2;bendingangle;eta",61,-30.5,30.5,20,1.5,2.2)
etahist_modelendcap2 = ROOT.RDF.TH2DModel(etahistendcap2)
etahistoddendcap2 = ROOT.TH2D("etaoddendcap2","eta odd chambers endcap2;bendingangle;eta",61,-30.5,30.5,20,1.5,2.2)
etahist_modeloddendcap2 = ROOT.RDF.TH2DModel(etahistoddendcap2)
alignmentendcap2 = rdf_tmp.Filter("(LCT_CSC_endcap==2)").Define("LCT_emu_charged_slope","LCT_emu_slope*(LCT_emu_bend ? 1 : -1)").Histo2D(hist2d_modelendcap2,"LCT_CSC_chamber","LCT_emu_charged_slope")
ptchargeplotendcap2 = rdf_tmp.Filter("(LCT_CSC_endcap==2)").Histo2D(ptcharge_modelendcap1,"LCT_CSC_chamber","chargedpt")


etaplotevenendcap2 = rdf_tmp.Filter("(LCT_CSC_endcap==2)").Define("LCT_emu_charged_slope","LCT_emu_slope*(LCT_emu_bend ? 1 : -1)").Define("abseta","abs(l1muon_match_eta)").Filter("LCT_CSC_chamber%2==0").Histo2D(etahist_modelendcap2,"LCT_emu_charged_slope","abseta")
etaplotoddendcap2 = rdf_tmp.Filter("(LCT_CSC_endcap==2)").Define("LCT_emu_charged_slope","LCT_emu_slope*(LCT_emu_bend ? 1 : -1)").Define("abseta","abs(l1muon_match_eta)").Filter("LCT_CSC_chamber%2==1").Histo2D(etahist_modeloddendcap2,"LCT_emu_charged_slope","abseta")
ptplotevenendcap2 = rdf_tmp.Filter("(LCT_CSC_endcap==2)").Filter("LCT_CSC_chamber%2==0").Histo1D(("muonptevenendcap2","muonpt even chambers endcap2;pt;count",50,0,100),"l1muon_match_pt")
ptplotoddendcap2 = rdf_tmp.Filter("(LCT_CSC_endcap==2)").Filter("LCT_CSC_chamber%2==1").Histo1D(("muonptoddendcap2","muonpt odd chambers endcap2;pt;count",50,0,100),"l1muon_match_pt")

c = ROOT.TCanvas("c")
alignmentendcap1 = alignmentendcap1.GetPtr()
alignmentendcap1.Draw("colz")
c.SaveAs(plotdir+"bendingangevschamberendcap1.png","png")
etaplotevenendcap1 = etaplotevenendcap1.GetPtr()
etaplotevenendcap1.Draw("colz")
c.SaveAs(plotdir+"etaevenchambersendcap1.png","png")
etaplotoddendcap1 = etaplotoddendcap1.GetPtr()
# etaplotoddendcap1.SetTitle("eta odd chambers;bendingangle;eta")
etaplotoddendcap1.Draw("colz")
c.SaveAs(plotdir+"etaoddchambersendcap1.png","png")
ptplotevenendcap1 = ptplotevenendcap1.GetPtr()
ptplotevenendcap1.Draw("colz")
c.SaveAs(plotdir+"ptplotevenendcap1.png","png")
ptplotoddendcap1 = ptplotoddendcap1.GetPtr()
ptplotoddendcap1.Draw("colz")
c.SaveAs(plotdir+"ptplotoddendcap1.png","png")

alignmentendcap2 = alignmentendcap2.GetPtr()
alignmentendcap2.Draw("colz")
c.SaveAs(plotdir+"bendingangevschamberendcap2.png","png")
etaplotevenendcap2 = etaplotevenendcap2.GetPtr()
etaplotevenendcap2.Draw("colz")
c.SaveAs(plotdir+"etaevenchambersendcap2.png","png")
etaplotoddendcap2 = etaplotoddendcap2.GetPtr()
# etaplotoddendcap2.SetTitle("eta odd chambers;bendingangle;eta")
etaplotoddendcap2.Draw("colz")
c.SaveAs(plotdir+"etaoddchambersendcap2.png","png")
ptplotevenendcap2 = ptplotevenendcap2.GetPtr()
ptplotevenendcap2.Draw("colz")
c.SaveAs(plotdir+"ptplotevenendcap2.png","png")
ptplotoddendcap2 = ptplotoddendcap2.GetPtr()
ptplotoddendcap2.Draw("colz")
c.SaveAs(plotdir+"ptplotoddendcap2.png","png")

ptchargeplotendcap1 = ptchargeplotendcap1.GetPtr()
ptchargeplotendcap1.Draw("colz")
c.SaveAs(plotdir+"chargexptendcap1.png","png")

ptchargeplotendcap2 = ptchargeplotendcap2.GetPtr()
ptchargeplotendcap2.Draw("colz")
c.SaveAs(plotdir+"chargexptendcap2.png","png")