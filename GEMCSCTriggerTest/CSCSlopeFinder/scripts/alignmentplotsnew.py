import ROOT
import os
  

# files = ["/eos/user/p/pflanaga/singlemu_rerun16june2026/"+folder+"/output.root"
#          for folder in os.listdir("/eos/user/p/pflanaga/singlemu_rerun16june2026/")]
# files = ["/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_7june2026noalignment/260609_035229/"+folder+"/"+file for folder in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_7june2026noalignment/260609_035229/") for file in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_7june2026noalignment/260609_035229/"+folder) if (file.find('log')==-1 and file.find('sys')==-1)]
files = ["/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_9june2026alignment/260610_072238/"+folder+"/"+file for folder in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_9june2026alignment/260610_072238/") for file in os.listdir("/eos/user/p/pflanaga/GEMCSCTrigger/2026_ZMu/Muon0/2026_ZMu_9june2026alignment/260610_072238/"+folder) if (file.find('log')==-1 and file.find('sys')==-1)]
  
  

ROOT.gROOT.SetBatch(1)
ROOT.EnableImplicitMT(8)
rdf = ROOT.RDataFrame("GEMCSCBendingAngleTester/OnlyRecos", files)
ROOT.RDF.Experimental.AddProgressBar(rdf)
  

if not os.path.exists("plots/"):
    os.makedirs("plots/")
  

plotdir = "plots/alignment_both_aligned_data/"
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
  

mode_cut = "((emtftrack_mode == 11) | (emtftrack_mode == 13) | (emtftrack_mode == 14) | (emtftrack_mode == 15))"
base_cut = "has_emtf_track_match"
# Cut WITHOUT the pt cut (used for chargedpt vs charged_slope plots)
cut_nopt = f"{base_cut}  & {mode_cut}"
cut_nopt += " & (LCT_emu_eighthstrip < 512) & (l1muon_match_pt>10)"
# Cut WITH the pt cut (used for all other plots)
cut = cut_nopt + " & (l1muon_match_pt>30) & (l1muon_match_pt<50)"
  

def define_columns(d):
    d = d.Define("chargedpt", "l1muon_match_pt*l1muon_match_charge")
    d = d.Define("LCT_emu_charged_slope", "LCT_emu_slope*(LCT_emu_bend ? 1 : -1)")
    d = d.Define("abseta", "abs(l1muon_match_eta)")
    return d

rdf_base = define_columns(rdf.Filter(cut))
rdf_base_nopt = define_columns(rdf.Filter(cut_nopt))
  
  

def make_plots(rdf_in, rdf_in_nopt, charge_label):
    """Generate all plots for a given charge selection. charge_label is appended to plot names."""

    c = ROOT.TCanvas("c_" + charge_label)

    # Storage for the new chargedpt vs charged_slope plots (need to keep refs)
    results = {}

    for endcap in [1, 2]:
        rdf_ec = rdf_in.Filter(f"(LCT_CSC_endcap=={endcap})")
        rdf_ec_nopt = rdf_in_nopt.Filter(f"(LCT_CSC_endcap=={endcap})")

        # ------ alignment (chamber vs bending angle) ------
        hist2d_model = ROOT.RDF.TH2DModel(
            f"alignmentendcap{endcap}_{charge_label}",
            f"alignment endcap{endcap} ({charge_label});chamber;bending angle",
            36, 0, 36, 61, -30.5, 30.5)
        alignment = rdf_ec.Histo2D(hist2d_model, "LCT_CSC_chamber", "LCT_emu_charged_slope")

        # ------ ptcharge vs chamber ------
        ptcharge_model = ROOT.RDF.TH2DModel(
            f"ptchargeendcap{endcap}_{charge_label}",
            f"muon pt*charge in endcap{endcap} ({charge_label});chamber;pt*charge",
            36, 0, 36, 100, -100, 100)
        ptchargeplot = rdf_ec.Histo2D(ptcharge_model, "LCT_CSC_chamber", "chargedpt")

        # ------ eta plots (even/odd) ------
        etahist_model_even = ROOT.RDF.TH2DModel(
            f"etaendcap{endcap}_{charge_label}",
            f"eta even chambers endcap{endcap} ({charge_label});bendingangle;eta",
            61, -30.5, 30.5, 20, 1.5, 2.2)
        etahist_model_odd = ROOT.RDF.TH2DModel(
            f"etaoddendcap{endcap}_{charge_label}",
            f"eta odd chambers endcap{endcap} ({charge_label});bendingangle;eta",
            61, -30.5, 30.5, 20, 1.5, 2.2)

        etaplot_even = rdf_ec.Filter("LCT_CSC_chamber%2==0").Histo2D(
            etahist_model_even, "LCT_emu_charged_slope", "abseta")
        etaplot_odd = rdf_ec.Filter("LCT_CSC_chamber%2==1").Histo2D(
            etahist_model_odd, "LCT_emu_charged_slope", "abseta")

        # ------ pt 1D (even/odd) ------
        ptplot_even = rdf_ec.Filter("LCT_CSC_chamber%2==0").Histo1D(
            (f"muonptevenendcap{endcap}_{charge_label}",
             f"muonpt even chambers endcap{endcap} ({charge_label});pt;count", 50, 0, 100),
            "l1muon_match_pt")
        ptplot_odd = rdf_ec.Filter("LCT_CSC_chamber%2==1").Histo1D(
            (f"muonptoddendcap{endcap}_{charge_label}",
             f"muonpt odd chambers endcap{endcap} ({charge_label});pt;count", 50, 0, 100),
            "l1muon_match_pt")

        # ------ NEW: chargedpt vs LCT_emu_charged_slope, split by even/odd ------
        # Uses rdf_ec_nopt (NO pt cut applied)
        chargedpt_vs_slope_model_even = ROOT.RDF.TH2DModel(
            f"chargedpt_vs_slope_even_endcap{endcap}_{charge_label}",
            f"chargedpt vs charged slope, even chambers, endcap{endcap} ({charge_label});LCT_emu_charged_slope;pt*charge",
            61, -30.5, 30.5, 100, -100, 100)
        chargedpt_vs_slope_model_odd = ROOT.RDF.TH2DModel(
            f"chargedpt_vs_slope_odd_endcap{endcap}_{charge_label}",
            f"chargedpt vs charged slope, odd chambers, endcap{endcap} ({charge_label});LCT_emu_charged_slope;pt*charge",
            61, -30.5, 30.5, 100, -100, 100)

        chargedpt_vs_slope_even = rdf_ec_nopt.Filter("LCT_CSC_chamber%2==0").Histo2D(
            chargedpt_vs_slope_model_even, "LCT_emu_charged_slope", "chargedpt")
        chargedpt_vs_slope_odd = rdf_ec_nopt.Filter("LCT_CSC_chamber%2==1").Histo2D(
            chargedpt_vs_slope_model_odd, "LCT_emu_charged_slope", "chargedpt")

        # ---- Draw and save ----
        alignment.GetPtr().Draw("colz")
        c.SaveAs(f"{plotdir}bendingangevschamberendcap{endcap}_{charge_label}.png", "png")

        etaplot_even.GetPtr().Draw("colz")
        c.SaveAs(f"{plotdir}etaevenchambersendcap{endcap}_{charge_label}.png", "png")

        etaplot_odd.GetPtr().Draw("colz")
        c.SaveAs(f"{plotdir}etaoddchambersendcap{endcap}_{charge_label}.png", "png")

        ptplot_even.GetPtr().Draw("colz")
        c.SaveAs(f"{plotdir}ptplotevenendcap{endcap}_{charge_label}.png", "png")

        ptplot_odd.GetPtr().Draw("colz")
        c.SaveAs(f"{plotdir}ptplotoddendcap{endcap}_{charge_label}.png", "png")

        ptchargeplot.GetPtr().Draw("colz")
        c.SaveAs(f"{plotdir}chargexptendcap{endcap}_{charge_label}.png", "png")

        # ---- New chargedpt vs charged_slope plots with correlation ----
        for parity_label, hist in [("even", chargedpt_vs_slope_even),
                                   ("odd", chargedpt_vs_slope_odd)]:
            h = hist.GetPtr()
            corr = h.GetCorrelationFactor()
            h.Draw("colz")

            # Place text in top-left for negative charge, bottom-left for positive
            if charge_label == "negative":
                x_ndc, y_ndc = 0.15, 0.85
            else:
                x_ndc, y_ndc = 0.15, 0.20

            # Print correlation on the plot
            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.04)
            latex.SetTextColor(ROOT.kRed + 1)
            latex.DrawLatex(x_ndc, y_ndc, f"Correlation = {corr:.4f}")

            c.SaveAs(f"{plotdir}chargedpt_vs_chargedslope_{parity_label}_endcap{endcap}_{charge_label}.png", "png")

            results[(endcap, parity_label)] = (h, corr)
    return results
  
  

# --- Run for negative and positive charges ---
rdf_neg = rdf_base.Filter("(l1muon_match_charge==-1)")
rdf_pos = rdf_base.Filter("(l1muon_match_charge==1)")
rdf_neg_nopt = rdf_base_nopt.Filter("(l1muon_match_charge==-1)")
rdf_pos_nopt = rdf_base_nopt.Filter("(l1muon_match_charge==1)")
  

print("Producing plots for negative charge muons...")
results_neg = make_plots(rdf_neg, rdf_neg_nopt, "negative")
  

print("Producing plots for positive charge muons...")
results_pos = make_plots(rdf_pos, rdf_pos_nopt, "positive")
  

# Print summary of correlations
print("\n=== Correlation summary (chargedpt vs LCT_emu_charged_slope) ===")
for label, results in [("negative", results_neg), ("positive", results_pos)]:
    for (endcap, parity), (_, corr) in results.items():
        print(f"  {label}, endcap {endcap}, {parity} chambers : correlation = {corr:.4f}")