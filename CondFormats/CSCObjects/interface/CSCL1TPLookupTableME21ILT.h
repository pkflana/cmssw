#ifndef CondFormats_CSCObjects_CSCL1TPLookupTableME21ILT_h
#define CondFormats_CSCObjects_CSCL1TPLookupTableME21ILT_h

#include "CSCL1TPLookupTableME11ILT.h"

class CSCL1TPLookupTableME21ILT {
public:
  using t_lut = CSCL1TPLookupTableUtils::t_lut;

  DECLARE_CSCL1TP_LUT(GEM_pad_CSC_es_ME21_even);
  DECLARE_CSCL1TP_LUT(GEM_pad_CSC_es_ME21_odd);
  DECLARE_CSCL1TP_LUT(GEM_roll_L1_CSC_min_wg_ME21_even);
  DECLARE_CSCL1TP_LUT(GEM_roll_L1_CSC_max_wg_ME21_even);
  DECLARE_CSCL1TP_LUT(GEM_roll_L1_CSC_min_wg_ME21_odd);
  DECLARE_CSCL1TP_LUT(GEM_roll_L1_CSC_max_wg_ME21_odd);
  DECLARE_CSCL1TP_LUT(GEM_roll_L2_CSC_min_wg_ME21_even);
  DECLARE_CSCL1TP_LUT(GEM_roll_L2_CSC_max_wg_ME21_even);
  DECLARE_CSCL1TP_LUT(GEM_roll_L2_CSC_min_wg_ME21_odd);
  DECLARE_CSCL1TP_LUT(GEM_roll_L2_CSC_max_wg_ME21_odd);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L1_ME21_even);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L1_ME21_odd);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L2_ME21_even);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L2_ME21_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_2to1_L1_ME21_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_2to1_L1_ME21_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_3to1_L1_ME21_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_3to1_L1_ME21_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_corr_L1_ME21_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_corr_L1_ME21_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L1_ME21_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L1_ME21_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L2_ME21_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L2_ME21_odd);

  unsigned es_diff_slope_bit_width() const;

  // GEM alignment corrections in 1/8-strip units, indexed by 16*(chamber-1)+(roll-1)
  void set_GEM_align_corr_es_ME21_positive_endcap(std::vector<int> lut);
  void set_GEM_align_corr_es_ME21_negative_endcap(std::vector<int> lut);
  int GEM_align_corr_es_ME21_positive_endcap(unsigned chamber, unsigned roll) const;
  int GEM_align_corr_es_ME21_negative_endcap(unsigned chamber, unsigned roll) const;

  COND_SERIALIZABLE;

private:
  std::vector<int> GEM_align_corr_es_ME21_positive_endcap_;
  std::vector<int> GEM_align_corr_es_ME21_negative_endcap_;
};

#endif
