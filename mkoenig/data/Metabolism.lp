\ Mycoplasma genitalium metabolism
\ Generated 2013-07-10 16:44:43 by edu.stanford.covert.cell.sim.util.FbaLPWriter

Maximize
\ mmol/gDCW
 obj: -3.442882e-06 ADix -3.442882e-06 ADPix -3.442882e-06 AHCYSix -3.442882e-06 AMPix -3.442882e-06 cmnm5s2UMPix -3.442882e-06 CMPix -3.442882e-06 DR5Pix -3.442882e-06 FMETix -3.442882e-06 FORix -3.442882e-06 GDPix -3.442882e-06 GmMPix -3.442882e-06 GMPix -3.442882e-06 GNix -3.442882e-06 Hix -3.442882e-06 k2CMPix -3.442882e-06 LIPOYLLYSix -3.442882e-06 m1GMPix -3.442882e-06 m2GMPix -3.442882e-06 m62AMPix -3.442882e-06 m7GMPix -3.442882e-06 NH3ix -3.442882e-06 NMNix -3.442882e-06 PIix -3.442882e-06 PPIix -3.442882e-06 pSERix -3.442882e-06 PSIURIMPix -3.442882e-06 pTHRix -3.442882e-06 pTYRix -3.442882e-06 s4UMPix -3.442882e-06 SNGLYPix -3.442882e-06 THFix -3.442882e-06 THYix -3.442882e-06 UmMPix -3.442882e-06 UMPix -3.442882e-06 URAix + 2.777778e-01 Growth
Subject To
\ Metabolites
\ Compartments
\ - c: Cytosol
\ - e: Extracellular Space
\ - m: Membrane
 A23CMP_c: -1 CpdB3 = 0
 A3MP_c: 1 CpdB3 = 0
 AC_c: 1 AckA + 1 TX_ACUP = 0
 ACAL_c: 1 DeoC + 1 GlyA4 + 1 PhnX -1 TX_ACAL = 0
 ACCOA_c: 1 AceE -1 Pta -1.513147e-03 Growth = 0
 acriflavine_c: -1 TX_NormM_acriflavine = 0
 ACTP_c: -1 AckA + 1 Pta = 0
 AD_c: -1 Apt + 1 DeoD2 -1 DeoD5 + 1 TX_AD -1 YicP -1 ADix = 0
 ADN_c: 1 AHCYSs + 1 DeoD5 + 1 UshA_AMP = 0
 ADP_c: -1 AckA + 2 Adk1 + 1 Adk2 + 1 Adk3 + 1 Adk4 + 1 ArcC -1 AtpA + 1 CholK + 1 CmkA1 + 1 CmkA2 + 1 CmkB3 + 1 CoaE + 1 Dak1 + 1 Dck + 1 Dgk1 + 1 FruK1 + 1 FruK2 + 1 FthC + 1 FTHFS + 1 GlpK + 1 Gmk1 + 1 Gmk2 + 1 NadK + 1 NADR2 -2 NrdEF_ADP + 1 PcbX + 1 PfkA1 -1 Pgk1 + 1 PncB + 1 PYDXK -1 Pyk_ADP + 1 RibF1 + 1 Tdk1 + 1 Tdk2 + 1 ThiL + 1 Tmk + 1 TX_Cbi_Co + 1 TX_Cbi_Ni + 1 TX_GL + 1 TX_GL3P + 1 TX_GLNHPQ + 1 TX_MdlB_doxorubicin + 1 TX_MGTA_CA + 1 TX_MGTA_CU + 1 TX_MGTA_FE2 + 1 TX_MGTA_FE3 + 1 TX_MGTA_MG + 1 TX_MGTA_MN + 1 TX_MGTA_MOBD + 1 TX_MGTA_ZN + 1 TX_OppB_Ala + 1 TX_OppB_Arg + 1 TX_OppB_Asn + 1 TX_OppB_Asp + 1 TX_OppB_Cys + 1 TX_OppB_Gln + 1 TX_OppB_Glu + 1 TX_OppB_Gly + 1 TX_OppB_His + 1 TX_OppB_Ile + 1 TX_OppB_Leu + 1 TX_OppB_Lys + 1 TX_OppB_Met + 1 TX_OppB_Phe + 1 TX_OppB_Pro + 1 TX_OppB_Ser + 1 TX_OppB_Thr + 1 TX_OppB_Trp + 1 TX_OppB_Tyr + 1 TX_OppB_Val + 1 TX_PHNT + 1 TX_PIUP1 + 1 TX_RIBUP + 1 TX_SalX_nisin + 1 TX_SO4 + 1 TX_SPMDUP1 + 1 TX_SPRDUP2 + 1 TX_SunT_nisin -1 ADPix + 2.652520e+01 Growth = 0
 AEPP_c: -1 PhnW + 1 TX_PHNT = 0
 AHCYS_c: -1 AHCYSs + 1 HCYSMT -1 AHCYSix = 0
 AKG_c: 1 AspC1 + 1 AspC2 + 1 AspC3 + 1 TX_AKG = 0
 ALA_c: 2 Apts_Ala + 1 CsdA -1 GlyA2 -3.329213e-01 Growth = 0
 AlaAla_c: -1 Apts_Ala + 1 TX_OppB_Ala = 0
 AMET_c: -1 HCYSMT + 1 MetK -1.018158e-02 Growth = 0
 aminodimethylaniline_c: 1 AzoR = 0
 AMP_c: 1 Aas1 + 1 Aas2 + 1 Aas3 + 1 Aas4 + 1 Aas5 + 1 Aas6 + 1 Aas7 -1 Adk1 -1 Adk2 -1 Adk3 + 1 Apt + 1 BPNT + 1 NadE + 1 PrsA + 1 RelA1 + 1 RelA2 + 1 TMDPK -1 UshA_AMP -1 AMPix = 0
 anthranilate_c: 1 AzoR = 0
 ARG_c: -1 Adi + 2 Apts_Arg -1.835731e-01 Growth = 0
 ArgArg_c: -1 Apts_Arg + 1 TX_OppB_Arg = 0
 ASN_c: 2 Apts_Asn -3.812638e-01 Growth = 0
 AsnAsn_c: -1 Apts_Asn + 1 TX_OppB_Asn = 0
 ASP_c: 2 Apts_Asp + 1 AspC1 -2.581155e-01 Growth = 0
 AspAsp_c: -1 Apts_Asp + 1 TX_OppB_Asp = 0
 ATP_c: -1 Aas1 -1 Aas2 -1 Aas3 -1 Aas4 -1 Aas5 -1 Aas6 -1 Aas7 + 1 AckA -1 Adk1 -1 Adk4 -1 ArcC + 1 AtpA -1 CholK -1 CmkA1 -1 CmkA2 -1 CmkB3 -1 CoaE -1 Dak1 -1 Dck -1 Dgk1 -2 diadenylate_cyclase -1 FruK1 -1 FruK2 -1 FthC -1 FTHFS -1 GlpK -1 Gmk1 -1 Gmk2 -1 LplA -1 MetK -1 NADD1 -1 NADD2 -1 NadE -1 NadK -1 NADR1 -1 NADR2 -1 PcbX -1 PfkA1 + 1 Pgk1 -1 PncB -1 PrsA -1 PYDXK + 1 Pyk_ADP -1 RelA1 -1 RelA2 -1 RibF1 -1 RibF2 -1 Tdk1 -1 Tdk2 -1 ThiL -1 TMDPK -1 Tmk -1 TX_Cbi_Co -1 TX_Cbi_Ni -1 TX_GL -1 TX_GL3P -1 TX_GLNHPQ -1 TX_MdlB_doxorubicin -1 TX_MGTA_CA -1 TX_MGTA_CU -1 TX_MGTA_FE2 -1 TX_MGTA_FE3 -1 TX_MGTA_MG -1 TX_MGTA_MN -1 TX_MGTA_MOBD -1 TX_MGTA_ZN -1 TX_OppB_Ala -1 TX_OppB_Arg -1 TX_OppB_Asn -1 TX_OppB_Asp -1 TX_OppB_Cys -1 TX_OppB_Gln -1 TX_OppB_Glu -1 TX_OppB_Gly -1 TX_OppB_His -1 TX_OppB_Ile -1 TX_OppB_Leu -1 TX_OppB_Lys -1 TX_OppB_Met -1 TX_OppB_Phe -1 TX_OppB_Pro -1 TX_OppB_Ser -1 TX_OppB_Thr -1 TX_OppB_Trp -1 TX_OppB_Tyr -1 TX_OppB_Val -1 TX_PHNT -1 TX_PIUP1 -1 TX_RIBUP -1 TX_SalX_nisin -1 TX_SO4 -1 TX_SPMDUP1 -1 TX_SPRDUP2 -1 TX_SunT_nisin -1 ATPix -2.661817e+01 Growth = 0
 aza6dT_c: 1 TX_aza6dT + 1 UshA_aza6dTMP = 0
 aza6dTMP_c: -1 UshA_aza6dTMP = 0
 aza8dA_c: 1 TX_aza8dA + 1 UshA_aza8dAMP = 0
 aza8dAMP_c: -1 UshA_aza8dAMP = 0
 benzalkonium_c: -1 TX_NormM_benzalkonium = 0
 berberine_c: -1 TX_NormM_berberine = 0
 C23CMP_c: -1 CpdB2 = 0
 C3MP_c: 1 CpdB2 = 0
 c5dU_c: 1 TX_c5dU + 1 UshA_c5dUMP = 0
 c5dUMP_c: -1 UshA_c5dUMP = 0
 CA_c: 1 TX_MGTA_CA -4.603684e-03 Growth = 0
 CAP_c: 1 ArcC -1 ArgF + 1 TX_CAP = 0
 carbamoylmethyl6dA_c: 1 TX_carbamoylmethyl6dA + 1 UshA_carbamoylmethyl6dAMP = 0
 carbamoylmethyl6dAMP_c: -1 UshA_carbamoylmethyl6dAMP = 0
 CDP_c: 1 CmkA2 -2 NrdEF_CDP -1 Pyk_CDP = 0
 chloroacetaldehyde_c: 1 vinyl_chloride_activation_2 = 0
 chloroethylene_oxide_c: 1 vinyl_chloride_activation_1 -1 vinyl_chloride_activation_2 = 0
 CHOL_c: -1 CholK + 1 GlpQ1 + 1 MG385r1 = 0
 CHP_c: -2 Ohr_CHP -2 OsmC_CHP = 0
 ciprofloxacin_c: -1 TX_NormM_ciprofloxacin = 0
 CITR_c: 1 Adi + 1 ArgF = 0
 CL_c: 1 TX_CL -4.603684e-03 Growth = 0
 cmnm5s2U_c: 1 TX_cmnm5s2U + 1 UshA_cmnm5s2UMP = 0
 cmnm5s2UMP_c: -1 UshA_cmnm5s2UMP -1 cmnm5s2UMPix = 0
 CMP_c: -1 CmkA2 + 1 PgsA1 + 1 PgsA2 + 1 PgsA3 + 1 Udk2 -1 CMPix = 0
 CO_c: 1 TX_Cbi_Co -3.069116e-03 Growth = 0
 CO2_c: 1 AceE -1 ArcC -1 PcbX -1 Ppc + 1 TX_CO2 = 0
 COA_c: -1 AceE -1 AcpS + 1 CoaE + 1 Pta + 1 TX_COA -9.078902e-04 Growth = 0
 crystalViolet_c: -1 TX_NormM_crystalViolet = 0
 CTP_c: -1 CdsA1 -1 CdsA2 -1 CdsA3 -1 CdsA4 -1 CdsA5 -1 CdsA6 -1 CdsA7 + 1 Pyk_CDP -1 CTPix -7.490119e-02 Growth = 0
 CU_c: 1 TX_MGTA_CU -3.069116e-03 Growth = 0
 cyclobutane_dCdC_c: 1 TX_cyclobutane_dCdC + 1 UshA_cyclobutane_dCMPdCMP = 0
 cyclobutane_dCdT_c: 1 TX_cyclobutane_dCdT + 1 UshA_cyclobutane_dCMPdTMP = 0
 cyclobutane_dCMPdCMP_c: -1 UshA_cyclobutane_dCMPdCMP = 0
 cyclobutane_dCMPdTMP_c: -1 UshA_cyclobutane_dCMPdTMP = 0
 cyclobutane_dTdC_c: 1 TX_cyclobutane_dTdC + 1 UshA_cyclobutane_dTMPdCMP = 0
 cyclobutane_dTdT_c: 1 TX_cyclobutane_dTdT + 1 UshA_cyclobutane_dTMPdTMP = 0
 cyclobutane_dTMPdCMP_c: -1 UshA_cyclobutane_dTMPdCMP = 0
 cyclobutane_dTMPdTMP_c: -1 UshA_cyclobutane_dTMPdTMP = 0
 cyclodA_c: 1 TX_cyclodA + 1 UshA_cyclodAMP = 0
 cyclodAMP_c: -1 UshA_cyclodAMP = 0
 cyclodG_c: 1 TX_cyclodG + 1 UshA_cyclodGMP = 0
 cyclodGMP_c: -1 UshA_cyclodGMP = 0
 CYS_c: 2 Apts_Cys -6.101639e-02 Growth = 0
 CysCys_c: -1 Apts_Cys + 1 TX_OppB_Cys = 0
 CYTD_c: -1 Cdd1 + 1 TX_CYTD -1 Udk2 = 0
 c_di_AMP_c: 1 diadenylate_cyclase = 0
 DA_c: -1 Dak1 -1 DeoD2 + 1 TX_DA = 0
 DADP_c: 1 Adk4 + 1 Gmk3 + 1 Gmk4 + 2 NrdEF_ADP + 1 PfkA3 -1 Pgk3 -1 Pyk_DADP = 0
 DALA_c: -1 GlyA1 = 0
 DAMP_c: -1 Adk4 + 1 Dak1 + 2 TatD_dAMP = 0
 dApdAp_c: -1 TatD_dAMP + 1 TX_ComEC_dApdAp = 0
 DATP_c: -1 Gmk3 -1 Gmk4 -1 PfkA3 + 1 Pgk3 + 1 Pyk_DADP -1.800849e-01 Growth = 0
 DC_c: -1 Cdd2 -1 Dck + 1 UshA_dCMP = 0
 dC64dC_c: 1 TX_dC64dC + 1 UshA_dCMP64dCMP = 0
 dC64dC_dewar_c: 1 TX_dC64dC_dewar + 1 UshA_dCMP64dCMP_dewar = 0
 dC64dT_c: 1 TX_dC64dT + 1 UshA_dCMP64dTMP = 0
 dC64dT_dewar_c: 1 TX_dC64dT_dewar + 1 UshA_dCMP64dTMP_dewar = 0
 DCDP_c: 1 CmkA1 + 2 DcdK + 2 NrdEF_CDP -1 Pyk_DCDP = 0
 DCMP_c: -1 CmkA1 -1 DcdK + 1 Dck + 2 TatD_dCMP -1 UshA_dCMP = 0
 dCMP64dCMP_c: -1 UshA_dCMP64dCMP = 0
 dCMP64dCMP_dewar_c: -1 UshA_dCMP64dCMP_dewar = 0
 dCMP64dTMP_c: -1 UshA_dCMP64dTMP = 0
 dCMP64dTMP_dewar_c: -1 UshA_dCMP64dTMP_dewar = 0
 dCMP_GLYC_c: -1 UshA_dCMP_GLYC = 0
 dCpdCp_c: -1 TatD_dCMP + 1 TX_ComEC_dCpdCp = 0
 DCTP_c: -1 DcdK + 1 Pyk_DCDP -8.354063e-02 Growth = 0
 dC_GLYC_c: 1 TX_dC_GLYC + 1 UshA_dCMP_GLYC = 0
 DDCA_c: -1 Aas1 + 1 TX_DDCA = 0
 DDECG3P_c: -1 PlsC1 + 1 PlsX1 = 0
 deoxyG6P_c: -1 YbiV_deoxyG6P = 0
 deoxyGLC_c: 1 YbiV_deoxyG6P = 0
 DG_c: -1 DeoD3 -1 Dgk1 = 0
 DGDP_c: 1 Gmk2 + 1 Gmk4 + 2 NrdEF_GDP -1 Pgk5 -1 Pyk_DGDP = 0
 DGMP_c: 1 Dgk1 -1 Gmk2 -1 Gmk4 + 2 TatD_dGMP = 0
 dGpdGp_c: -1 TatD_dGMP + 1 TX_ComEC_dGpdGp = 0
 DGTP_c: 1 Pgk5 + 1 Pyk_DGDP -8.354063e-02 Growth = 0
 DHA_c: -1 DhaK = 0
 DHdT_c: 1 TX_DHdT + 1 UshA_DHdTMP = 0
 DHdTMP_c: -1 UshA_DHdTMP = 0
 DHdU_c: 1 TX_DHdU + 1 UshA_DHdUMP = 0
 DHdUMP_c: -1 UshA_DHdUMP = 0
 DHF_c: -1 DhfR1 + 1 ThyA = 0
 dhpdU_c: 1 TX_dhpdU + 1 UshA_dhpdUMP = 0
 dhpdUMP_c: -1 UshA_dhpdUMP = 0
 diacylglycerolCys_c: -1 TX_diacylglycerolCys = 0
 dINS_c: -1 DeoD1 + 1 TX_dINS + 1 UshA_dINSMP = 0
 dINSMP_c: -1 UshA_dINSMP = 0
 DMBA_c: 2 Ohr_CHP + 2 OsmC_CHP = 0
 doxorubicin_c: -1 TX_MdlB_doxorubicin -1 TX_NormM_doxorubicin = 0
 DPCOA_c: -1 CoaE = 0
 DPG_c: 1 GapA -1 Pgk1 -1 Pgk2 -1 Pgk3 -1 Pgk4 -1 Pgk5 = 0
 DR1P_c: 1 DeoD1 + 1 DeoD2 + 1 DeoD3 + 1 DeoD8 + 1 Pdp1 + 1 Pdp2 = 0
 DR5P_c: -1 DeoC -1 DR5Pix = 0
 DT_c: -1 Pdp2 -1 Tdk2 + 1 UshA_dTMP = 0
 dT64dC_c: 1 TX_dT64dC + 1 UshA_dTMP64dCMP = 0
 dT64dC_dewar_c: 1 TX_dT64dC_dewar + 1 UshA_dTMP64dCMP_dewar = 0
 dT64dT_c: 1 TX_dT64dT + 1 UshA_dTMP64dTMP = 0
 dT64dT_dewar_c: 1 TX_dT64dT_dewar + 1 UshA_dTMP64dTMP_dewar = 0
 DTDP_c: -1 Pyk_DTDP + 1 Tmk = 0
 DTMP_c: 2 TatD_dTMP + 1 Tdk2 + 1 ThyA -1 Tmk -1 UshA_dTMP = 0
 dTMP64dCMP_c: -1 UshA_dTMP64dCMP = 0
 dTMP64dCMP_dewar_c: -1 UshA_dTMP64dCMP_dewar = 0
 dTMP64dTMP_c: -1 UshA_dTMP64dTMP = 0
 dTMP64dTMP_dewar_c: -1 UshA_dTMP64dTMP_dewar = 0
 dTMP_GLYC_c: -1 UshA_dTMP_GLY = 0
 dTpdTp_c: -1 TatD_dTMP + 1 TX_ComEC_dTpdTp = 0
 DTTP_c: 1 Pyk_DTDP -1.800849e-01 Growth = 0
 dT_GLYC_c: 1 TX_dT_GLYC + 1 UshA_dTMP_GLY = 0
 DU_c: 1 Cdd2 -1 DeoD8 -1 Pdp1 -1 Tdk1 + 1 UshA_dUMP = 0
 DUDP_c: 2 NrdEF_UDP -1 Pyk_DUDP = 0
 DUMP_c: 1 Tdk1 -1 ThyA -1 UshA_dUMP = 0
 DUTP_c: 1 Pyk_DUDP = 0
 dXAN_c: 1 TX_dXAN + 1 UshA_dXANMP = 0
 dXANMP_c: -1 UshA_dXANMP = 0
 e1dA_c: 1 TX_e1dA + 1 UshA_e1dAMP = 0
 e1dAMP_c: -1 UshA_e1dAMP = 0
 e1dG_c: 1 TX_e1dG + 1 UshA_e1dGMP = 0
 e1dGMP_c: -1 UshA_e1dGMP = 0
 e3dA_c: 1 TX_e3dA + 1 UshA_e3dAMP = 0
 e3dAMP_c: -1 UshA_e3dAMP = 0
 e3dC_c: 1 TX_e3dC + 1 UshA_e3dCMP = 0
 e3dCMP_c: -1 UshA_e3dCMP = 0
 e3dG_c: 1 TX_e3dG + 1 UshA_e3dGMP = 0
 e3dGMP_c: -1 UshA_e3dGMP = 0
 e3dT_c: 1 TX_e3dT + 1 UshA_e3dTMP = 0
 e3dTMP_c: -1 UshA_e3dTMP = 0
 e4dC_c: 1 TX_e4dC + 1 UshA_e4dCMP = 0
 e4dCMP_c: -1 UshA_e4dCMP = 0
 E4P_c: 1 TalA -1 Tkt2 = 0
 e5dC_c: 1 TX_e5dC + 1 UshA_e5dCMP = 0
 e5dCMP_c: -1 UshA_e5dCMP = 0
 e6dA_c: 1 TX_e6dA + 1 UshA_e6dAMP = 0
 e6dAMP_c: -1 UshA_e6dAMP = 0
 e7dA_c: 1 TX_e7dA + 1 UshA_e7dAMP = 0
 e7dAMP_c: -1 UshA_e7dAMP = 0
 e7dG_c: 1 TX_e7dG + 1 UshA_e7dGMP = 0
 e7dGMP_c: -1 UshA_e7dGMP = 0
 e8dG_c: 1 TX_e8dG + 1 UshA_e8dGMP = 0
 e8dGMP_c: -1 UshA_e8dGMP = 0
 enoxacin_c: -1 TX_NormM_enoxacin = 0
 epsilon1N2dG_c: 1 TX_epsilon1N2dG + 1 UshA_epsilon1N2dGMP = 0
 epsilon1N2dGMP_c: -1 UshA_epsilon1N2dGMP = 0
 epsilon1N6dA_c: 1 TX_epsilon1N6dA + 1 UshA_epsilon1N6dAMP = 0
 epsilon1N6dAMP_c: -1 UshA_epsilon1N6dAMP = 0
 epsilon3N4dC_c: 1 TX_epsilon3N4dC + 1 UshA_epsilon3N4dCMP = 0
 epsilon3N4dCMP_c: -1 UshA_epsilon3N4dCMP = 0
 epsilonN23dG_c: 1 TX_epsilonN23dG + 1 UshA_epsilonN23dGMP = 0
 epsilonN23dGMP_c: -1 UshA_epsilonN23dGMP = 0
 ethidium_c: -1 TX_NormM_ethidium = 0
 ETN_c: 1 GlpQ2 + 1 MG385r2 = 0
 F1P_c: -1 FruK1 + 1 TX_FRUPTS = 0
 f5dU_c: 1 TX_f5dU + 1 UshA_f5dUMP = 0
 f5dUMP_c: -1 UshA_f5dUMP = 0
 F6P_c: -1 PfkA1 -1 PfkA2 -1 PfkA3 -1 PfkA4 -1 PfkA5 + 1 Pgi + 1 TalA + 1 Tkt2 -1 YbiV_F6P = 0
 FAD_c: -1 Fre1 + 1 Fre2 + 1 GpsA + 1 RibF2 -1.210526e-03 Growth = 0
 FADH2_c: 1 Fre1 -1 Fre2 -1 GpsA = 0
 FALD_c: 1 TX_FALD = 0
 FAPydA_c: 1 TX_FAPydA + 1 UshA_FAPydAMP = 0
 FAPydAMP_c: -1 UshA_FAPydAMP = 0
 FAPydG_c: 1 TX_FAPydG + 1 UshA_FAPydGMP = 0
 FAPydGMP_c: -1 UshA_FAPydGMP = 0
 FDP_c: -1 Fba + 1 FruK1 + 1 PfkA1 + 1 PfkA2 + 1 PfkA3 + 1 PfkA4 + 1 PfkA5 = 0
 FE2_c: 2 Fre2 + 1 TX_MGTA_FE2 -7.016000e-03 Growth = 0
 FE3_c: -2 Fre2 + 1 TX_MGTA_FE3 -6.905516e-03 Growth = 0
 FMET_c: -1 FMETHs -1 FMETix = 0
 FMN_c: -1 Fre5 -1 Fre6 + 2 NrdEF_maintenance + 1 RibF1 -1 RibF2 = 0
 FMNH2_c: 1 Fre5 + 1 Fre6 -2 NrdEF_maintenance = 0
 FOR_c: 1 FMETHs -1 FTHFS + 1 TX_FOR -1 FORix = 0
 FRU_c: 1 YbiV_F6P = 0
 FTHF10_c: 1 FolD2 -1 FthC + 1 FTHFS -1 FTHF10ix -1.210526e-03 Growth = 0
 FTHF5_c: 1 METHFHYD -1.638103e-03 Growth = 0
 FUC1P_c: 1 TX_UHPT2 = 0
 G1P_c: -1 GalU1 -1 Pgm = 0
 G23CMP_c: -1 CpdB4 = 0
 G2P_c: -1 Eno + 1 GpmA = 0
 G3MP_c: 1 CpdB4 = 0
 G3P_c: -1 GpmA + 1 Pgk1 + 1 Pgk2 + 1 Pgk3 + 1 Pgk4 + 1 Pgk5 = 0
 G6P_c: -1 Pgi + 1 Pgm + 1 TX_GLCPTS + 1 TX_UHPT1 -1 YbiV_G6P = 0
 GDP_c: 1 Adk2 + 1 Gmk1 + 1 Gmk3 -2 NrdEF_GDP + 1 PfkA2 -1 Pgk2 -1 Pyk_GDP -1 RelA1 + 1 Udk1 + 1 Udk2 -1 GDPix = 0
 GDPTP_c: 1 RelA2 = 0
 GL_c: 1 Cls1 + 1 Cls2 + 1 Cls3 -1 GlpK + 1 GlpQ4 + 1 MG385r4 + 1 TX_GL + 1 TX_GLUP + 1 YbiV_GL2P + 1 YbiV_SNGLYP = 0
 GL2P_c: -1 YbiV_GL2P = 0
 GL3P_c: 1 GlpK + 1 GlpQ1 + 1 GlpQ2 + 1 GlpQ3 + 1 GlpQ4 + 1 GlpQ5 + 1 glycosyltransferase_PG160_DAG160 + 1 glycosyltransferase_PG161_DAG161 + 1 glycosyltransferase_PG181_DAG181 + 1 GpsA + 1 MG385r1 + 1 MG385r2 + 1 MG385r3 + 1 MG385r4 + 1 MG385r5 -1 PgsA1 -1 PgsA2 -1 PgsA3 -1 PlsX1 -1 PlsX2 -1 PlsX3 -1 PlsX4 -1 PlsX5 -1 PlsX6 -1 PlsX7 + 1 TX_GL3P = 0
 GLC_c: 1 YbiV_G6P = 0
 glchm5dC_c: 1 TX_glchm5dC + 1 UshA_glchm5dCMP = 0
 glchm5dCMP_c: -1 UshA_glchm5dCMP = 0
 glchm5dU_c: 1 TX_glchm5dU + 1 UshA_glchm5dUMP = 0
 glchm5dUMP_c: -1 UshA_glchm5dUMP = 0
 GLN_c: 2 Apts_Gln -1 GLNHs + 1 TX_GLNHPQ -2.419680e-01 Growth = 0
 GlnGln_c: -1 Apts_Gln + 1 TX_OppB_Gln = 0
 GLU_c: 2 Apts_Glu -1 AspC1 -1 AspC2 -1 AspC3 + 1 GLNHs -3.183253e-01 Growth = 0
 gludT_c: 1 TX_gludT + 1 UshA_gludTMP = 0
 gludTMP_c: -1 UshA_gludTMP = 0
 GluGlu_c: -1 Apts_Glu + 1 TX_OppB_Glu = 0
 GLY_c: 2 Apts_Gly + 1 GlyA3 + 1 GlyA4 -2.918973e-01 Growth = 0
 GlyGly_c: -1 Apts_Gly + 1 TX_OppB_Gly = 0
 Gm_c: 1 TX_Gm + 1 UshA_GmMP = 0
 GmMP_c: -1 UshA_GmMP -1 GmMPix = 0
 GMP_c: -1 Gmk1 -1 Gmk3 + 1 Hpt2 -1 UshA_GMP -1 GMPix = 0
 GN_c: 1 DeoD3 -1 DeoD6 -1 Hpt2 + 1 TX_GN -1 GNix = 0
 GPCH_c: -1 GlpQ1 -1 MG385r1 = 0
 GPE_c: -1 GlpQ2 -1 MG385r2 = 0
 GPG_c: -1 GlpQ4 -1 MG385r4 = 0
 GPI_c: -1 GlpQ5 -1 MG385r5 = 0
 GPS_c: -1 GlpQ3 -1 MG385r3 = 0
 GSN_c: 1 DeoD6 + 1 UshA_GMP = 0
 GTP_c: -1 Adk2 -1 PfkA2 + 1 Pgk2 + 1 Pyk_GDP -1 RelA2 -1 Udk1 -1 Udk2 -1 GTPix -8.095696e-02 Growth = 0
 H_c: 1 AcpS + 1 Adi + 1 ArcC + 1 ArgF + 3 AtpA -2 AzoR -1 CdsA1 -1 CdsA2 -1 CdsA3 -1 CdsA4 -1 CdsA5 -1 CdsA6 -1 CdsA7 + 1 CholK + 1 CoaE + 1 CpdB1 + 1 CpdB2 + 1 CpdB3 + 1 CpdB4 -2 CsdA + 1 Dak1 + 1 Dck + 2 DeoD4 + 2 DeoD5 + 2 DeoD6 + 2 DeoD7 + 1 Dgk1 -1 DhfR1 + 1 FMETHs + 2 FolD1 + 1 FolD2 -1 Fre1 + 2 Fre2 -1 Fre3 -1 Fre4 -1 Fre5 -1 Fre6 + 1 FruK1 + 1 FruK2 -3 FthC -1 GalU1 + 1 GapA + 1 GLNHs + 1 GlpQ2 -2 GlyA3 + 1 glycosyltransferase_Ceramide16_GalCer16 + 1 glycosyltransferase_Ceramide16_GlcCer16 + 1 glycosyltransferase_DAG160_GalDAG160 + 1 glycosyltransferase_DAG160_GlcDAG160 + 1 glycosyltransferase_DAG161_GalDAG161 + 1 glycosyltransferase_DAG161_GlcDAG161 + 1 glycosyltransferase_DAG181_GalDAG181 + 1 glycosyltransferase_DAG181_GlcDAG181 + 1 glycosyltransferase_GalCer16_GalGalCer16 + 1 glycosyltransferase_GalCer16_GlcGalCer16 + 1 glycosyltransferase_GalDAG160_GalGalDAG160 + 1 glycosyltransferase_GalDAG160_GlcGalDAG160 + 1 glycosyltransferase_GalDAG161_GalGalDAG161 + 1 glycosyltransferase_GalDAG161_GlcGalDAG161 + 1 glycosyltransferase_GalDAG181_GalGalDAG181 + 1 glycosyltransferase_GalDAG181_GlcGalDAG181 + 1 glycosyltransferase_GalGalDAG160_GalGalGalDAG160 + 1 glycosyltransferase_GalGalDAG160_GlcGalGalDAG160 + 1 glycosyltransferase_GalGalDAG161_GalGalGalDAG161 + 1 glycosyltransferase_GalGalDAG161_GlcGalGalDAG161 + 1 glycosyltransferase_GalGalDAG181_GalGalGalDAG181 + 1 glycosyltransferase_GalGalDAG181_GlcGalGalDAG181 + 1 glycosyltransferase_GalGlcCer16_GalGalGlcCer16 + 1 glycosyltransferase_GalGlcDAG160_GalGalGlcDAG160 + 1 glycosyltransferase_GalGlcDAG160_GlcGalGlcDAG160 + 1 glycosyltransferase_GalGlcDAG161_GalGalGlcDAG161 + 1 glycosyltransferase_GalGlcDAG161_GlcGalGlcDAG161 + 1 glycosyltransferase_GalGlcDAG181_GalGalGlcDAG181 + 1 glycosyltransferase_GalGlcDAG181_GlcGalGlcDAG181 + 1 glycosyltransferase_GlcCer16_GalGlcCer16 + 1 glycosyltransferase_GlcDAG160_GalGlcDAG160 + 1 glycosyltransferase_GlcDAG160_GlcGlcDAG160 + 1 glycosyltransferase_GlcDAG161_GalGlcDAG161 + 1 glycosyltransferase_GlcDAG161_GlcGlcDAG161 + 1 glycosyltransferase_GlcDAG181_GalGlcDAG181 + 1 glycosyltransferase_GlcDAG181_GlcGlcDAG181 + 1 glycosyltransferase_GlcGlcDAG160_GlcGlcGlcDAG160 + 1 glycosyltransferase_GlcGlcDAG161_GlcGlcGlcDAG161 + 1 glycosyltransferase_GlcGlcDAG181_GlcGlcGlcDAG181 + 2 glycosyltransferase_PG160_DAG160 -1 GpsA + 1 HCYSMT + 2 hydroxyl_radical_generation_1 -1 LdhA -1 LplA + 1 Mdh -1 METHFHYD + 1 MetS + 1 MG385r2 -1 NADD1 -1 NADD2 + 3 NadE -1 NADR1 -1 NH3eq -1 NOX + 2 PcbX -2 Pdp3 + 1 PfkA1 + 1 PfkA2 + 1 PfkA3 + 1 PfkA4 + 1 PfkA5 + 2 PgsA1 + 2 PgsA2 + 2 PgsA3 + 2 PhnW + 2 PhnX -2 PlsC4 + 2 PlsC5 + 1 PlsX1 + 1 PlsX2 + 1 PlsX3 + 3 PlsX4 -1 PlsX5 + 1 PlsX6 + 1 PlsX7 -3 PncB + 1 Ppa + 1 Ppc + 1 PrsA + 1 PYDXK -1 Pyk_ADP -1 Pyk_CDP -1 Pyk_DADP -1 Pyk_DCDP -1 Pyk_DGDP -1 Pyk_DTDP -1 Pyk_DUDP -1 Pyk_GDP -1 Pyk_IDP -1 Pyk_TDP -1 Pyk_UDP + 1 RelA1 + 1 RelA2 + 1 RibF1 -1 RibF2 + 1 SNGLYPOr + 1 TatD_dAMP + 1 TatD_dCMP + 1 TatD_dGMP + 1 TatD_dTMP + 1 Tdk1 + 1 Tdk2 + 2 ThyA + 1 TMDPK -1 TrxB + 1 TXTHM + 1 TX_ACUP + 1 TX_AD + 1 TX_AKG + 1 TX_AROP11 + 1 TX_AROP12 + 1 TX_AROP13 + 1 TX_AROP21 + 1 TX_AROP22 + 1 TX_AROP23 + 1 TX_Cbi_Co + 1 TX_Cbi_Ni -2 TX_CL + 1 TX_CYTD + 1 TX_DA + 1 TX_GL + 1 TX_GL3P + 1 TX_GLNHPQ + 1 TX_GN + 1 TX_MdlB_doxorubicin + 1 TX_MGTA_CA + 1 TX_MGTA_CU + 1 TX_MGTA_FE2 + 1 TX_MGTA_FE3 + 1 TX_MGTA_MG + 1 TX_MGTA_MN + 1 TX_MGTA_MOBD + 1 TX_MGTA_ZN -1 TX_NAUP2 + 1 TX_OppB_Ala + 1 TX_OppB_Arg + 1 TX_OppB_Asn + 1 TX_OppB_Asp + 1 TX_OppB_Cys + 1 TX_OppB_Gln + 1 TX_OppB_Glu + 1 TX_OppB_Gly + 1 TX_OppB_His + 1 TX_OppB_Ile + 1 TX_OppB_Leu + 1 TX_OppB_Lys + 1 TX_OppB_Met + 1 TX_OppB_Phe + 1 TX_OppB_Pro + 1 TX_OppB_Ser + 1 TX_OppB_Thr + 1 TX_OppB_Trp + 1 TX_OppB_Tyr + 1 TX_OppB_Val + 1 TX_PHNT + 1 TX_PIUP1 + 1 TX_RIBUP + 1 TX_SalX_nisin + 1 TX_SO4 + 1 TX_SPMDUP1 + 1 TX_SPRDUP2 + 1 TX_SunT_nisin -1 TX_THY + 1 TX_TRK + 1 TX_URA + 1 Udk1 + 1 Udk2 + 2 vinyl_chloride_activation_1 + 2 YbiV_deoxyG6P + 1 YbiV_SBT6P -1 Hix + 2.652531e+01 Growth = 0
 H2O_c: -1 Adi -1 AHCYSs -1 Apts_Ala -1 Apts_Arg -1 Apts_Asn -1 Apts_Asp -1 Apts_Cys -1 Apts_Gln -1 Apts_Glu -1 Apts_Gly -1 Apts_His -1 Apts_Ile -1 Apts_Leu -1 Apts_Lys -1 Apts_Met -1 Apts_Phe -1 Apts_Pro -1 Apts_Ser -1 Apts_Thr -1 Apts_Trp -1 Apts_Tyr -1 Apts_Val + 1 AtpA -1 BPNT -1 Cdd1 -1 Cdd2 -1 CpdB1 -1 CpdB2 -1 CpdB3 -1 CpdB4 + 1 Eno -1 FMETHs -1 FolD2 -1 GLNHs -1 GlpQ1 -1 GlpQ2 -1 GlpQ3 -1 GlpQ4 -1 GlpQ5 + 1 GlyA3 -1 glycosyltransferase_PG160_DAG160 -1 glycosyltransferase_PG161_DAG161 -1 glycosyltransferase_PG181_DAG181 -2 hydroxyl_radical_generation_1 -1 hydroxyl_radical_generation_2 -1 LIPOYLLYSHs -1 METHFHYD -1 MetK -1 MG385r1 -1 MG385r2 -1 MG385r3 -1 MG385r4 -1 MG385r5 -1 MsrA -1 MsrB + 2 NrdEF_ADP + 2 NrdEF_CDP + 2 NrdEF_GDP + 2 NrdEF_maintenance + 2 NrdEF_UDP + 2 Ohr_CHP + 4 Ohr_H2O2 + 2 Ohr_tBHP + 2 OsmC_CHP + 4 OsmC_H2O2 + 2 OsmC_tBHP -1 PcbX -1 PgpA1 -1 PgpA2 -1 PgpA3 -1 PhnW -1 PhnX -1 PncB -1 Ppa -1 Ppc -1 pSERpp -1 pTHRpp -1 pTYRpp -1 TatD_dAMP -1 TatD_dCMP -1 TatD_dGMP -1 TatD_dTMP -1 TX_Cbi_Co -1 TX_Cbi_Ni -1 TX_GL -1 TX_GL3P -1 TX_GLNHPQ + 1 TX_H2O -1 TX_MdlB_doxorubicin -1 TX_MGTA_CA -1 TX_MGTA_CU -1 TX_MGTA_FE2 -1 TX_MGTA_FE3 -1 TX_MGTA_MG -1 TX_MGTA_MN -1 TX_MGTA_MOBD -1 TX_MGTA_ZN -1 TX_OppB_Ala -1 TX_OppB_Arg -1 TX_OppB_Asn -1 TX_OppB_Asp -1 TX_OppB_Cys -1 TX_OppB_Gln -1 TX_OppB_Glu -1 TX_OppB_Gly -1 TX_OppB_His -1 TX_OppB_Ile -1 TX_OppB_Leu -1 TX_OppB_Lys -1 TX_OppB_Met -1 TX_OppB_Phe -1 TX_OppB_Pro -1 TX_OppB_Ser -1 TX_OppB_Thr -1 TX_OppB_Trp -1 TX_OppB_Tyr -1 TX_OppB_Val -1 TX_PHNT -1 TX_PIUP1 -1 TX_RIBUP -1 TX_SalX_nisin -1 TX_SO4 -1 TX_SPMDUP1 -1 TX_SPRDUP2 -1 TX_SunT_nisin -1 UshA_AMP -1 UshA_aza6dTMP -1 UshA_aza8dAMP -1 UshA_c5dUMP -1 UshA_carbamoylmethyl6dAMP -1 UshA_cmnm5s2UMP -2 UshA_cyclobutane_dCMPdCMP -2 UshA_cyclobutane_dCMPdTMP -2 UshA_cyclobutane_dTMPdCMP -2 UshA_cyclobutane_dTMPdTMP -1 UshA_cyclodAMP -1 UshA_cyclodGMP -1 UshA_dCMP -2 UshA_dCMP64dCMP -2 UshA_dCMP64dCMP_dewar -2 UshA_dCMP64dTMP -2 UshA_dCMP64dTMP_dewar -1 UshA_dCMP_GLYC -1 UshA_DHdTMP -1 UshA_DHdUMP -1 UshA_dhpdUMP -1 UshA_dINSMP -1 UshA_dTMP -2 UshA_dTMP64dCMP -2 UshA_dTMP64dCMP_dewar -2 UshA_dTMP64dTMP -2 UshA_dTMP64dTMP_dewar -1 UshA_dTMP_GLY -1 UshA_dUMP -1 UshA_dXANMP -1 UshA_e1dAMP -1 UshA_e1dGMP -1 UshA_e3dAMP -1 UshA_e3dCMP -1 UshA_e3dGMP -1 UshA_e3dTMP -1 UshA_e4dCMP -1 UshA_e5dCMP -1 UshA_e6dAMP -1 UshA_e7dAMP -1 UshA_e7dGMP -1 UshA_e8dGMP -1 UshA_epsilon1N2dGMP -1 UshA_epsilon1N6dAMP -1 UshA_epsilon3N4dCMP -1 UshA_epsilonN23dGMP -1 UshA_f5dUMP -1 UshA_FAPydAMP -1 UshA_FAPydGMP -1 UshA_glchm5dCMP -1 UshA_glchm5dUMP -1 UshA_gludTMP -1 UshA_GmMP -1 UshA_GMP -1 UshA_hm3dCMP -1 UshA_hm5dCMP -1 UshA_hm5dUMP -1 UshA_hm5hydantoin_dRiboseMP -1 UshA_ho5dCMP -1 UshA_ho5dUMP -1 UshA_ho5ho6dTMP -1 UshA_ho5ho6dUMP -1 UshA_ho5hydantoin_dRiboseMP -1 UshA_IMP -1 UshA_k2CMP -1 UshA_m1dAMP -1 UshA_m1dGMP -1 UshA_m1GMP -1 UshA_m2GMP -1 UshA_m3dAMP -1 UshA_m3dCMP -1 UshA_m3dGMP -1 UshA_m3dTMP -1 UshA_m4dCMP -1 UshA_m5CMP -1 UshA_m5dCMP -1 UshA_m62AMP -1 UshA_m6AMP -1 UshA_m6dAMP -1 UshA_m7dAMP -1 UshA_m7dGMP -1 UshA_m7GMP -1 UshA_m8dGMP -1 UshA_n2dAMP -1 UshA_O2edCMP -1 UshA_O2edTMP -1 UshA_O2mdCMP -1 UshA_O2mdTMP -1 UshA_O4edTMP -1 UshA_O4mdTMP -1 UshA_O6edGMP -1 UshA_O6mdGMP -1 UshA_oxo2dAMP -1 UshA_oxo8dAMP -1 UshA_oxo8dGMP -1 UshA_oxoethyl7dGMP -1 UshA_PSIURIMP -1 UshA_putdTMP -1 UshA_s2UMP -1 UshA_s4UMP -1 UshA_UmMP -1 vinyl_chloride_activation_1 -1 YbiV_deoxyG6P -1 YbiV_F6P -1 YbiV_G6P -1 YbiV_GL2P -1 YbiV_MAN6P -1 YbiV_R5P -1 YbiV_SBT6P -1 YbiV_SNGLYP -1 YicP -1 H2Oix -1.507214e+02 Growth = 0
 H2O2_c: 3 HemG + 1 NOX -2 Ohr_H2O2 -2 OsmC_H2O2 + 1 TX_H2O2 = 0
 HCYS_c: 1 AHCYSs -1 HCYSMT -1 MetS = 0
 HDCA_c: -1 Aas4 + 1 TX_HDCA = 0
 HDCEA_c: -1 Aas5 + 1 TX_HDCEA = 0
 HDECEG3P_c: -1 PlsC5 + 1 PlsX5 = 0
 HDECG3P_c: -1 PlsC4 + 1 PlsX4 = 0
 HIS_c: 2 Apts_His -8.916427e-02 Growth = 0
 HisHis_c: -1 Apts_His + 1 TX_OppB_His = 0
 hm3dC_c: 1 TX_hm3dC + 1 UshA_hm3dCMP = 0
 hm3dCMP_c: -1 UshA_hm3dCMP = 0
 hm5dC_c: 1 TX_hm5dC + 1 UshA_hm5dCMP = 0
 hm5dCMP_c: -1 UshA_hm5dCMP = 0
 hm5dU_c: 1 TX_hm5dU + 1 UshA_hm5dUMP = 0
 hm5dUMP_c: -1 UshA_hm5dUMP = 0
 hm5hydantoin_dRibose_c: 1 TX_hm5hydantoin_dRibose + 1 UshA_hm5hydantoin_dRiboseMP = 0
 hm5hydantoin_dRiboseMP_c: -1 UshA_hm5hydantoin_dRiboseMP = 0
 ho5dC_c: 1 TX_ho5dC + 1 UshA_ho5dCMP = 0
 ho5dCMP_c: -1 UshA_ho5dCMP = 0
 ho5dU_c: 1 TX_ho5dU + 1 UshA_ho5dUMP = 0
 ho5dUMP_c: -1 UshA_ho5dUMP = 0
 ho5ho6dT_c: 1 TX_ho5ho6dT + 1 UshA_ho5ho6dTMP = 0
 ho5ho6dTMP_c: -1 UshA_ho5ho6dTMP = 0
 ho5ho6dU_c: 1 TX_ho5ho6dU + 1 UshA_ho5ho6dUMP = 0
 ho5ho6dUMP_c: -1 UshA_ho5ho6dUMP = 0
 ho5hydantoin_dRibose_c: 1 TX_ho5hydantoin_dRibose + 1 UshA_ho5hydantoin_dRiboseMP = 0
 ho5hydantoin_dRiboseMP_c: -1 UshA_ho5hydantoin_dRiboseMP = 0
 HPHPYR_c: -1 AspC2 = 0
 hydride_c: 1 hydroxyl_radical_generation_2 = 0
 hydroxyl_radical_c: 1 hydroxyl_radical_generation_1 + 1 hydroxyl_radical_generation_2 = 0
 HYXN_c: 1 DeoD1 -1 DeoD4 -1 Hpt1 + 1 YicP = 0
 IDP_c: 1 Adk3 + 1 PfkA4 -1 Pgk4 -1 Pyk_IDP = 0
 ILE_c: 2 Apts_Ile -4.498328e-01 Growth = 0
 IleIle_c: -1 Apts_Ile + 1 TX_OppB_Ile = 0
 IMP_c: 1 Hpt1 -1 UshA_IMP = 0
 INOS_c: 1 GlpQ5 + 1 MG385r5 = 0
 INS_c: 1 DeoD4 + 1 UshA_IMP = 0
 ITP_c: -1 Adk3 -1 PfkA4 + 1 Pgk4 + 1 Pyk_IDP = 0
 K_c: 1 TX_TRK -1.728266e-01 Growth = 0
 k2C_c: 1 TX_k2C + 1 UshA_k2CMP = 0
 k2CMP_c: -1 UshA_k2CMP -1 k2CMPix = 0
 LAC_c: 1 LdhA = 0
 LEU_c: 2 Apts_Leu -5.507228e-01 Growth = 0
 LeuLeu_c: -1 Apts_Leu + 1 TX_OppB_Leu = 0
 LIPOATE_c: 1 LIPOYLLYSHs -1 LplA + 1 TX_LIPOATE = 0
 LIPOYLAMP_c: 1 LplA -4.410944e-04 Growth = 0
 LIPOYLLYS_c: -1 LIPOYLLYSHs -1 LIPOYLLYSix = 0
 lomefloxacin_c: -1 TX_NormM_lomefloxacin = 0
 LYS_c: 2 Apts_Lys + 1 LIPOYLLYSHs -5.272188e-01 Growth = 0
 LysLys_c: -1 Apts_Lys + 1 TX_OppB_Lys = 0
 m1dA_c: 1 TX_m1dA + 1 UshA_m1dAMP = 0
 m1dAMP_c: -1 UshA_m1dAMP = 0
 m1dG_c: 1 TX_m1dG + 1 UshA_m1dGMP = 0
 m1dGMP_c: -1 UshA_m1dGMP = 0
 m1G_c: 1 TX_m1G + 1 UshA_m1GMP = 0
 m1GMP_c: -1 UshA_m1GMP -1 m1GMPix = 0
 m2G_c: 1 TX_m2G + 1 UshA_m2GMP = 0
 m2GMP_c: -1 UshA_m2GMP -1 m2GMPix = 0
 m3dA_c: 1 TX_m3dA + 1 UshA_m3dAMP = 0
 m3dAMP_c: -1 UshA_m3dAMP = 0
 m3dC_c: 1 TX_m3dC + 1 UshA_m3dCMP = 0
 m3dCMP_c: -1 UshA_m3dCMP = 0
 m3dG_c: 1 TX_m3dG + 1 UshA_m3dGMP = 0
 m3dGMP_c: -1 UshA_m3dGMP = 0
 m3dT_c: 1 TX_m3dT + 1 UshA_m3dTMP = 0
 m3dTMP_c: -1 UshA_m3dTMP = 0
 m4dC_c: 1 TX_m4dC + 1 UshA_m4dCMP = 0
 m4dCMP_c: -1 UshA_m4dCMP = 0
 m5C_c: 1 TX_m5C + 1 UshA_m5CMP = 0
 m5CMP_c: -1 UshA_m5CMP = 0
 m5dC_c: 1 TX_m5dC + 1 UshA_m5dCMP = 0
 m5dCMP_c: -1 UshA_m5dCMP = 0
 m62A_c: 1 TX_m62A + 1 UshA_m62AMP = 0
 m62AMP_c: -1 UshA_m62AMP -1 m62AMPix = 0
 m6A_c: 1 TX_m6A + 1 UshA_m6AMP = 0
 m6AMP_c: -1 UshA_m6AMP = 0
 m6dA_c: 1 TX_m6dA + 1 UshA_m6dAMP = 0
 m6dAMP_c: -1 UshA_m6dAMP = 0
 m7dA_c: 1 TX_m7dA + 1 UshA_m7dAMP = 0
 m7dAMP_c: -1 UshA_m7dAMP = 0
 m7dG_c: 1 TX_m7dG + 1 UshA_m7dGMP = 0
 m7dGMP_c: -1 UshA_m7dGMP = 0
 m7G_c: 1 TX_m7G + 1 UshA_m7GMP = 0
 m7GMP_c: -1 UshA_m7GMP -1 m7GMPix = 0
 m8dG_c: 1 TX_m8dG + 1 UshA_m8dGMP = 0
 m8dGMP_c: -1 UshA_m8dGMP = 0
 MAL_c: -1 Mdh = 0
 MAN_c: 1 YbiV_MAN6P = 0
 MAN1P_c: 1 CpsG = 0
 MAN6P_c: -1 CpsG + 1 TX_UHPT3 -1 YbiV_MAN6P = 0
 MET_c: 2 Apts_Met + 1 FMETHs + 1 HCYSMT -1 MetK + 1 MetS -1 MsrA -1 MsrB -1 METix -9.236000e-02 Growth = 0
 METHF_c: 1 FolD1 -1 FolD2 -1 METHFHYD = 0
 methylRed_c: -1 AzoR = 0
 MetMet_c: -1 Apts_Met + 1 TX_OppB_Met = 0
 METSOXRL_c: 1 MsrB = 0
 METSOXSL_c: 1 MsrA = 0
 METTHF_c: -1 FolD1 + 1 FthC + 1 GlyA3 -1 MetF -1 ThyA -1.210526e-03 Growth = 0
 MG_c: 1 TX_MGTA_MG -1.015642e-02 Growth = 0
 MN_c: 1 TX_MGTA_MN -3.522324e-03 Growth = 0
 MO_c: 1 MOrx -3.069116e-03 Growth = 0
 MOBD_c: -1 MOrx + 1 TX_MGTA_MOBD = 0
 MTHF_c: 1 MetF -1 MetS -1.210526e-03 Growth = 0
 n2dA_c: 1 TX_n2dA + 1 UshA_n2dAMP = 0
 n2dAMP_c: -1 UshA_n2dAMP = 0
 NA_c: 1 TX_NAUP2 + 1 TX_NormM_acriflavine + 1 TX_NormM_benzalkonium + 1 TX_NormM_berberine + 1 TX_NormM_ciprofloxacin + 1 TX_NormM_crystalViolet + 1 TX_NormM_doxorubicin + 1 TX_NormM_enoxacin + 1 TX_NormM_ethidium + 1 TX_NormM_lomefloxacin + 1 TX_NormM_norfloxacin + 1 TX_NormM_novobiocin + 1 TX_NormM_proflavine + 1 TX_NormM_rhodamine6G + 1 TX_NormM_tetraphenylphosphonium -3.941782e-03 Growth = 0
 NAAD_c: 1 NADD1 -1 NadE = 0
 NAC_c: -1 PncB + 1 TX_NAC = 0
 NAD_c: -1 AceE + 2 AzoR + 1 Fre1 + 1 Fre4 + 1 Fre5 -1 GapA + 1 LdhA -1 Mdh + 1 MetF + 1 NADD2 + 1 NadE -1 NadK + 1 NADR1 + 1 NOX -9.848224e-03 Growth = 0
 NADH_c: 1 AceE -2 AzoR -1 Fre1 -1 Fre4 -1 Fre5 + 1 GapA -1 LdhA + 1 Mdh -1 MetF -1 NOX -2.421041e-04 Growth = 0
 NADP_c: 1 DhfR1 -1 FolD1 + 1 Fre3 + 1 Fre6 + 1 NadK -1 SNGLYPOr + 1 TrxB -6.052609e-04 Growth = 0
 NADPH_c: -1 DhfR1 + 1 FolD1 -1 Fre3 -1 Fre6 + 1 SNGLYPOr -1 TrxB -1.815779e-03 Growth = 0
 NAMN_c: -1 NADD1 + 1 PncB = 0
 NH3_c: 1 Adi -1 ArcC + 1 Cdd1 + 1 Cdd2 + 1 GLNHs -1 NadE -1 NADR2 -1 NH3eq + 1 PhnW + 1 YicP -1 NH3ix = 0
 NH4_c: 1 NH3eq + 1 TX_NH4 -1.150921e-02 Growth = 0
 NI_c: 1 TX_Cbi_Ni -7.596073e-04 Growth = 0
 nisin_c: -1 TX_SalX_nisin -1 TX_SunT_nisin = 0
 NMN_c: -1 NADD2 -1 NADR1 + 1 NADR2 -1 NMNix = 0
 norfloxacin_c: -1 TX_NormM_norfloxacin = 0
 novobiocin_c: -1 TX_NormM_novobiocin = 0
 NR_c: -1 NADR2 = 0
 O2_c: -3 HemG + 2 MOrx -1 NOX -1 NrdEF_ADP -1 NrdEF_CDP -1 NrdEF_GDP -1 NrdEF_UDP + 1 TX_O2 = 0
 O2edC_c: 1 TX_O2edC + 1 UshA_O2edCMP = 0
 O2edCMP_c: -1 UshA_O2edCMP = 0
 O2edT_c: 1 TX_O2edT + 1 UshA_O2edTMP = 0
 O2edTMP_c: -1 UshA_O2edTMP = 0
 O2mdC_c: 1 TX_O2mdC + 1 UshA_O2mdCMP = 0
 O2mdCMP_c: -1 UshA_O2mdCMP = 0
 O2mdT_c: 1 TX_O2mdT + 1 UshA_O2mdTMP = 0
 O2mdTMP_c: -1 UshA_O2mdTMP = 0
 O4edT_c: 1 TX_O4edT + 1 UshA_O4edTMP = 0
 O4edTMP_c: -1 UshA_O4edTMP = 0
 O4mdT_c: 1 TX_O4mdT + 1 UshA_O4mdTMP = 0
 O4mdTMP_c: -1 UshA_O4mdTMP = 0
 O6edG_c: 1 TX_O6edG + 1 UshA_O6edGMP = 0
 O6edGMP_c: -1 UshA_O6edGMP = 0
 O6mdG_c: 1 TX_O6mdG + 1 UshA_O6mdGMP = 0
 O6mdGMP_c: -1 UshA_O6mdGMP = 0
 OA_c: -1 AspC1 + 1 Mdh + 1 PcbX + 1 Ppc + 1 TX_OA = 0
 OCDCA_c: -1 Aas6 + 1 TX_OCDCA = 0
 OCDCEA_c: -1 Aas7 + 1 TX_OCDCEA = 0
 ODEC11GP_c: -1 PlsC7 + 1 PlsX7 = 0
 ODECG3P_c: -1 PlsC6 + 1 PlsX6 = 0
 OH_c: 1 hydroxyl_radical_generation_1 = 0
 ORN_c: -1 ArgF = 0
 oxo2dA_c: 1 TX_oxo2dA + 1 UshA_oxo2dAMP = 0
 oxo2dAMP_c: -1 UshA_oxo2dAMP = 0
 oxo8dA_c: 1 TX_oxo8dA + 1 UshA_oxo8dAMP = 0
 oxo8dAMP_c: -1 UshA_oxo8dAMP = 0
 oxo8dG_c: 1 TX_oxo8dG + 1 UshA_oxo8dGMP = 0
 oxo8dGMP_c: -1 UshA_oxo8dGMP = 0
 oxoethyl7dG_c: 1 TX_oxoethyl7dG + 1 UshA_oxoethyl7dGMP = 0
 oxoethyl7dGMP_c: -1 UshA_oxoethyl7dGMP = 0
 PAP_c: 1 AcpS -1 BPNT = 0
 pCHOL_c: 1 CholK = 0
 PDLA5P_c: 1 GlyA1 + 1 GlyA2 = 0
 PEP_c: -1 DhaK + 1 Eno -1 Ppc -1 Pyk_ADP -1 Pyk_CDP -1 Pyk_DADP -1 Pyk_DCDP -1 Pyk_DGDP -1 Pyk_DTDP -1 Pyk_DUDP -1 Pyk_GDP -1 Pyk_IDP -1 Pyk_TDP -1 Pyk_UDP -1 TX_FRUPTS -1 TX_GLCPTS = 0
 PHE_c: 2 Apts_Phe + 1 AspC3 + 1 TX_AROP12 + 1 TX_AROP22 -3.046751e-01 Growth = 0
 PhePhe_c: -1 Apts_Phe + 1 TX_OppB_Phe = 0
 PHPYR_c: -1 AspC3 = 0
 PI_c: 1 ArgF -1 AtpA + 1 BPNT -1 DeoD1 -1 DeoD2 -1 DeoD3 + 1 DeoD4 + 1 DeoD5 + 1 DeoD6 + 1 DeoD7 -1 DeoD8 + 1 FthC + 1 FTHFS -1 GapA + 1 MetK + 1 NADR2 + 1 PcbX -1 Pdp1 -1 Pdp2 -1 Pdp3 + 1 PgpA1 + 1 PgpA2 + 1 PgpA3 + 1 PhnX + 1 PncB + 2 Ppa + 1 Ppc + 1 pSERpp -1 Pta + 1 pTHRpp + 1 pTYRpp + 1 TX_Cbi_Co + 1 TX_Cbi_Ni + 1 TX_GL + 1 TX_GL3P + 1 TX_GLNHPQ + 1 TX_MdlB_doxorubicin + 1 TX_MGTA_CA + 1 TX_MGTA_CU + 1 TX_MGTA_FE2 + 1 TX_MGTA_FE3 + 1 TX_MGTA_MG + 1 TX_MGTA_MN + 1 TX_MGTA_MOBD + 1 TX_MGTA_ZN + 1 TX_OppB_Ala + 1 TX_OppB_Arg + 1 TX_OppB_Asn + 1 TX_OppB_Asp + 1 TX_OppB_Cys + 1 TX_OppB_Gln + 1 TX_OppB_Glu + 1 TX_OppB_Gly + 1 TX_OppB_His + 1 TX_OppB_Ile + 1 TX_OppB_Leu + 1 TX_OppB_Lys + 1 TX_OppB_Met + 1 TX_OppB_Phe + 1 TX_OppB_Pro + 1 TX_OppB_Ser + 1 TX_OppB_Thr + 1 TX_OppB_Trp + 1 TX_OppB_Tyr + 1 TX_OppB_Val + 1 TX_PHNT + 2 TX_PIUP1 + 1 TX_RIBUP + 1 TX_SalX_nisin + 1 TX_SO4 + 1 TX_SPMDUP1 + 1 TX_SPRDUP2 + 1 TX_SunT_nisin -2 TX_UHPT1 -2 TX_UHPT2 -2 TX_UHPT3 + 1 UshA_AMP + 1 UshA_aza6dTMP + 1 UshA_aza8dAMP + 1 UshA_c5dUMP + 1 UshA_carbamoylmethyl6dAMP + 1 UshA_cmnm5s2UMP + 2 UshA_cyclobutane_dCMPdCMP + 2 UshA_cyclobutane_dCMPdTMP + 2 UshA_cyclobutane_dTMPdCMP + 2 UshA_cyclobutane_dTMPdTMP + 1 UshA_cyclodAMP + 1 UshA_cyclodGMP + 1 UshA_dCMP + 2 UshA_dCMP64dCMP + 2 UshA_dCMP64dCMP_dewar + 2 UshA_dCMP64dTMP + 2 UshA_dCMP64dTMP_dewar + 1 UshA_dCMP_GLYC + 1 UshA_DHdTMP + 1 UshA_DHdUMP + 1 UshA_dhpdUMP + 1 UshA_dINSMP + 1 UshA_dTMP + 2 UshA_dTMP64dCMP + 2 UshA_dTMP64dCMP_dewar + 2 UshA_dTMP64dTMP + 2 UshA_dTMP64dTMP_dewar + 1 UshA_dTMP_GLY + 1 UshA_dUMP + 1 UshA_dXANMP + 1 UshA_e1dAMP + 1 UshA_e1dGMP + 1 UshA_e3dAMP + 1 UshA_e3dCMP + 1 UshA_e3dGMP + 1 UshA_e3dTMP + 1 UshA_e4dCMP + 1 UshA_e5dCMP + 1 UshA_e6dAMP + 1 UshA_e7dAMP + 1 UshA_e7dGMP + 1 UshA_e8dGMP + 1 UshA_epsilon1N2dGMP + 1 UshA_epsilon1N6dAMP + 1 UshA_epsilon3N4dCMP + 1 UshA_epsilonN23dGMP + 1 UshA_f5dUMP + 1 UshA_FAPydAMP + 1 UshA_FAPydGMP + 1 UshA_glchm5dCMP + 1 UshA_glchm5dUMP + 1 UshA_gludTMP + 1 UshA_GmMP + 1 UshA_GMP + 1 UshA_hm3dCMP + 1 UshA_hm5dCMP + 1 UshA_hm5dUMP + 1 UshA_hm5hydantoin_dRiboseMP + 1 UshA_ho5dCMP + 1 UshA_ho5dUMP + 1 UshA_ho5ho6dTMP + 1 UshA_ho5ho6dUMP + 1 UshA_ho5hydantoin_dRiboseMP + 1 UshA_IMP + 1 UshA_k2CMP + 1 UshA_m1dAMP + 1 UshA_m1dGMP + 1 UshA_m1GMP + 1 UshA_m2GMP + 1 UshA_m3dAMP + 1 UshA_m3dCMP + 1 UshA_m3dGMP + 1 UshA_m3dTMP + 1 UshA_m4dCMP + 1 UshA_m5CMP + 1 UshA_m5dCMP + 1 UshA_m62AMP + 1 UshA_m6AMP + 1 UshA_m6dAMP + 1 UshA_m7dAMP + 1 UshA_m7dGMP + 1 UshA_m7GMP + 1 UshA_m8dGMP + 1 UshA_n2dAMP + 1 UshA_O2edCMP + 1 UshA_O2edTMP + 1 UshA_O2mdCMP + 1 UshA_O2mdTMP + 1 UshA_O4edTMP + 1 UshA_O4mdTMP + 1 UshA_O6edGMP + 1 UshA_O6mdGMP + 1 UshA_oxo2dAMP + 1 UshA_oxo8dAMP + 1 UshA_oxo8dGMP + 1 UshA_oxoethyl7dGMP + 1 UshA_PSIURIMP + 1 UshA_putdTMP + 1 UshA_s2UMP + 1 UshA_s4UMP + 1 UshA_UmMP + 1 YbiV_deoxyG6P + 1 YbiV_F6P + 1 YbiV_G6P + 1 YbiV_GL2P + 1 YbiV_MAN6P + 1 YbiV_R5P + 1 YbiV_SBT6P + 1 YbiV_SNGLYP -1 PIix + 2.652520e+01 Growth = 0
 PL5P_c: -1 GlyA1 -1 GlyA2 + 1 PYDXK -1.210526e-03 Growth = 0
 PPAA_c: 1 PhnW -1 PhnX = 0
 PPGPP_c: 1 RelA1 -8.780461e-04 Growth = 0
 PPHG_c: -1 HemG = 0
 PPI_c: 1 Aas1 + 1 Aas2 + 1 Aas3 + 1 Aas4 + 1 Aas5 + 1 Aas6 + 1 Aas7 + 1 Apt + 1 CdsA1 + 1 CdsA2 + 1 CdsA3 + 1 CdsA4 + 1 CdsA5 + 1 CdsA6 + 1 CdsA7 + 2 diadenylate_cyclase + 1 GalU1 + 1 Hpt1 + 1 Hpt2 + 1 LplA + 1 MetK + 1 NADD1 + 1 NADD2 + 1 NadE + 1 NADR1 + 1 PncB -1 Ppa + 1 RibF2 + 1 Upp -1 PPIix = 0
 PPIX_c: 1 HemG = 0
 PRO_c: 2 Apts_Pro -1.626517e-01 Growth = 0
 proflavine_c: -1 TX_NormM_proflavine = 0
 ProPro_c: -1 Apts_Pro + 1 TX_OppB_Pro = 0
 PRPP_c: -1 Apt -1 Hpt1 -1 Hpt2 -1 PncB + 1 PrsA -1 Upp = 0
 pSER_c: -1 pSERpp -1 pSERix = 0
 PSIURI_c: 1 TX_PSIURI + 1 UshA_PSIURIMP = 0
 PSIURIMP_c: -1 UshA_PSIURIMP -1 PSIURIMPix = 0
 pTHR_c: -1 pTHRpp -1 pTHRix = 0
 PTRC_c: 1 TX_SPRDUP2 -1.100900e-01 Growth = 0
 pTYR_c: -1 pTYRpp -1 pTYRix = 0
 putdT_c: 1 TX_putdT + 1 UshA_putdTMP = 0
 putdTMP_c: -1 UshA_putdTMP = 0
 PYDX_c: -1 PYDXK + 1 TXPYDX = 0
 PYR_c: -1 AceE + 1 DhaK + 1 GlyA1 + 1 GlyA2 -1 LdhA -1 PcbX + 1 Pyk_ADP + 1 Pyk_CDP + 1 Pyk_DADP + 1 Pyk_DCDP + 1 Pyk_DGDP + 1 Pyk_DTDP + 1 Pyk_DUDP + 1 Pyk_GDP + 1 Pyk_IDP + 1 Pyk_TDP + 1 Pyk_UDP + 1 TX_FRUPTS + 1 TX_GLCPTS = 0
 R1P_c: -1 DeoD4 -1 DeoD5 -1 DeoD6 -1 DeoD7 + 1 Pdp3 = 0
 R5P_c: -1 PrsA + 1 RpiB -1 Tkt1 -1 YbiV_R5P = 0
 RBFLVRD_c: 1 Fre3 + 1 Fre4 = 0
 rhodamine6G_c: -1 TX_NormM_rhodamine6G = 0
 RIB_c: 1 TX_RIBUP + 1 YbiV_R5P = 0
 RIBFLV_c: -1 Fre3 -1 Fre4 -1 RibF1 + 1 TX_RIBFLV -1.210526e-03 Growth = 0
 RL5P_c: -1 Rpe -1 RpiB = 0
 s2U_c: 1 TX_s2U + 1 UshA_s2UMP = 0
 s2UMP_c: -1 UshA_s2UMP = 0
 s4U_c: 1 TX_s4U + 1 UshA_s4UMP = 0
 s4UMP_c: -1 UshA_s4UMP -1 s4UMPix = 0
 S7P_c: -1 TalA + 1 Tkt1 = 0
 SALA_c: -1 CsdA = 0
 SBT_c: 1 YbiV_SBT6P = 0
 SBT6P_c: -1 YbiV_SBT6P = 0
 SER_c: 2 Apts_Ser + 1 GlpQ3 -1 GlyA3 + 1 MG385r3 + 1 pSERpp -3.199191e-01 Growth = 0
 SerSer_c: -1 Apts_Ser + 1 TX_OppB_Ser = 0
 SNGLYP_c: -1 SNGLYPOr -1 YbiV_SNGLYP -1 SNGLYPix = 0
 SO2_c: 1 CsdA = 0
 SO4_c: 1 TX_SO4 -3.836400e-03 Growth = 0
 SPMD_c: 1 TX_SPMDUP1 -2.231641e-02 Growth = 0
 T3P1_c: 1 DeoC + 1 Fba -1 GapA -1 TalA + 1 Tkt1 + 1 Tkt2 + 1 TpiA = 0
 T3P2_c: 1 DhaK + 1 Fba -1 GpsA + 1 SNGLYPOr -1 TpiA = 0
 TAG16P_c: 1 FruK2 = 0
 TAG6P_c: -1 FruK2 = 0
 tBHP_c: -2 Ohr_tBHP -2 OsmC_tBHP = 0
 tBOH_c: 2 Ohr_tBHP + 2 OsmC_tBHP = 0
 TDECEG3P_c: -1 PlsC3 + 1 PlsX3 = 0
 TDECG3P_c: -1 PlsC2 + 1 PlsX2 = 0
 TDP_c: -1 Pyk_TDP = 0
 tetraphenylphosphonium_c: -1 TX_NormM_tetraphenylphosphonium = 0
 THF_c: 1 DhfR1 -1 FTHFS -1 GlyA3 + 1 MetS + 1 TX_THF -1 THFix = 0
 THM_c: -1 TMDPK + 1 TXTHM = 0
 THMP_c: -1 ThiL = 0
 THR_c: 2 Apts_Thr -1 GlyA4 + 1 pTHRpp -3.021227e-01 Growth = 0
 ThrThr_c: -1 Apts_Thr + 1 TX_OppB_Thr = 0
 THY_c: 1 Pdp2 + 1 TX_THY -1 THYix = 0
 TPP_c: 1 ThiL + 1 TMDPK -1.210526e-03 Growth = 0
 TRP_c: 2 Apts_Trp + 1 TX_AROP11 + 1 TX_AROP21 -4.515552e-02 Growth = 0
 TrpTrp_c: -1 Apts_Trp + 1 TX_OppB_Trp = 0
 TTDCA_c: -1 Aas2 + 1 TX_TTDCA = 0
 TTDCEA_c: -1 Aas3 + 1 TX_TTDCEA = 0
 TTP_c: 1 Pyk_TDP = 0
 TYR_c: 2 Apts_Tyr + 1 AspC2 + 1 pTYRpp + 1 TX_AROP13 + 1 TX_AROP23 -1.698999e-01 Growth = 0
 TyrTyr_c: -1 Apts_Tyr + 1 TX_OppB_Tyr = 0
 U23CMP_c: -1 CpdB1 = 0
 U3MP_c: 1 CpdB1 = 0
 UDP_c: 1 CmkB3 + 1 glycosyltransferase_Ceramide16_GalCer16 + 1 glycosyltransferase_Ceramide16_GlcCer16 + 1 glycosyltransferase_DAG160_GalDAG160 + 1 glycosyltransferase_DAG160_GlcDAG160 + 1 glycosyltransferase_DAG161_GalDAG161 + 1 glycosyltransferase_DAG161_GlcDAG161 + 1 glycosyltransferase_DAG181_GalDAG181 + 1 glycosyltransferase_DAG181_GlcDAG181 + 1 glycosyltransferase_GalCer16_GalGalCer16 + 1 glycosyltransferase_GalCer16_GlcGalCer16 + 1 glycosyltransferase_GalDAG160_GalGalDAG160 + 1 glycosyltransferase_GalDAG160_GlcGalDAG160 + 1 glycosyltransferase_GalDAG161_GalGalDAG161 + 1 glycosyltransferase_GalDAG161_GlcGalDAG161 + 1 glycosyltransferase_GalDAG181_GalGalDAG181 + 1 glycosyltransferase_GalDAG181_GlcGalDAG181 + 1 glycosyltransferase_GalGalDAG160_GalGalGalDAG160 + 1 glycosyltransferase_GalGalDAG160_GlcGalGalDAG160 + 1 glycosyltransferase_GalGalDAG161_GalGalGalDAG161 + 1 glycosyltransferase_GalGalDAG161_GlcGalGalDAG161 + 1 glycosyltransferase_GalGalDAG181_GalGalGalDAG181 + 1 glycosyltransferase_GalGalDAG181_GlcGalGalDAG181 + 1 glycosyltransferase_GalGlcCer16_GalGalGlcCer16 + 1 glycosyltransferase_GalGlcDAG160_GalGalGlcDAG160 + 1 glycosyltransferase_GalGlcDAG160_GlcGalGlcDAG160 + 1 glycosyltransferase_GalGlcDAG161_GalGalGlcDAG161 + 1 glycosyltransferase_GalGlcDAG161_GlcGalGlcDAG161 + 1 glycosyltransferase_GalGlcDAG181_GalGalGlcDAG181 + 1 glycosyltransferase_GalGlcDAG181_GlcGalGlcDAG181 + 1 glycosyltransferase_GlcCer16_GalGlcCer16 + 1 glycosyltransferase_GlcDAG160_GalGlcDAG160 + 1 glycosyltransferase_GlcDAG160_GlcGlcDAG160 + 1 glycosyltransferase_GlcDAG161_GalGlcDAG161 + 1 glycosyltransferase_GlcDAG161_GlcGlcDAG161 + 1 glycosyltransferase_GlcDAG181_GalGlcDAG181 + 1 glycosyltransferase_GlcDAG181_GlcGlcDAG181 + 1 glycosyltransferase_GlcGlcDAG160_GlcGlcGlcDAG160 + 1 glycosyltransferase_GlcGlcDAG161_GlcGlcGlcDAG161 + 1 glycosyltransferase_GlcGlcDAG181_GlcGlcGlcDAG181 -2 NrdEF_UDP + 1 PfkA5 -1 Pyk_UDP = 0
 UDPG_c: 1 GalE + 1 GalU1 -1 glycosyltransferase_Ceramide16_GlcCer16 -1 glycosyltransferase_DAG160_GlcDAG160 -1 glycosyltransferase_DAG161_GlcDAG161 -1 glycosyltransferase_DAG181_GlcDAG181 -1 glycosyltransferase_GalCer16_GlcGalCer16 -1 glycosyltransferase_GalDAG160_GlcGalDAG160 -1 glycosyltransferase_GalDAG161_GlcGalDAG161 -1 glycosyltransferase_GalDAG181_GlcGalDAG181 -1 glycosyltransferase_GalGalDAG160_GlcGalGalDAG160 -1 glycosyltransferase_GalGalDAG161_GlcGalGalDAG161 -1 glycosyltransferase_GalGalDAG181_GlcGalGalDAG181 -1 glycosyltransferase_GalGlcDAG160_GlcGalGlcDAG160 -1 glycosyltransferase_GalGlcDAG161_GlcGalGlcDAG161 -1 glycosyltransferase_GalGlcDAG181_GlcGalGlcDAG181 -1 glycosyltransferase_GlcDAG160_GlcGlcDAG160 -1 glycosyltransferase_GlcDAG161_GlcGlcDAG161 -1 glycosyltransferase_GlcDAG181_GlcGlcDAG181 -1 glycosyltransferase_GlcGlcDAG160_GlcGlcGlcDAG160 -1 glycosyltransferase_GlcGlcDAG161_GlcGlcGlcDAG161 -1 glycosyltransferase_GlcGlcDAG181_GlcGlcGlcDAG181 = 0
 UDPGAL_c: -1 GalE -1 Glf -1 glycosyltransferase_Ceramide16_GalCer16 -1 glycosyltransferase_DAG160_GalDAG160 -1 glycosyltransferase_DAG161_GalDAG161 -1 glycosyltransferase_DAG181_GalDAG181 -1 glycosyltransferase_GalCer16_GalGalCer16 -1 glycosyltransferase_GalDAG160_GalGalDAG160 -1 glycosyltransferase_GalDAG161_GalGalDAG161 -1 glycosyltransferase_GalDAG181_GalGalDAG181 -1 glycosyltransferase_GalGalDAG160_GalGalGalDAG160 -1 glycosyltransferase_GalGalDAG161_GalGalGalDAG161 -1 glycosyltransferase_GalGalDAG181_GalGalGalDAG181 -1 glycosyltransferase_GalGlcCer16_GalGalGlcCer16 -1 glycosyltransferase_GalGlcDAG160_GalGalGlcDAG160 -1 glycosyltransferase_GalGlcDAG161_GalGalGlcDAG161 -1 glycosyltransferase_GalGlcDAG181_GalGalGlcDAG181 -1 glycosyltransferase_GlcCer16_GalGlcCer16 -1 glycosyltransferase_GlcDAG160_GalGlcDAG160 -1 glycosyltransferase_GlcDAG161_GalGlcDAG161 -1 glycosyltransferase_GlcDAG181_GalGlcDAG181 = 0
 UDPGALF_c: 1 Glf -8.222598e-03 Growth = 0
 Um_c: 1 TX_Um + 1 UshA_UmMP = 0
 UmMP_c: -1 UshA_UmMP -1 UmMPix = 0
 UMP_c: -1 CmkB3 + 1 Udk1 + 1 Upp -1 UMPix = 0
 URA_c: 1 DeoD8 + 1 Pdp1 + 1 Pdp3 + 1 TX_URA -1 Upp -1 URAix = 0
 URI_c: 1 Cdd1 -1 Pdp3 -1 Udk1 = 0
 UTP_c: -1 GalU1 -1 PfkA5 + 1 Pyk_UDP -1 UTPix -1.047582e-01 Growth = 0
 VAL_c: 2 Apts_Val -3.642338e-01 Growth = 0
 ValVal_c: -1 Apts_Val + 1 TX_OppB_Val = 0
 vinyl_chloride_c: -1 vinyl_chloride_activation_1 = 0
 X5P_c: 1 Rpe -1 Tkt1 -1 Tkt2 = 0
 XAN_c: -1 DeoD7 = 0
 XTSN_c: 1 DeoD7 = 0
 ZN_c: 1 TX_MGTA_ZN -4.253442e-03 Growth = 0
 MG_124_MONOMER_c: 1 MsrA + 1 MsrB -2 NrdEF_maintenance -2 Ohr_rx -2 OsmC_rx + 1 TrxB = 0
 MG_287_MONOMER_c: -1 AcpS = 0
 MG_124_MONOMER_ox_c: -1 MsrA -1 MsrB + 2 NrdEF_maintenance + 2 Ohr_rx + 2 OsmC_rx -1 TrxB = 0
 MG_229_231_TETRAMER_c: -1 NrdEF_ADP -1 NrdEF_CDP -1 NrdEF_GDP + 1 NrdEF_maintenance -1 NrdEF_UDP = 0
 MG_229_231_TETRAMER_ox_c: 1 NrdEF_ADP + 1 NrdEF_CDP + 1 NrdEF_GDP -1 NrdEF_maintenance + 1 NrdEF_UDP = 0
 MG_287_MONOMER_ACP_c: -1 Aas1 -1 Aas2 -1 Aas3 -1 Aas4 -1 Aas5 -1 Aas6 -1 Aas7 + 1 AcpS + 1 PlsC1 + 1 PlsC2 + 1 PlsC3 + 1 PlsC4 + 1 PlsC5 + 1 PlsC6 + 1 PlsC7 + 1 PlsX1 + 1 PlsX2 + 1 PlsX3 + 1 PlsX4 + 1 PlsX5 + 1 PlsX6 + 1 PlsX7 = 0
 MG_287_MONOMER_ddcaACP_c: 1 Aas1 -1 PlsC1 -1 PlsX1 = 0
 MG_287_MONOMER_hdeACP_c: 1 Aas5 -1 PlsC5 -1 PlsX5 = 0
 MG_287_MONOMER_myrsACP_c: 1 Aas2 -1 PlsC2 -1 PlsX2 = 0
 MG_287_MONOMER_ocdcaACP_c: 1 Aas6 -1 PlsC6 -1 PlsX6 = 0
 MG_287_MONOMER_octeACP_c: 1 Aas7 -1 PlsC7 -1 PlsX7 = 0
 MG_287_MONOMER_palmACP_c: 1 Aas4 -1 PlsC4 -1 PlsX4 = 0
 MG_287_MONOMER_tdeACP_c: 1 Aas3 -1 PlsC3 -1 PlsX3 = 0
 MG_427_DIMER_c: -1 OsmC_CHP -1 OsmC_H2O2 + 1 OsmC_rx -1 OsmC_tBHP = 0
 MG_427_DIMER_ox_c: 1 OsmC_CHP + 1 OsmC_H2O2 -1 OsmC_rx + 1 OsmC_tBHP = 0
 MG_454_DIMER_c: -1 Ohr_CHP -1 Ohr_H2O2 + 1 Ohr_rx -1 Ohr_tBHP = 0
 MG_454_DIMER_ox_c: 1 Ohr_CHP + 1 Ohr_H2O2 -1 Ohr_rx + 1 Ohr_tBHP = 0
 Internal_Exchange_Constraint_01: 1 ADPix + 1 AMPix + 1 ATPix = 0
 Internal_Exchange_Constraint_02: 1 CMPix + 1 CTPix = 0
 Internal_Exchange_Constraint_03: 1 GDPix + 1 GMPix + 1 GTPix = 0
 Internal_Exchange_Constraint_04: 1 UMPix + 1 UTPix = 0
 Internal_Exchange_Constraint_05: 1 FTHF10ix + 1 THFix = 0
 Internal_Exchange_Constraint_06: 1 FMETix + 1 FORix + 1 FTHF10ix = 0
 Internal_Exchange_Constraint_07: 1 FMETix + 1 METix = 0
 Biomass: 1 Growth + 1 Bmssex = 0
 gamma_radiation_e: -1 hydroxyl_radical_generation_1 -1 hydroxyl_radical_generation_2 + 1 gamma_radiationex = 0
 A23CMP_e: 1 A23CMPex = 0
 A3MP_e: 1 A3MPex = 0
 AC_e: -1 TX_ACUP + 1 ACex = 0
 ACAL_e: 1 TX_ACAL + 1 ACALex = 0
 ACCOA_e: 1 ACCOAex = 0
 acriflavine_e: 1 TX_NormM_acriflavine + 1 acriflavineex = 0
 ACTP_e: 1 ACTPex = 0
 AD_e: -1 TX_AD + 1 ADex = 0
 ADN_e: 1 ADNex = 0
 ADP_e: 1 ADPex = 0
 AEPP_e: -1 TX_PHNT + 1 AEPPex = 0
 AHCYS_e: 1 AHCYSex = 0
 AKG_e: -1 TX_AKG + 1 AKGex = 0
 ALA_e: 1 ALAex = 0
 AlaAla_e: -1 TX_OppB_Ala + 1 AlaAlaex = 0
 AMET_e: 1 AMETex = 0
 aminodimethylaniline_e: 1 aminodimethylanilineex = 0
 AMP_e: 1 HinT_AMP_Lysine_MG132 + 1 HinT_AMP_Mor_MG132 + 1 HinT_AMP_NH2_MG132 + 1 AMPex = 0
 AMP_Lysine_e: -1 HinT_AMP_Lysine_MG132 + 1 AMP_Lysineex = 0
 AMP_Mor_e: -1 HinT_AMP_Mor_MG132 + 1 AMP_Morex = 0
 AMP_NH2_e: -1 HinT_AMP_NH2_MG132 + 1 AMP_NH2ex = 0
 anthranilate_e: 1 anthranilateex = 0
 ARG_e: 1 ARGex = 0
 ArgArg_e: -1 TX_OppB_Arg + 1 ArgArgex = 0
 ASN_e: 1 ASNex = 0
 AsnAsn_e: -1 TX_OppB_Asn + 1 AsnAsnex = 0
 ASP_e: 1 ASPex = 0
 AspAsp_e: -1 TX_OppB_Asp + 1 AspAspex = 0
 ATP_e: 1 ATPex = 0
 aza6dT_e: -1 TX_aza6dT + 1 aza6dTex = 0
 aza6dTMP_e: 1 aza6dTMPex = 0
 aza8dA_e: -1 TX_aza8dA + 1 aza8dAex = 0
 aza8dAMP_e: 1 aza8dAMPex = 0
 benzalkonium_e: 1 TX_NormM_benzalkonium + 1 benzalkoniumex = 0
 berberine_e: 1 TX_NormM_berberine + 1 berberineex = 0
 BUT_e: 2 LIPASE_DIBUTYRIN_MG310 + 2 LIPASE_DIBUTYRIN_MG327 + 2 LIPASE_DIBUTYRIN_MG344 + 1 LIPASE_MONOBUTYRIN_MG310 + 1 LIPASE_MONOBUTYRIN_MG327 + 1 LIPASE_MONOBUTYRIN_MG344 + 3 LIPASE_TRIBUTYRIN_MG310 + 3 LIPASE_TRIBUTYRIN_MG327 + 3 LIPASE_TRIBUTYRIN_MG344 + 1 BUTex = 0
 C23CMP_e: 1 C23CMPex = 0
 C3MP_e: 1 C3MPex = 0
 c5dU_e: -1 TX_c5dU + 1 c5dUex = 0
 c5dUMP_e: 1 c5dUMPex = 0
 CA_e: -1 TX_MGTA_CA + 1 CAex = 0
 CAP_e: -1 TX_CAP + 1 CAPex = 0
 carbamoylmethyl6dA_e: -1 TX_carbamoylmethyl6dA + 1 carbamoylmethyl6dAex = 0
 carbamoylmethyl6dAMP_e: 1 carbamoylmethyl6dAMPex = 0
 CDP_e: 1 CDPex = 0
 CDPDG120_e: 1 CDPDG120ex = 0
 CDPDG140_e: 1 CDPDG140ex = 0
 CDPDG141_e: 1 CDPDG141ex = 0
 CDPDG160_e: 1 CDPDG160ex = 0
 CDPDG161_e: 1 CDPDG161ex = 0
 CDPDG180_e: 1 CDPDG180ex = 0
 CDPDG181_e: 1 CDPDG181ex = 0
 Ceramide16_e: 1 Ceramide16ex = 0
 CH3OH_e: 1 LIPASE_METHYL_OCDCEA_MG310 + 1 LIPASE_METHYL_OCDCEA_MG327 + 1 LIPASE_METHYL_OCDCEA_MG344 + 1 CH3OHex = 0
 chloroacetaldehyde_e: 1 chloroacetaldehydeex = 0
 chloroethylene_oxide_e: 1 chloroethylene_oxideex = 0
 CHOL_e: 1 CHOLex = 0
 CHP_e: 1 CHPex = 0
 ciprofloxacin_e: 1 TX_NormM_ciprofloxacin + 1 ciprofloxacinex = 0
 CITR_e: 1 CITRex = 0
 CL_e: -1 TX_CL + 1 CLex = 0
 CL160_e: 1 CL160ex = 0
 CL161_e: 1 CL161ex = 0
 CL181_e: 1 CL181ex = 0
 cmnm5s2U_e: -1 TX_cmnm5s2U + 1 cmnm5s2Uex = 0
 cmnm5s2UMP_e: 1 cmnm5s2UMPex = 0
 CMP_e: 1 CMPex = 0
 CO_e: -1 TX_Cbi_Co + 1 COex = 0
 CO2_e: -1 TX_CO2 + 1 CO2ex = 0
 COA_e: -1 TX_COA + 1 COAex = 0
 crystalViolet_e: 1 TX_NormM_crystalViolet + 1 crystalVioletex = 0
 CTP_e: 1 CTPex = 0
 CU_e: -1 TX_MGTA_CU + 1 CUex = 0
 cyclobutane_dCdC_e: -1 TX_cyclobutane_dCdC + 1 cyclobutane_dCdCex = 0
 cyclobutane_dCdT_e: -1 TX_cyclobutane_dCdT + 1 cyclobutane_dCdTex = 0
 cyclobutane_dCMPdCMP_e: 1 cyclobutane_dCMPdCMPex = 0
 cyclobutane_dCMPdTMP_e: 1 cyclobutane_dCMPdTMPex = 0
 cyclobutane_dTdC_e: -1 TX_cyclobutane_dTdC + 1 cyclobutane_dTdCex = 0
 cyclobutane_dTdT_e: -1 TX_cyclobutane_dTdT + 1 cyclobutane_dTdTex = 0
 cyclobutane_dTMPdCMP_e: 1 cyclobutane_dTMPdCMPex = 0
 cyclobutane_dTMPdTMP_e: 1 cyclobutane_dTMPdTMPex = 0
 cyclodA_e: -1 TX_cyclodA + 1 cyclodAex = 0
 cyclodAMP_e: 1 cyclodAMPex = 0
 cyclodG_e: -1 TX_cyclodG + 1 cyclodGex = 0
 cyclodGMP_e: 1 cyclodGMPex = 0
 CYS_e: 1 CYSex = 0
 CysCys_e: -1 TX_OppB_Cys + 1 CysCysex = 0
 CYTD_e: -1 TX_CYTD + 1 CYTDex = 0
 c_di_AMP_e: 1 c_di_AMPex = 0
 DA_e: -1 TX_DA + 1 DAex = 0
 DADP_e: 1 DADPex = 0
 DAG160_e: 1 DAG160ex = 0
 DAG161_e: 1 DAG161ex = 0
 DAG181_e: 1 DAG181ex = 0
 DALA_e: 1 DALAex = 0
 DAMP_e: 2 TX_ComEC_dApdAp + 1 DAMPex = 0
 dApdAp_e: -2 TX_ComEC_dApdAp + 1 dApdApex = 0
 DATP_e: 1 DATPex = 0
 DC_e: 1 DCex = 0
 dC64dC_e: -1 TX_dC64dC + 1 dC64dCex = 0
 dC64dC_dewar_e: -1 TX_dC64dC_dewar + 1 dC64dC_dewarex = 0
 dC64dT_e: -1 TX_dC64dT + 1 dC64dTex = 0
 dC64dT_dewar_e: -1 TX_dC64dT_dewar + 1 dC64dT_dewarex = 0
 DCDP_e: 1 DCDPex = 0
 DCMP_e: 2 TX_ComEC_dCpdCp + 1 DCMPex = 0
 dCMP64dCMP_e: 1 dCMP64dCMPex = 0
 dCMP64dCMP_dewar_e: 1 dCMP64dCMP_dewarex = 0
 dCMP64dTMP_e: 1 dCMP64dTMPex = 0
 dCMP64dTMP_dewar_e: 1 dCMP64dTMP_dewarex = 0
 dCMP_GLYC_e: 1 dCMP_GLYCex = 0
 dCpdCp_e: -2 TX_ComEC_dCpdCp + 1 dCpdCpex = 0
 DCTP_e: 1 DCTPex = 0
 dC_GLYC_e: -1 TX_dC_GLYC + 1 dC_GLYCex = 0
 DDCA_e: 3 LIPASE_TRILAURIN_MG310 + 3 LIPASE_TRILAURIN_MG327 + 3 LIPASE_TRILAURIN_MG344 -1 TX_DDCA + 1 DDCAex = 0
 DDECG3P_e: 1 DDECG3Pex = 0
 deoxyG6P_e: 1 deoxyG6Pex = 0
 deoxyGLC_e: 1 deoxyGLCex = 0
 DG_e: 1 DGex = 0
 DGDP_e: 1 DGDPex = 0
 DGMP_e: 2 TX_ComEC_dGpdGp + 1 DGMPex = 0
 dGpdGp_e: -2 TX_ComEC_dGpdGp + 1 dGpdGpex = 0
 DGTP_e: 1 DGTPex = 0
 DHA_e: 1 DHAex = 0
 DHdT_e: -1 TX_DHdT + 1 DHdTex = 0
 DHdTMP_e: 1 DHdTMPex = 0
 DHdU_e: -1 TX_DHdU + 1 DHdUex = 0
 DHdUMP_e: 1 DHdUMPex = 0
 DHF_e: 1 DHFex = 0
 dhpdU_e: -1 TX_dhpdU + 1 dhpdUex = 0
 dhpdUMP_e: 1 dhpdUMPex = 0
 diacylglycerolCys_e: 1 TX_diacylglycerolCys + 1 diacylglycerolCysex = 0
 DIBUTYRIN_e: -1 LIPASE_DIBUTYRIN_MG310 -1 LIPASE_DIBUTYRIN_MG327 -1 LIPASE_DIBUTYRIN_MG344 + 1 DIBUTYRINex = 0
 dINS_e: -1 TX_dINS + 1 dINSex = 0
 dINSMP_e: 1 dINSMPex = 0
 DMBA_e: 1 DMBAex = 0
 doxorubicin_e: 1 TX_MdlB_doxorubicin + 1 TX_NormM_doxorubicin + 1 doxorubicinex = 0
 DPCOA_e: 1 DPCOAex = 0
 DPG_e: 1 DPGex = 0
 DR1P_e: 1 DR1Pex = 0
 DR5P_e: 1 DR5Pex = 0
 DT_e: 1 DTex = 0
 dT64dC_e: -1 TX_dT64dC + 1 dT64dCex = 0
 dT64dC_dewar_e: -1 TX_dT64dC_dewar + 1 dT64dC_dewarex = 0
 dT64dT_e: -1 TX_dT64dT + 1 dT64dTex = 0
 dT64dT_dewar_e: -1 TX_dT64dT_dewar + 1 dT64dT_dewarex = 0
 DTDP_e: 1 DTDPex = 0
 DTMP_e: 2 TX_ComEC_dTpdTp + 1 DTMPex = 0
 dTMP64dCMP_e: 1 dTMP64dCMPex = 0
 dTMP64dCMP_dewar_e: 1 dTMP64dCMP_dewarex = 0
 dTMP64dTMP_e: 1 dTMP64dTMPex = 0
 dTMP64dTMP_dewar_e: 1 dTMP64dTMP_dewarex = 0
 dTMP_GLYC_e: 1 dTMP_GLYCex = 0
 dTpdTp_e: -2 TX_ComEC_dTpdTp + 1 dTpdTpex = 0
 DTTP_e: 1 DTTPex = 0
 dT_GLYC_e: -1 TX_dT_GLYC + 1 dT_GLYCex = 0
 DU_e: 1 DUex = 0
 DUDP_e: 1 DUDPex = 0
 DUMP_e: 1 DUMPex = 0
 DUTP_e: 1 DUTPex = 0
 dXAN_e: -1 TX_dXAN + 1 dXANex = 0
 dXANMP_e: 1 dXANMPex = 0
 e1dA_e: -1 TX_e1dA + 1 e1dAex = 0
 e1dAMP_e: 1 e1dAMPex = 0
 e1dG_e: -1 TX_e1dG + 1 e1dGex = 0
 e1dGMP_e: 1 e1dGMPex = 0
 e3dA_e: -1 TX_e3dA + 1 e3dAex = 0
 e3dAMP_e: 1 e3dAMPex = 0
 e3dC_e: -1 TX_e3dC + 1 e3dCex = 0
 e3dCMP_e: 1 e3dCMPex = 0
 e3dG_e: -1 TX_e3dG + 1 e3dGex = 0
 e3dGMP_e: 1 e3dGMPex = 0
 e3dT_e: -1 TX_e3dT + 1 e3dTex = 0
 e3dTMP_e: 1 e3dTMPex = 0
 e4dC_e: -1 TX_e4dC + 1 e4dCex = 0
 e4dCMP_e: 1 e4dCMPex = 0
 E4P_e: 1 E4Pex = 0
 e5dC_e: -1 TX_e5dC + 1 e5dCex = 0
 e5dCMP_e: 1 e5dCMPex = 0
 e6dA_e: -1 TX_e6dA + 1 e6dAex = 0
 e6dAMP_e: 1 e6dAMPex = 0
 e7dA_e: -1 TX_e7dA + 1 e7dAex = 0
 e7dAMP_e: 1 e7dAMPex = 0
 e7dG_e: -1 TX_e7dG + 1 e7dGex = 0
 e7dGMP_e: 1 e7dGMPex = 0
 e8dG_e: -1 TX_e8dG + 1 e8dGex = 0
 e8dGMP_e: 1 e8dGMPex = 0
 enoxacin_e: 1 TX_NormM_enoxacin + 1 enoxacinex = 0
 epsilon1N2dG_e: -1 TX_epsilon1N2dG + 1 epsilon1N2dGex = 0
 epsilon1N2dGMP_e: 1 epsilon1N2dGMPex = 0
 epsilon1N6dA_e: -1 TX_epsilon1N6dA + 1 epsilon1N6dAex = 0
 epsilon1N6dAMP_e: 1 epsilon1N6dAMPex = 0
 epsilon3N4dC_e: -1 TX_epsilon3N4dC + 1 epsilon3N4dCex = 0
 epsilon3N4dCMP_e: 1 epsilon3N4dCMPex = 0
 epsilonN23dG_e: -1 TX_epsilonN23dG + 1 epsilonN23dGex = 0
 epsilonN23dGMP_e: 1 epsilonN23dGMPex = 0
 ethidium_e: 1 TX_NormM_ethidium + 1 ethidiumex = 0
 ETN_e: 1 ETNex = 0
 F1P_e: 1 F1Pex = 0
 f5dU_e: -1 TX_f5dU + 1 f5dUex = 0
 f5dUMP_e: 1 f5dUMPex = 0
 F6P_e: 1 F6Pex = 0
 FAD_e: 1 FADex = 0
 FADH2_e: 1 FADH2ex = 0
 FALD_e: -1 TX_FALD + 1 FALDex = 0
 FAPydA_e: -1 TX_FAPydA + 1 FAPydAex = 0
 FAPydAMP_e: 1 FAPydAMPex = 0
 FAPydG_e: -1 TX_FAPydG + 1 FAPydGex = 0
 FAPydGMP_e: 1 FAPydGMPex = 0
 FDP_e: 1 FDPex = 0
 FE2_e: -1 TX_MGTA_FE2 + 1 FE2ex = 0
 FE3_e: -1 TX_MGTA_FE3 + 1 FE3ex = 0
 FMET_e: 1 FMETex = 0
 FMN_e: 1 FMNex = 0
 FMNH2_e: 1 FMNH2ex = 0
 FOR_e: -1 TX_FOR + 1 FORex = 0
 FRU_e: -1 TX_FRUPTS + 1 FRUex = 0
 FTHF10_e: 1 FTHF10ex = 0
 FTHF5_e: 1 FTHF5ex = 0
 FUC1P_e: -1 TX_UHPT2 + 1 FUC1Pex = 0
 G1P_e: 1 G1Pex = 0
 G23CMP_e: 1 G23CMPex = 0
 G2P_e: 1 G2Pex = 0
 G3MP_e: 1 G3MPex = 0
 G3P_e: 1 G3Pex = 0
 G6P_e: -1 TX_UHPT1 + 1 G6Pex = 0
 GalCer16_e: 1 GalCer16ex = 0
 GalDAG160_e: 1 GalDAG160ex = 0
 GalDAG161_e: 1 GalDAG161ex = 0
 GalDAG181_e: 1 GalDAG181ex = 0
 GalGalCer16_e: 1 GalGalCer16ex = 0
 GalGalDAG160_e: 1 GalGalDAG160ex = 0
 GalGalDAG161_e: 1 GalGalDAG161ex = 0
 GalGalDAG181_e: 1 GalGalDAG181ex = 0
 GalGalGalDAG160_e: 1 GalGalGalDAG160ex = 0
 GalGalGalDAG161_e: 1 GalGalGalDAG161ex = 0
 GalGalGalDAG181_e: 1 GalGalGalDAG181ex = 0
 GalGalGlcCer16_e: 1 GalGalGlcCer16ex = 0
 GalGalGlcDAG160_e: 1 GalGalGlcDAG160ex = 0
 GalGalGlcDAG161_e: 1 GalGalGlcDAG161ex = 0
 GalGalGlcDAG181_e: 1 GalGalGlcDAG181ex = 0
 GalGlcCer16_e: 1 GalGlcCer16ex = 0
 GalGlcDAG160_e: 1 GalGlcDAG160ex = 0
 GalGlcDAG161_e: 1 GalGlcDAG161ex = 0
 GalGlcDAG181_e: 1 GalGlcDAG181ex = 0
 GDP_e: 1 GDPex = 0
 GDPTP_e: 1 GDPTPex = 0
 GL_e: 1 LIPASE_DIBUTYRIN_MG310 + 1 LIPASE_DIBUTYRIN_MG327 + 1 LIPASE_DIBUTYRIN_MG344 + 1 LIPASE_MONOBUTYRIN_MG310 + 1 LIPASE_MONOBUTYRIN_MG327 + 1 LIPASE_MONOBUTYRIN_MG344 + 1 LIPASE_TRIBUTYRIN_MG310 + 1 LIPASE_TRIBUTYRIN_MG327 + 1 LIPASE_TRIBUTYRIN_MG344 + 1 LIPASE_TRILAURIN_MG310 + 1 LIPASE_TRILAURIN_MG327 + 1 LIPASE_TRILAURIN_MG344 + 1 LIPASE_TRIMYRISTIN_MG310 + 1 LIPASE_TRIMYRISTIN_MG327 + 1 LIPASE_TRIMYRISTIN_MG344 + 1 LIPASE_TRIOLEIN_MG310 + 1 LIPASE_TRIOLEIN_MG327 + 1 LIPASE_TRIOLEIN_MG344 + 1 LIPASE_TRIPALMITIN_MG310 + 1 LIPASE_TRIPALMITIN_MG327 + 1 LIPASE_TRIPALMITIN_MG344 + 1 LIPASE_TRISTEARIN_MG310 + 1 LIPASE_TRISTEARIN_MG327 + 1 LIPASE_TRISTEARIN_MG344 + 1 LIPASE_TRI_HDCEA_IN_MG310 + 1 LIPASE_TRI_HDCEA_IN_MG327 + 1 LIPASE_TRI_HDCEA_IN_MG344 + 1 LIPASE_TRI_TTDCEA_IN_MG310 + 1 LIPASE_TRI_TTDCEA_IN_MG327 + 1 LIPASE_TRI_TTDCEA_IN_MG344 -1 TX_GL -1 TX_GLUP + 1 GLex = 0
 GL2P_e: 1 GL2Pex = 0
 GL3P_e: -1 TX_GL3P + 1 GL3Pex = 0
 GLC_e: -1 TX_GLCPTS + 1 GLCex = 0
 GlcCer16_e: 1 GlcCer16ex = 0
 GlcDAG160_e: 1 GlcDAG160ex = 0
 GlcDAG161_e: 1 GlcDAG161ex = 0
 GlcDAG181_e: 1 GlcDAG181ex = 0
 GlcGalCer16_e: 1 GlcGalCer16ex = 0
 GlcGalDAG160_e: 1 GlcGalDAG160ex = 0
 GlcGalDAG161_e: 1 GlcGalDAG161ex = 0
 GlcGalDAG181_e: 1 GlcGalDAG181ex = 0
 GlcGalGalDAG160_e: 1 GlcGalGalDAG160ex = 0
 GlcGalGalDAG161_e: 1 GlcGalGalDAG161ex = 0
 GlcGalGalDAG181_e: 1 GlcGalGalDAG181ex = 0
 GlcGalGlcDAG160_e: 1 GlcGalGlcDAG160ex = 0
 GlcGalGlcDAG161_e: 1 GlcGalGlcDAG161ex = 0
 GlcGalGlcDAG181_e: 1 GlcGalGlcDAG181ex = 0
 GlcGlcDAG160_e: 1 GlcGlcDAG160ex = 0
 GlcGlcDAG161_e: 1 GlcGlcDAG161ex = 0
 GlcGlcDAG181_e: 1 GlcGlcDAG181ex = 0
 GlcGlcGlcDAG160_e: 1 GlcGlcGlcDAG160ex = 0
 GlcGlcGlcDAG161_e: 1 GlcGlcGlcDAG161ex = 0
 GlcGlcGlcDAG181_e: 1 GlcGlcGlcDAG181ex = 0
 glchm5dC_e: -1 TX_glchm5dC + 1 glchm5dCex = 0
 glchm5dCMP_e: 1 glchm5dCMPex = 0
 glchm5dU_e: -1 TX_glchm5dU + 1 glchm5dUex = 0
 glchm5dUMP_e: 1 glchm5dUMPex = 0
 GLN_e: -1 TX_GLNHPQ + 1 GLNex = 0
 GlnGln_e: -1 TX_OppB_Gln + 1 GlnGlnex = 0
 GLU_e: 1 GLUex = 0
 gludT_e: -1 TX_gludT + 1 gludTex = 0
 gludTMP_e: 1 gludTMPex = 0
 GluGlu_e: -1 TX_OppB_Glu + 1 GluGluex = 0
 GLY_e: 1 GLYex = 0
 GlyGly_e: -1 TX_OppB_Gly + 1 GlyGlyex = 0
 Gm_e: -1 TX_Gm + 1 Gmex = 0
 GmMP_e: 1 GmMPex = 0
 GMP_e: 1 HinT_GMP_Lysine_MG132 + 1 HinT_GMP_Mor_MG132 + 1 GMPex = 0
 GMP_Lysine_e: -1 HinT_GMP_Lysine_MG132 + 1 GMP_Lysineex = 0
 GMP_Mor_e: -1 HinT_GMP_Mor_MG132 + 1 GMP_Morex = 0
 GN_e: -1 TX_GN + 1 GNex = 0
 GPCH_e: 1 GPCHex = 0
 GPE_e: 1 GPEex = 0
 GPG_e: 1 GPGex = 0
 GPI_e: 1 GPIex = 0
 GPS_e: 1 GPSex = 0
 GSN_e: 1 GSNex = 0
 GTP_e: 1 GTPex = 0
 H_e: -4 AtpA + 1 HinT_AMP_Lysine_MG132 + 2 HinT_AMP_Mor_MG132 + 2 HinT_AMP_NH2_MG132 + 1 HinT_GMP_Lysine_MG132 + 2 HinT_GMP_Mor_MG132 + 1 LIPASE_METHYL_OCDCEA_MG310 + 1 LIPASE_METHYL_OCDCEA_MG327 + 1 LIPASE_METHYL_OCDCEA_MG344 + 3 LIPASE_TRILAURIN_MG310 + 3 LIPASE_TRILAURIN_MG327 + 3 LIPASE_TRILAURIN_MG344 + 3 LIPASE_TRIMYRISTIN_MG310 + 3 LIPASE_TRIMYRISTIN_MG327 + 3 LIPASE_TRIMYRISTIN_MG344 + 3 LIPASE_TRIOLEIN_MG310 + 3 LIPASE_TRIOLEIN_MG327 + 3 LIPASE_TRIOLEIN_MG344 + 3 LIPASE_TRIPALMITIN_MG310 + 3 LIPASE_TRIPALMITIN_MG327 + 3 LIPASE_TRIPALMITIN_MG344 + 3 LIPASE_TRISTEARIN_MG310 + 3 LIPASE_TRISTEARIN_MG327 + 3 LIPASE_TRISTEARIN_MG344 + 3 LIPASE_TRI_HDCEA_IN_MG310 + 3 LIPASE_TRI_HDCEA_IN_MG327 + 3 LIPASE_TRI_HDCEA_IN_MG344 + 1 LIPASE_TRI_TTDCEA_IN_MG310 + 1 LIPASE_TRI_TTDCEA_IN_MG327 + 1 LIPASE_TRI_TTDCEA_IN_MG344 + 1 LIPASE_TWEEN80_MG310 + 1 LIPASE_TWEEN80_MG327 + 1 LIPASE_TWEEN80_MG344 -1 TXTHM -1 TX_ACUP -1 TX_AD -1 TX_AKG -1 TX_AROP11 -1 TX_AROP12 -1 TX_AROP13 -1 TX_AROP21 -1 TX_AROP22 -1 TX_AROP23 + 2 TX_CL + 1 TX_ComEC_dApdAp + 1 TX_ComEC_dCpdCp + 1 TX_ComEC_dGpdGp + 1 TX_ComEC_dTpdTp -1 TX_CYTD -1 TX_DA -1 TX_GN + 1 TX_NAUP2 + 1 TX_THY -1 TX_TRK -1 TX_URA + 1 Hex = 0
 H2O_e: -1 HinT_AMP_Lysine_MG132 -1 HinT_AMP_Mor_MG132 -1 HinT_AMP_NH2_MG132 -1 HinT_GMP_Lysine_MG132 -1 HinT_GMP_Mor_MG132 -2 LIPASE_DIBUTYRIN_MG310 -2 LIPASE_DIBUTYRIN_MG327 -2 LIPASE_DIBUTYRIN_MG344 -1 LIPASE_METHYL_OCDCEA_MG310 -1 LIPASE_METHYL_OCDCEA_MG327 -1 LIPASE_METHYL_OCDCEA_MG344 -1 LIPASE_MONOBUTYRIN_MG310 -1 LIPASE_MONOBUTYRIN_MG327 -1 LIPASE_MONOBUTYRIN_MG344 -3 LIPASE_TRIBUTYRIN_MG310 -3 LIPASE_TRIBUTYRIN_MG327 -3 LIPASE_TRIBUTYRIN_MG344 -3 LIPASE_TRILAURIN_MG310 -3 LIPASE_TRILAURIN_MG327 -3 LIPASE_TRILAURIN_MG344 -3 LIPASE_TRIMYRISTIN_MG310 -3 LIPASE_TRIMYRISTIN_MG327 -3 LIPASE_TRIMYRISTIN_MG344 -3 LIPASE_TRIOLEIN_MG310 -3 LIPASE_TRIOLEIN_MG327 -3 LIPASE_TRIOLEIN_MG344 -3 LIPASE_TRIPALMITIN_MG310 -3 LIPASE_TRIPALMITIN_MG327 -3 LIPASE_TRIPALMITIN_MG344 -3 LIPASE_TRISTEARIN_MG310 -3 LIPASE_TRISTEARIN_MG327 -3 LIPASE_TRISTEARIN_MG344 -3 LIPASE_TRI_HDCEA_IN_MG310 -3 LIPASE_TRI_HDCEA_IN_MG327 -3 LIPASE_TRI_HDCEA_IN_MG344 -3 LIPASE_TRI_TTDCEA_IN_MG310 -3 LIPASE_TRI_TTDCEA_IN_MG327 -3 LIPASE_TRI_TTDCEA_IN_MG344 -1 LIPASE_TWEEN80_MG310 -1 LIPASE_TWEEN80_MG327 -1 LIPASE_TWEEN80_MG344 -1 TX_ComEC_dApdAp -1 TX_ComEC_dCpdCp -1 TX_ComEC_dGpdGp -1 TX_ComEC_dTpdTp -1 TX_H2O + 1 H2Oex = 0
 H2O2_e: -1 TX_H2O2 + 1 H2O2ex = 0
 HCYS_e: 1 HCYSex = 0
 HDCA_e: 3 LIPASE_TRIPALMITIN_MG310 + 3 LIPASE_TRIPALMITIN_MG327 + 3 LIPASE_TRIPALMITIN_MG344 -1 TX_HDCA + 1 HDCAex = 0
 HDCEA_e: 3 LIPASE_TRI_HDCEA_IN_MG310 + 3 LIPASE_TRI_HDCEA_IN_MG327 + 3 LIPASE_TRI_HDCEA_IN_MG344 -1 TX_HDCEA + 1 HDCEAex = 0
 HDECEG3P_e: 1 HDECEG3Pex = 0
 HDECG3P_e: 1 HDECG3Pex = 0
 HIS_e: 1 HISex = 0
 HisHis_e: -1 TX_OppB_His + 1 HisHisex = 0
 hm3dC_e: -1 TX_hm3dC + 1 hm3dCex = 0
 hm3dCMP_e: 1 hm3dCMPex = 0
 hm5dC_e: -1 TX_hm5dC + 1 hm5dCex = 0
 hm5dCMP_e: 1 hm5dCMPex = 0
 hm5dU_e: -1 TX_hm5dU + 1 hm5dUex = 0
 hm5dUMP_e: 1 hm5dUMPex = 0
 hm5hydantoin_dRibose_e: -1 TX_hm5hydantoin_dRibose + 1 hm5hydantoin_dRiboseex = 0
 hm5hydantoin_dRiboseMP_e: 1 hm5hydantoin_dRiboseMPex = 0
 ho5dC_e: -1 TX_ho5dC + 1 ho5dCex = 0
 ho5dCMP_e: 1 ho5dCMPex = 0
 ho5dU_e: -1 TX_ho5dU + 1 ho5dUex = 0
 ho5dUMP_e: 1 ho5dUMPex = 0
 ho5ho6dT_e: -1 TX_ho5ho6dT + 1 ho5ho6dTex = 0
 ho5ho6dTMP_e: 1 ho5ho6dTMPex = 0
 ho5ho6dU_e: -1 TX_ho5ho6dU + 1 ho5ho6dUex = 0
 ho5ho6dUMP_e: 1 ho5ho6dUMPex = 0
 ho5hydantoin_dRibose_e: -1 TX_ho5hydantoin_dRibose + 1 ho5hydantoin_dRiboseex = 0
 ho5hydantoin_dRiboseMP_e: 1 ho5hydantoin_dRiboseMPex = 0
 HPHPYR_e: 1 HPHPYRex = 0
 hydride_e: 1 hydrideex = 0
 hydroxyl_radical_e: 1 hydroxyl_radicalex = 0
 HYXN_e: 1 HYXNex = 0
 IDP_e: 1 IDPex = 0
 ILE_e: 1 ILEex = 0
 IleIle_e: -1 TX_OppB_Ile + 1 IleIleex = 0
 IMP_e: 1 IMPex = 0
 INOS_e: 1 INOSex = 0
 INS_e: 1 INSex = 0
 ITP_e: 1 ITPex = 0
 K_e: -1 TX_TRK + 1 Kex = 0
 k2C_e: -1 TX_k2C + 1 k2Cex = 0
 k2CMP_e: 1 k2CMPex = 0
 LAC_e: 1 LACex = 0
 LEU_e: 1 LEUex = 0
 LeuLeu_e: -1 TX_OppB_Leu + 1 LeuLeuex = 0
 LIPOATE_e: -1 TX_LIPOATE + 1 LIPOATEex = 0
 LIPOYLAMP_e: 1 LIPOYLAMPex = 0
 LIPOYLLYS_e: 1 LIPOYLLYSex = 0
 lomefloxacin_e: 1 TX_NormM_lomefloxacin + 1 lomefloxacinex = 0
 LYS_e: 1 HinT_AMP_Lysine_MG132 + 1 HinT_GMP_Lysine_MG132 + 1 LYSex = 0
 LysLys_e: -1 TX_OppB_Lys + 1 LysLysex = 0
 m1dA_e: -1 TX_m1dA + 1 m1dAex = 0
 m1dAMP_e: 1 m1dAMPex = 0
 m1dG_e: -1 TX_m1dG + 1 m1dGex = 0
 m1dGMP_e: 1 m1dGMPex = 0
 m1G_e: -1 TX_m1G + 1 m1Gex = 0
 m1GMP_e: 1 m1GMPex = 0
 m2G_e: -1 TX_m2G + 1 m2Gex = 0
 m2GMP_e: 1 m2GMPex = 0
 m3dA_e: -1 TX_m3dA + 1 m3dAex = 0
 m3dAMP_e: 1 m3dAMPex = 0
 m3dC_e: -1 TX_m3dC + 1 m3dCex = 0
 m3dCMP_e: 1 m3dCMPex = 0
 m3dG_e: -1 TX_m3dG + 1 m3dGex = 0
 m3dGMP_e: 1 m3dGMPex = 0
 m3dT_e: -1 TX_m3dT + 1 m3dTex = 0
 m3dTMP_e: 1 m3dTMPex = 0
 m4dC_e: -1 TX_m4dC + 1 m4dCex = 0
 m4dCMP_e: 1 m4dCMPex = 0
 m5C_e: -1 TX_m5C + 1 m5Cex = 0
 m5CMP_e: 1 m5CMPex = 0
 m5dC_e: -1 TX_m5dC + 1 m5dCex = 0
 m5dCMP_e: 1 m5dCMPex = 0
 m62A_e: -1 TX_m62A + 1 m62Aex = 0
 m62AMP_e: 1 m62AMPex = 0
 m6A_e: -1 TX_m6A + 1 m6Aex = 0
 m6AMP_e: 1 m6AMPex = 0
 m6dA_e: -1 TX_m6dA + 1 m6dAex = 0
 m6dAMP_e: 1 m6dAMPex = 0
 m7dA_e: -1 TX_m7dA + 1 m7dAex = 0
 m7dAMP_e: 1 m7dAMPex = 0
 m7dG_e: -1 TX_m7dG + 1 m7dGex = 0
 m7dGMP_e: 1 m7dGMPex = 0
 m7G_e: -1 TX_m7G + 1 m7Gex = 0
 m7GMP_e: 1 m7GMPex = 0
 m8dG_e: -1 TX_m8dG + 1 m8dGex = 0
 m8dGMP_e: 1 m8dGMPex = 0
 MAL_e: 1 MALex = 0
 MAN_e: 1 MANex = 0
 MAN1P_e: 1 MAN1Pex = 0
 MAN6P_e: -1 TX_UHPT3 + 1 MAN6Pex = 0
 MET_e: 1 METex = 0
 METHF_e: 1 METHFex = 0
 methylRed_e: 1 methylRedex = 0
 METHYL_OCDCEA_e: -1 LIPASE_METHYL_OCDCEA_MG310 -1 LIPASE_METHYL_OCDCEA_MG327 -1 LIPASE_METHYL_OCDCEA_MG344 + 1 METHYL_OCDCEAex = 0
 MetMet_e: -1 TX_OppB_Met + 1 MetMetex = 0
 METSOXRL_e: 1 METSOXRLex = 0
 METSOXSL_e: 1 METSOXSLex = 0
 METTHF_e: 1 METTHFex = 0
 MG_e: -1 TX_MGTA_MG + 1 MGex = 0
 MN_e: -1 TX_MGTA_MN + 1 MNex = 0
 MO_e: 1 MOex = 0
 MOBD_e: -1 TX_MGTA_MOBD + 1 MOBDex = 0
 MONOBUTYRIN_e: -1 LIPASE_MONOBUTYRIN_MG310 -1 LIPASE_MONOBUTYRIN_MG327 -1 LIPASE_MONOBUTYRIN_MG344 + 1 MONOBUTYRINex = 0
 MOR_e: 1 HinT_AMP_Mor_MG132 + 1 HinT_GMP_Mor_MG132 + 1 MORex = 0
 MTHF_e: 1 MTHFex = 0
 n2dA_e: -1 TX_n2dA + 1 n2dAex = 0
 n2dAMP_e: 1 n2dAMPex = 0
 NA_e: -1 TX_NAUP2 -1 TX_NormM_acriflavine -1 TX_NormM_benzalkonium -1 TX_NormM_berberine -1 TX_NormM_ciprofloxacin -1 TX_NormM_crystalViolet -1 TX_NormM_doxorubicin -1 TX_NormM_enoxacin -1 TX_NormM_ethidium -1 TX_NormM_lomefloxacin -1 TX_NormM_norfloxacin -1 TX_NormM_novobiocin -1 TX_NormM_proflavine -1 TX_NormM_rhodamine6G -1 TX_NormM_tetraphenylphosphonium + 1 NAex = 0
 NAAD_e: 1 NAADex = 0
 NAC_e: -1 TX_NAC + 1 NACex = 0
 NAD_e: 1 NADex = 0
 NADH_e: 1 NADHex = 0
 NADP_e: 1 NADPex = 0
 NADPH_e: 1 NADPHex = 0
 NAMN_e: 1 NAMNex = 0
 NH3_e: 1 HinT_AMP_NH2_MG132 + 1 NH3ex = 0
 NH4_e: -1 TX_NH4 + 1 NH4ex = 0
 NI_e: -1 TX_Cbi_Ni + 1 NIex = 0
 nisin_e: 1 TX_SalX_nisin + 1 TX_SunT_nisin + 1 nisinex = 0
 NMN_e: 1 NMNex = 0
 norfloxacin_e: 1 TX_NormM_norfloxacin + 1 norfloxacinex = 0
 novobiocin_e: 1 TX_NormM_novobiocin + 1 novobiocinex = 0
 NR_e: 1 NRex = 0
 O2_e: -1 TX_O2 + 1 O2ex = 0
 O2edC_e: -1 TX_O2edC + 1 O2edCex = 0
 O2edCMP_e: 1 O2edCMPex = 0
 O2edT_e: -1 TX_O2edT + 1 O2edTex = 0
 O2edTMP_e: 1 O2edTMPex = 0
 O2mdC_e: -1 TX_O2mdC + 1 O2mdCex = 0
 O2mdCMP_e: 1 O2mdCMPex = 0
 O2mdT_e: -1 TX_O2mdT + 1 O2mdTex = 0
 O2mdTMP_e: 1 O2mdTMPex = 0
 O4edT_e: -1 TX_O4edT + 1 O4edTex = 0
 O4edTMP_e: 1 O4edTMPex = 0
 O4mdT_e: -1 TX_O4mdT + 1 O4mdTex = 0
 O4mdTMP_e: 1 O4mdTMPex = 0
 O6edG_e: -1 TX_O6edG + 1 O6edGex = 0
 O6edGMP_e: 1 O6edGMPex = 0
 O6mdG_e: -1 TX_O6mdG + 1 O6mdGex = 0
 O6mdGMP_e: 1 O6mdGMPex = 0
 OA_e: -1 TX_OA + 1 OAex = 0
 OCDCA_e: 3 LIPASE_TRISTEARIN_MG310 + 3 LIPASE_TRISTEARIN_MG327 + 3 LIPASE_TRISTEARIN_MG344 -1 TX_OCDCA + 1 OCDCAex = 0
 OCDCEA_e: 1 LIPASE_METHYL_OCDCEA_MG310 + 1 LIPASE_METHYL_OCDCEA_MG327 + 1 LIPASE_METHYL_OCDCEA_MG344 + 3 LIPASE_TRIOLEIN_MG310 + 3 LIPASE_TRIOLEIN_MG327 + 3 LIPASE_TRIOLEIN_MG344 + 1 LIPASE_TWEEN80_MG310 + 1 LIPASE_TWEEN80_MG327 + 1 LIPASE_TWEEN80_MG344 -1 TX_OCDCEA + 1 OCDCEAex = 0
 ODEC11GP_e: 1 ODEC11GPex = 0
 ODECG3P_e: 1 ODECG3Pex = 0
 OH_e: 1 OHex = 0
 ORN_e: 1 ORNex = 0
 oxo2dA_e: -1 TX_oxo2dA + 1 oxo2dAex = 0
 oxo2dAMP_e: 1 oxo2dAMPex = 0
 oxo8dA_e: -1 TX_oxo8dA + 1 oxo8dAex = 0
 oxo8dAMP_e: 1 oxo8dAMPex = 0
 oxo8dG_e: -1 TX_oxo8dG + 1 oxo8dGex = 0
 oxo8dGMP_e: 1 oxo8dGMPex = 0
 oxoethyl7dG_e: -1 TX_oxoethyl7dG + 1 oxoethyl7dGex = 0
 oxoethyl7dGMP_e: 1 oxoethyl7dGMPex = 0
 PA120_e: 1 PA120ex = 0
 PA140_e: 1 PA140ex = 0
 PA141_e: 1 PA141ex = 0
 PA160_e: 1 PA160ex = 0
 PA161_e: 1 PA161ex = 0
 PA180_e: 1 PA180ex = 0
 PA181_e: 1 PA181ex = 0
 PAP_e: 1 PAPex = 0
 pCHOL_e: 1 pCHOLex = 0
 PDLA5P_e: 1 PDLA5Pex = 0
 PEP_e: 1 PEPex = 0
 PG160_e: 1 PG160ex = 0
 PG161_e: 1 PG161ex = 0
 PG181_e: 1 PG181ex = 0
 PGP160_e: 1 PGP160ex = 0
 PGP161_e: 1 PGP161ex = 0
 PGP181_e: 1 PGP181ex = 0
 PHE_e: -1 TX_AROP12 -1 TX_AROP22 + 1 PHEex = 0
 PhePhe_e: -1 TX_OppB_Phe + 1 PhePheex = 0
 PHPYR_e: 1 PHPYRex = 0
 PI_e: -1 TX_PIUP1 + 2 TX_UHPT1 + 2 TX_UHPT2 + 2 TX_UHPT3 + 1 PIex = 0
 PL5P_e: 1 PL5Pex = 0
 POE_SORBITAN_e: 1 LIPASE_TWEEN80_MG310 + 1 LIPASE_TWEEN80_MG327 + 1 LIPASE_TWEEN80_MG344 + 1 POE_SORBITANex = 0
 PPAA_e: 1 PPAAex = 0
 PPGPP_e: 1 PPGPPex = 0
 PPHG_e: 1 PPHGex = 0
 PPI_e: 1 PPIex = 0
 PPIX_e: 1 PPIXex = 0
 PRO_e: 1 PROex = 0
 proflavine_e: 1 TX_NormM_proflavine + 1 proflavineex = 0
 ProPro_e: -1 TX_OppB_Pro + 1 ProProex = 0
 PRPP_e: 1 PRPPex = 0
 pSER_e: 1 pSERex = 0
 PSIURI_e: -1 TX_PSIURI + 1 PSIURIex = 0
 PSIURIMP_e: 1 PSIURIMPex = 0
 pTHR_e: 1 pTHRex = 0
 PTRC_e: -1 TX_SPRDUP2 + 1 PTRCex = 0
 pTYR_e: 1 pTYRex = 0
 putdT_e: -1 TX_putdT + 1 putdTex = 0
 putdTMP_e: 1 putdTMPex = 0
 PYDX_e: -1 TXPYDX + 1 PYDXex = 0
 PYR_e: 1 PYRex = 0
 R1P_e: 1 R1Pex = 0
 R5P_e: 1 R5Pex = 0
 RBFLVRD_e: 1 RBFLVRDex = 0
 rhodamine6G_e: 1 TX_NormM_rhodamine6G + 1 rhodamine6Gex = 0
 RIB_e: -1 TX_RIBUP + 1 RIBex = 0
 RIBFLV_e: -1 TX_RIBFLV + 1 RIBFLVex = 0
 RL5P_e: 1 RL5Pex = 0
 s2U_e: -1 TX_s2U + 1 s2Uex = 0
 s2UMP_e: 1 s2UMPex = 0
 s4U_e: -1 TX_s4U + 1 s4Uex = 0
 s4UMP_e: 1 s4UMPex = 0
 S7P_e: 1 S7Pex = 0
 SALA_e: 1 SALAex = 0
 SBT_e: 1 SBTex = 0
 SBT6P_e: 1 SBT6Pex = 0
 SER_e: 1 SERex = 0
 SerSer_e: -1 TX_OppB_Ser + 1 SerSerex = 0
 SNGLYP_e: 1 SNGLYPex = 0
 SO2_e: 1 SO2ex = 0
 SO4_e: -1 TX_SO4 + 1 SO4ex = 0
 SPMD_e: -1 TX_SPMDUP1 + 1 SPMDex = 0
 STEROL_e: -1 TX_STEROL + 1 STEROLex = 0
 T3P1_e: 1 T3P1ex = 0
 T3P2_e: 1 T3P2ex = 0
 TAG16P_e: 1 TAG16Pex = 0
 TAG6P_e: 1 TAG6Pex = 0
 tBHP_e: 1 tBHPex = 0
 tBOH_e: 1 tBOHex = 0
 TDECEG3P_e: 1 TDECEG3Pex = 0
 TDECG3P_e: 1 TDECG3Pex = 0
 TDP_e: 1 TDPex = 0
 tetraphenylphosphonium_e: 1 TX_NormM_tetraphenylphosphonium + 1 tetraphenylphosphoniumex = 0
 THF_e: -1 TX_THF + 1 THFex = 0
 THM_e: -1 TXTHM + 1 THMex = 0
 THMP_e: 1 THMPex = 0
 THR_e: 1 THRex = 0
 ThrThr_e: -1 TX_OppB_Thr + 1 ThrThrex = 0
 THY_e: -1 TX_THY + 1 THYex = 0
 TPP_e: 1 TPPex = 0
 TRIBUTYRIN_e: -1 LIPASE_TRIBUTYRIN_MG310 -1 LIPASE_TRIBUTYRIN_MG327 -1 LIPASE_TRIBUTYRIN_MG344 + 1 TRIBUTYRINex = 0
 TRILAURIN_e: -1 LIPASE_TRILAURIN_MG310 -1 LIPASE_TRILAURIN_MG327 -1 LIPASE_TRILAURIN_MG344 + 1 TRILAURINex = 0
 TRIMYRISTIN_e: -1 LIPASE_TRIMYRISTIN_MG310 -1 LIPASE_TRIMYRISTIN_MG327 -1 LIPASE_TRIMYRISTIN_MG344 + 1 TRIMYRISTINex = 0
 TRIOLEIN_e: -1 LIPASE_TRIOLEIN_MG310 -1 LIPASE_TRIOLEIN_MG327 -1 LIPASE_TRIOLEIN_MG344 + 1 TRIOLEINex = 0
 TRIPALMITIN_e: -1 LIPASE_TRIPALMITIN_MG310 -1 LIPASE_TRIPALMITIN_MG327 -1 LIPASE_TRIPALMITIN_MG344 + 1 TRIPALMITINex = 0
 TRISTEARIN_e: -1 LIPASE_TRISTEARIN_MG310 -1 LIPASE_TRISTEARIN_MG327 -1 LIPASE_TRISTEARIN_MG344 + 1 TRISTEARINex = 0
 TRI_HDCEA_IN_e: -1 LIPASE_TRI_HDCEA_IN_MG310 -1 LIPASE_TRI_HDCEA_IN_MG327 -1 LIPASE_TRI_HDCEA_IN_MG344 + 1 TRI_HDCEA_INex = 0
 TRI_TTDCEA_IN_e: -1 LIPASE_TRI_TTDCEA_IN_MG310 -1 LIPASE_TRI_TTDCEA_IN_MG327 -1 LIPASE_TRI_TTDCEA_IN_MG344 + 1 TRI_TTDCEA_INex = 0
 TRP_e: -1 TX_AROP11 -1 TX_AROP21 + 1 TRPex = 0
 TrpTrp_e: -1 TX_OppB_Trp + 1 TrpTrpex = 0
 TTDCA_e: 3 LIPASE_TRIMYRISTIN_MG310 + 3 LIPASE_TRIMYRISTIN_MG327 + 3 LIPASE_TRIMYRISTIN_MG344 -1 TX_TTDCA + 1 TTDCAex = 0
 TTDCEA_e: 3 LIPASE_TRI_TTDCEA_IN_MG310 + 3 LIPASE_TRI_TTDCEA_IN_MG327 + 3 LIPASE_TRI_TTDCEA_IN_MG344 -1 TX_TTDCEA + 1 TTDCEAex = 0
 TTP_e: 1 TTPex = 0
 TWEEN80_e: -1 LIPASE_TWEEN80_MG310 -1 LIPASE_TWEEN80_MG327 -1 LIPASE_TWEEN80_MG344 + 1 TWEEN80ex = 0
 TYR_e: -1 TX_AROP13 -1 TX_AROP23 + 1 TYRex = 0
 TyrTyr_e: -1 TX_OppB_Tyr + 1 TyrTyrex = 0
 U23CMP_e: 1 U23CMPex = 0
 U3MP_e: 1 U3MPex = 0
 UDP_e: 1 UDPex = 0
 UDPG_e: 1 UDPGex = 0
 UDPGAL_e: 1 UDPGALex = 0
 UDPGALF_e: 1 UDPGALFex = 0
 Um_e: -1 TX_Um + 1 Umex = 0
 UmMP_e: 1 UmMPex = 0
 UMP_e: 1 UMPex = 0
 URA_e: -1 TX_URA + 1 URAex = 0
 URI_e: 1 URIex = 0
 UTP_e: 1 UTPex = 0
 VAL_e: 1 VALex = 0
 ValVal_e: -1 TX_OppB_Val + 1 ValValex = 0
 vinyl_chloride_e: 1 vinyl_chlorideex = 0
 X5P_e: 1 X5Pex = 0
 XAN_e: 1 XANex = 0
 XTSN_e: 1 XTSNex = 0
 ZN_e: -1 TX_MGTA_ZN + 1 ZNex = 0
 MG_124_MONOMER_e: 1 MG_124_MONOMERex = 0
 MG_287_MONOMER_e: 1 MG_287_MONOMERex = 0
 MG_124_MONOMER_ox_e: 1 MG_124_MONOMER_oxex = 0
 MG_229_231_TETRAMER_e: 1 MG_229_231_TETRAMERex = 0
 MG_229_231_TETRAMER_ox_e: 1 MG_229_231_TETRAMER_oxex = 0
 MG_287_MONOMER_ACP_e: 1 MG_287_MONOMER_ACPex = 0
 MG_287_MONOMER_ddcaACP_e: 1 MG_287_MONOMER_ddcaACPex = 0
 MG_287_MONOMER_hdeACP_e: 1 MG_287_MONOMER_hdeACPex = 0
 MG_287_MONOMER_myrsACP_e: 1 MG_287_MONOMER_myrsACPex = 0
 MG_287_MONOMER_ocdcaACP_e: 1 MG_287_MONOMER_ocdcaACPex = 0
 MG_287_MONOMER_octeACP_e: 1 MG_287_MONOMER_octeACPex = 0
 MG_287_MONOMER_palmACP_e: 1 MG_287_MONOMER_palmACPex = 0
 MG_287_MONOMER_tdeACP_e: 1 MG_287_MONOMER_tdeACPex = 0
 MG_427_DIMER_e: 1 MG_427_DIMERex = 0
 MG_427_DIMER_ox_e: 1 MG_427_DIMER_oxex = 0
 MG_454_DIMER_e: 1 MG_454_DIMERex = 0
 MG_454_DIMER_ox_e: 1 MG_454_DIMER_oxex = 0
 CDPDG120_m: 1 CdsA1 = 0
 CDPDG140_m: 1 CdsA2 = 0
 CDPDG141_m: 1 CdsA3 = 0
 CDPDG160_m: 1 CdsA4 -1 PgsA1 = 0
 CDPDG161_m: 1 CdsA5 -1 PgsA2 = 0
 CDPDG180_m: 1 CdsA6 = 0
 CDPDG181_m: 1 CdsA7 -1 PgsA3 = 0
 Ceramide16_m: -1 glycosyltransferase_Ceramide16_GalCer16 -1 glycosyltransferase_Ceramide16_GlcCer16 = 0
 CL160_m: 1 Cls1 -2.710410e-03 Growth = 0
 CL161_m: 1 Cls2 -2.092557e-03 Growth = 0
 CL181_m: 1 Cls3 -1.404570e-03 Growth = 0
 DAG160_m: -1 glycosyltransferase_DAG160_GalDAG160 -1 glycosyltransferase_DAG160_GlcDAG160 + 1 glycosyltransferase_PG160_DAG160 -2.154425e-03 Growth = 0
 DAG161_m: -1 glycosyltransferase_DAG161_GalDAG161 -1 glycosyltransferase_DAG161_GlcDAG161 + 1 glycosyltransferase_PG161_DAG161 -1.659300e-03 Growth = 0
 DAG181_m: -1 glycosyltransferase_DAG181_GalDAG181 -1 glycosyltransferase_DAG181_GlcDAG181 + 1 glycosyltransferase_PG181_DAG181 -1.097734e-03 Growth = 0
 GalCer16_m: 1 glycosyltransferase_Ceramide16_GalCer16 -1 glycosyltransferase_GalCer16_GalGalCer16 -1 glycosyltransferase_GalCer16_GlcGalCer16 = 0
 GalDAG160_m: 1 glycosyltransferase_DAG160_GalDAG160 -1 glycosyltransferase_GalDAG160_GalGalDAG160 -1 glycosyltransferase_GalDAG160_GlcGalDAG160 = 0
 GalDAG161_m: 1 glycosyltransferase_DAG161_GalDAG161 -1 glycosyltransferase_GalDAG161_GalGalDAG161 -1 glycosyltransferase_GalDAG161_GlcGalDAG161 = 0
 GalDAG181_m: 1 glycosyltransferase_DAG181_GalDAG181 -1 glycosyltransferase_GalDAG181_GalGalDAG181 -1 glycosyltransferase_GalDAG181_GlcGalDAG181 = 0
 GalGalCer16_m: 1 glycosyltransferase_GalCer16_GalGalCer16 = 0
 GalGalDAG160_m: 1 glycosyltransferase_GalDAG160_GalGalDAG160 -1 glycosyltransferase_GalGalDAG160_GalGalGalDAG160 -1 glycosyltransferase_GalGalDAG160_GlcGalGalDAG160 = 0
 GalGalDAG161_m: 1 glycosyltransferase_GalDAG161_GalGalDAG161 -1 glycosyltransferase_GalGalDAG161_GalGalGalDAG161 -1 glycosyltransferase_GalGalDAG161_GlcGalGalDAG161 = 0
 GalGalDAG181_m: 1 glycosyltransferase_GalDAG181_GalGalDAG181 -1 glycosyltransferase_GalGalDAG181_GalGalGalDAG181 -1 glycosyltransferase_GalGalDAG181_GlcGalGalDAG181 = 0
 GalGalGalDAG160_m: 1 glycosyltransferase_GalGalDAG160_GalGalGalDAG160 = 0
 GalGalGalDAG161_m: 1 glycosyltransferase_GalGalDAG161_GalGalGalDAG161 = 0
 GalGalGalDAG181_m: 1 glycosyltransferase_GalGalDAG181_GalGalGalDAG181 = 0
 GalGalGlcCer16_m: 1 glycosyltransferase_GalGlcCer16_GalGalGlcCer16 = 0
 GalGalGlcDAG160_m: 1 glycosyltransferase_GalGlcDAG160_GalGalGlcDAG160 = 0
 GalGalGlcDAG161_m: 1 glycosyltransferase_GalGlcDAG161_GalGalGlcDAG161 = 0
 GalGalGlcDAG181_m: 1 glycosyltransferase_GalGlcDAG181_GalGalGlcDAG181 = 0
 GalGlcCer16_m: -1 glycosyltransferase_GalGlcCer16_GalGalGlcCer16 + 1 glycosyltransferase_GlcCer16_GalGlcCer16 = 0
 GalGlcDAG160_m: -1 glycosyltransferase_GalGlcDAG160_GalGalGlcDAG160 -1 glycosyltransferase_GalGlcDAG160_GlcGalGlcDAG160 + 1 glycosyltransferase_GlcDAG160_GalGlcDAG160 = 0
 GalGlcDAG161_m: -1 glycosyltransferase_GalGlcDAG161_GalGalGlcDAG161 -1 glycosyltransferase_GalGlcDAG161_GlcGalGlcDAG161 + 1 glycosyltransferase_GlcDAG161_GalGlcDAG161 = 0
 GalGlcDAG181_m: -1 glycosyltransferase_GalGlcDAG181_GalGalGlcDAG181 -1 glycosyltransferase_GalGlcDAG181_GlcGalGlcDAG181 + 1 glycosyltransferase_GlcDAG181_GalGlcDAG181 = 0
 GlcCer16_m: 1 glycosyltransferase_Ceramide16_GlcCer16 -1 glycosyltransferase_GlcCer16_GalGlcCer16 = 0
 GlcDAG160_m: 1 glycosyltransferase_DAG160_GlcDAG160 -1 glycosyltransferase_GlcDAG160_GalGlcDAG160 -1 glycosyltransferase_GlcDAG160_GlcGlcDAG160 = 0
 GlcDAG161_m: 1 glycosyltransferase_DAG161_GlcDAG161 -1 glycosyltransferase_GlcDAG161_GalGlcDAG161 -1 glycosyltransferase_GlcDAG161_GlcGlcDAG161 = 0
 GlcDAG181_m: 1 glycosyltransferase_DAG181_GlcDAG181 -1 glycosyltransferase_GlcDAG181_GalGlcDAG181 -1 glycosyltransferase_GlcDAG181_GlcGlcDAG181 = 0
 GlcGalCer16_m: 1 glycosyltransferase_GalCer16_GlcGalCer16 = 0
 GlcGalDAG160_m: 1 glycosyltransferase_GalDAG160_GlcGalDAG160 = 0
 GlcGalDAG161_m: 1 glycosyltransferase_GalDAG161_GlcGalDAG161 = 0
 GlcGalDAG181_m: 1 glycosyltransferase_GalDAG181_GlcGalDAG181 = 0
 GlcGalGalDAG160_m: 1 glycosyltransferase_GalGalDAG160_GlcGalGalDAG160 = 0
 GlcGalGalDAG161_m: 1 glycosyltransferase_GalGalDAG161_GlcGalGalDAG161 = 0
 GlcGalGalDAG181_m: 1 glycosyltransferase_GalGalDAG181_GlcGalGalDAG181 = 0
 GlcGalGlcDAG160_m: 1 glycosyltransferase_GalGlcDAG160_GlcGalGlcDAG160 = 0
 GlcGalGlcDAG161_m: 1 glycosyltransferase_GalGlcDAG161_GlcGalGlcDAG161 = 0
 GlcGalGlcDAG181_m: 1 glycosyltransferase_GalGlcDAG181_GlcGalGlcDAG181 = 0
 GlcGlcDAG160_m: 1 glycosyltransferase_GlcDAG160_GlcGlcDAG160 -1 glycosyltransferase_GlcGlcDAG160_GlcGlcGlcDAG160 = 0
 GlcGlcDAG161_m: 1 glycosyltransferase_GlcDAG161_GlcGlcDAG161 -1 glycosyltransferase_GlcGlcDAG161_GlcGlcGlcDAG161 = 0
 GlcGlcDAG181_m: 1 glycosyltransferase_GlcDAG181_GlcGlcDAG181 -1 glycosyltransferase_GlcGlcDAG181_GlcGlcGlcDAG181 = 0
 GlcGlcGlcDAG160_m: 1 glycosyltransferase_GlcGlcDAG160_GlcGlcGlcDAG160 = 0
 GlcGlcGlcDAG161_m: 1 glycosyltransferase_GlcGlcDAG161_GlcGlcGlcDAG161 = 0
 GlcGlcGlcDAG181_m: 1 glycosyltransferase_GlcGlcDAG181_GlcGlcGlcDAG181 = 0
 PA120_m: -1 CdsA1 + 1 PlsC1 = 0
 PA140_m: -1 CdsA2 + 1 PlsC2 = 0
 PA141_m: -1 CdsA3 + 1 PlsC3 = 0
 PA160_m: -1 CdsA4 + 1 PlsC4 = 0
 PA161_m: -1 CdsA5 + 1 PlsC5 = 0
 PA180_m: -1 CdsA6 + 1 PlsC6 = 0
 PA181_m: -1 CdsA7 + 1 PlsC7 = 0
 PG160_m: -2 Cls1 -1 glycosyltransferase_PG160_DAG160 + 1 PgpA1 -2.634973e-04 Growth = 0
 PG161_m: -2 Cls2 -1 glycosyltransferase_PG161_DAG161 + 1 PgpA2 = 0
 PG181_m: -2 Cls3 -1 glycosyltransferase_PG181_DAG181 + 1 PgpA3 = 0
 PGP160_m: -1 PgpA1 + 1 PgsA1 -1.832168e-02 Growth = 0
 PGP161_m: -1 PgpA2 + 1 PgsA2 -1.413206e-02 Growth = 0
 PGP181_m: -1 PgpA3 + 1 PgsA3 -9.601050e-03 Growth = 0
 STEROL_m: 1 TX_STEROL -2.938354e-02 Growth = 0
Bounds
\ mmol/gDCW/h
0 <= Aas1 <= 10000000000
0 <= Aas2 <= 10000000000
0 <= Aas3 <= 10000000000
0 <= Aas4 <= 10000000000
0 <= Aas5 <= 10000000000
0 <= Aas6 <= 10000000000
0 <= Aas7 <= 10000000000
-10000000000 <= AceE <= 1.914148e+01
-8.623237e+00 <= AckA <= 8.376859e+00
0 <= AcpS <= 5.140821e-02
0 <= Adi <= 10000000000
-10000000000 <= Adk1 <= 2.568424e+01
-10000000000 <= Adk2 <= 10000000000
-10000000000 <= Adk3 <= 10000000000
-1.802222e+01 <= Adk4 <= 10000000000
0 <= AHCYSs <= 10000000000
0 <= Apt <= 2.167548e-01
0 <= Apts_Ala <= 10000000000
0 <= Apts_Arg <= 10000000000
0 <= Apts_Asn <= 10000000000
0 <= Apts_Asp <= 10000000000
0 <= Apts_Cys <= 10000000000
0 <= Apts_Gln <= 10000000000
0 <= Apts_Glu <= 10000000000
0 <= Apts_Gly <= 10000000000
0 <= Apts_His <= 10000000000
0 <= Apts_Ile <= 10000000000
0 <= Apts_Leu <= 10000000000
0 <= Apts_Lys <= 10000000000
0 <= Apts_Met <= 10000000000
0 <= Apts_Phe <= 10000000000
0 <= Apts_Pro <= 10000000000
0 <= Apts_Ser <= 10000000000
0 <= Apts_Thr <= 10000000000
0 <= Apts_Trp <= 10000000000
0 <= Apts_Tyr <= 10000000000
0 <= Apts_Val <= 10000000000
-10000000000 <= ArcC <= 10000000000
-10000000000 <= ArgF <= 10000000000
-10000000000 <= AspC1 <= 10000000000
-10000000000 <= AspC2 <= 10000000000
-10000000000 <= AspC3 <= 10000000000
-6.045606e-02 <= AtpA <= 10000000000
0 <= AzoR <= 10000000000
0 <= BPNT <= 10000000000
0 <= Cdd1 <= 1.916789e+00
0 <= Cdd2 <= 1.916789e+00
-10000000000 <= CdsA1 <= 10000000000
-10000000000 <= CdsA2 <= 10000000000
-10000000000 <= CdsA3 <= 10000000000
-10000000000 <= CdsA4 <= 10000000000
-10000000000 <= CdsA5 <= 10000000000
-10000000000 <= CdsA6 <= 10000000000
-10000000000 <= CdsA7 <= 10000000000
0 <= CholK <= 10000000000
-10000000000 <= Cls1 <= 10000000000
-10000000000 <= Cls2 <= 10000000000
-10000000000 <= Cls3 <= 10000000000
-10000000000 <= CmkA1 <= 4.355513e+00
-10000000000 <= CmkA2 <= 4.140221e+00
-10000000000 <= CmkB3 <= 1.974435e+01
0 <= CoaE <= 3.916359e-02
0 <= CpdB1 <= 10000000000
0 <= CpdB2 <= 10000000000
0 <= CpdB3 <= 10000000000
0 <= CpdB4 <= 10000000000
-3.594763e+00 <= CpsG <= 10000000000
0 <= CsdA <= 7.671657e-01
0 <= Dak1 <= 3.123714e-01
-10000000000 <= DcdK <= 10000000000
0 <= Dck <= 10000000000
0 <= DeoC <= 5.813219e+00
-10000000000 <= DeoD1 <= 1.383478e-03
-10000000000 <= DeoD2 <= 4.464415e-01
-10000000000 <= DeoD3 <= 1.383478e-03
-1.383478e-03 <= DeoD4 <= 10000000000
-1.221459e-01 <= DeoD5 <= 10000000000
-1.383478e-03 <= DeoD6 <= 10000000000
-10000000000 <= DeoD7 <= 10000000000
-10000000000 <= DeoD8 <= 3.764826e-04
0 <= Dgk1 <= 6.263218e-02
0 <= DhaK <= 10000000000
-10000000000 <= DhfR1 <= 1.334352e-01
0 <= diadenylate_cyclase <= 10000000000
-10000000000 <= Eno <= 1.345633e+01
-10000000000 <= Fba <= 3.887461e+00
-10000000000 <= FMETHs <= 10000000000
-10000000000 <= FolD1 <= 4.259660e+00
-10000000000 <= FolD2 <= 8.253090e-01
0 <= Fre1 <= 1.376455e+01
0 <= Fre2 <= 10000000000
0 <= Fre3 <= 1.489810e+01
0 <= Fre4 <= 10000000000
0 <= Fre5 <= 1.511401e+00
0 <= Fre6 <= 2.240582e+00
0 <= FruK1 <= 1.394486e+00
0 <= FruK2 <= 6.089088e-02
-10000000000 <= FthC <= 6.452725e-01
0 <= FTHFS <= 1.594066e-01
-10000000000 <= GalE <= 7.711232e+00
-10000000000 <= GalU1 <= 10000000000
-1.085690e+01 <= GapA <= 1.197452e+01
0 <= Glf <= 10000000000
0 <= GLNHs <= 10000000000
0 <= GlpK <= 1.222161e+01
0 <= GlpQ1 <= 10000000000
0 <= GlpQ2 <= 10000000000
0 <= GlpQ3 <= 10000000000
0 <= GlpQ4 <= 10000000000
0 <= GlpQ5 <= 10000000000
0 <= GlyA1 <= 10000000000
0 <= GlyA2 <= 10000000000
0 <= GlyA3 <= 8.148904e-01
0 <= GlyA4 <= 1.018613e-01
0 <= glycosyltransferase_Ceramide16_GalCer16 <= 10000000000
0 <= glycosyltransferase_Ceramide16_GlcCer16 <= 10000000000
0 <= glycosyltransferase_DAG160_GalDAG160 <= 10000000000
0 <= glycosyltransferase_DAG160_GlcDAG160 <= 10000000000
0 <= glycosyltransferase_DAG161_GalDAG161 <= 10000000000
0 <= glycosyltransferase_DAG161_GlcDAG161 <= 10000000000
0 <= glycosyltransferase_DAG181_GalDAG181 <= 10000000000
0 <= glycosyltransferase_DAG181_GlcDAG181 <= 10000000000
0 <= glycosyltransferase_GalCer16_GalGalCer16 <= 10000000000
0 <= glycosyltransferase_GalCer16_GlcGalCer16 <= 10000000000
0 <= glycosyltransferase_GalDAG160_GalGalDAG160 <= 10000000000
0 <= glycosyltransferase_GalDAG160_GlcGalDAG160 <= 10000000000
0 <= glycosyltransferase_GalDAG161_GalGalDAG161 <= 10000000000
0 <= glycosyltransferase_GalDAG161_GlcGalDAG161 <= 10000000000
0 <= glycosyltransferase_GalDAG181_GalGalDAG181 <= 10000000000
0 <= glycosyltransferase_GalDAG181_GlcGalDAG181 <= 10000000000
0 <= glycosyltransferase_GalGalDAG160_GalGalGalDAG160 <= 10000000000
0 <= glycosyltransferase_GalGalDAG160_GlcGalGalDAG160 <= 10000000000
0 <= glycosyltransferase_GalGalDAG161_GalGalGalDAG161 <= 10000000000
0 <= glycosyltransferase_GalGalDAG161_GlcGalGalDAG161 <= 10000000000
0 <= glycosyltransferase_GalGalDAG181_GalGalGalDAG181 <= 10000000000
0 <= glycosyltransferase_GalGalDAG181_GlcGalGalDAG181 <= 10000000000
0 <= glycosyltransferase_GalGlcCer16_GalGalGlcCer16 <= 10000000000
0 <= glycosyltransferase_GalGlcDAG160_GalGalGlcDAG160 <= 10000000000
0 <= glycosyltransferase_GalGlcDAG160_GlcGalGlcDAG160 <= 10000000000
0 <= glycosyltransferase_GalGlcDAG161_GalGalGlcDAG161 <= 10000000000
0 <= glycosyltransferase_GalGlcDAG161_GlcGalGlcDAG161 <= 10000000000
0 <= glycosyltransferase_GalGlcDAG181_GalGalGlcDAG181 <= 10000000000
0 <= glycosyltransferase_GalGlcDAG181_GlcGalGlcDAG181 <= 10000000000
0 <= glycosyltransferase_GlcCer16_GalGlcCer16 <= 10000000000
0 <= glycosyltransferase_GlcDAG160_GalGlcDAG160 <= 10000000000
0 <= glycosyltransferase_GlcDAG160_GlcGlcDAG160 <= 10000000000
0 <= glycosyltransferase_GlcDAG161_GalGlcDAG161 <= 10000000000
0 <= glycosyltransferase_GlcDAG161_GlcGlcDAG161 <= 10000000000
0 <= glycosyltransferase_GlcDAG181_GalGlcDAG181 <= 10000000000
0 <= glycosyltransferase_GlcDAG181_GlcGlcDAG181 <= 10000000000
0 <= glycosyltransferase_GlcGlcDAG160_GlcGlcGlcDAG160 <= 10000000000
0 <= glycosyltransferase_GlcGlcDAG161_GlcGlcGlcDAG161 <= 10000000000
0 <= glycosyltransferase_GlcGlcDAG181_GlcGlcGlcDAG181 <= 10000000000
0 <= glycosyltransferase_PG160_DAG160 <= 10000000000
0 <= glycosyltransferase_PG161_DAG161 <= 10000000000
0 <= glycosyltransferase_PG181_DAG181 <= 10000000000
-10000000000 <= Gmk1 <= 2.701077e+01
-10000000000 <= Gmk2 <= 3.194822e+00
-10000000000 <= Gmk3 <= 3.194822e+00
-10000000000 <= Gmk4 <= 1.879307e-01
-10000000000 <= GpmA <= 1.157456e+02
-1.446500e-01 <= GpsA <= 10000000000
0 <= HCYSMT <= 10000000000
0 <= HemG <= 2.097877e+00
-10000000000 <= HinT_AMP_Lysine_MG132 <= 10000000000
-10000000000 <= HinT_AMP_Mor_MG132 <= 10000000000
-10000000000 <= HinT_AMP_NH2_MG132 <= 10000000000
-10000000000 <= HinT_GMP_Lysine_MG132 <= 10000000000
-10000000000 <= HinT_GMP_Mor_MG132 <= 10000000000
0 <= Hpt1 <= 4.659126e-01
0 <= Hpt2 <= 2.391099e-01
0 <= hydroxyl_radical_generation_1 <= 10000000000
0 <= hydroxyl_radical_generation_2 <= 10000000000
-10000000000 <= LdhA <= 1.831618e+02
-10000000000 <= LIPASE_DIBUTYRIN_MG310 <= 10000000000
-10000000000 <= LIPASE_DIBUTYRIN_MG327 <= 10000000000
-10000000000 <= LIPASE_DIBUTYRIN_MG344 <= 10000000000
-10000000000 <= LIPASE_METHYL_OCDCEA_MG310 <= 10000000000
-10000000000 <= LIPASE_METHYL_OCDCEA_MG327 <= 10000000000
-10000000000 <= LIPASE_METHYL_OCDCEA_MG344 <= 10000000000
-10000000000 <= LIPASE_MONOBUTYRIN_MG310 <= 10000000000
-10000000000 <= LIPASE_MONOBUTYRIN_MG327 <= 10000000000
-10000000000 <= LIPASE_MONOBUTYRIN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRIBUTYRIN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRIBUTYRIN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRIBUTYRIN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRILAURIN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRILAURIN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRILAURIN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRIMYRISTIN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRIMYRISTIN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRIMYRISTIN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRIOLEIN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRIOLEIN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRIOLEIN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRIPALMITIN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRIPALMITIN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRIPALMITIN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRISTEARIN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRISTEARIN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRISTEARIN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRI_HDCEA_IN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRI_HDCEA_IN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRI_HDCEA_IN_MG344 <= 10000000000
-10000000000 <= LIPASE_TRI_TTDCEA_IN_MG310 <= 10000000000
-10000000000 <= LIPASE_TRI_TTDCEA_IN_MG327 <= 10000000000
-10000000000 <= LIPASE_TRI_TTDCEA_IN_MG344 <= 10000000000
-10000000000 <= LIPASE_TWEEN80_MG310 <= 10000000000
-10000000000 <= LIPASE_TWEEN80_MG327 <= 10000000000
-10000000000 <= LIPASE_TWEEN80_MG344 <= 10000000000
0 <= LIPOYLLYSHs <= 10000000000
0 <= LplA <= 10000000000
-2.655846e+00 <= Mdh <= 1.007390e+00
0 <= MetF <= 10000000000
0 <= METHFHYD <= 10000000000
0 <= MetK <= 2.053136e-03
0 <= MetS <= 10000000000
0 <= MG385r1 <= 10000000000
0 <= MG385r2 <= 10000000000
0 <= MG385r3 <= 10000000000
0 <= MG385r4 <= 10000000000
0 <= MG385r5 <= 10000000000
0 <= MOrx <= 10000000000
-10000000000 <= MsrA <= 2.396651e-01
-10000000000 <= MsrB <= 4.919766e-02
0 <= NADD1 <= 10000000000
0 <= NADD2 <= 10000000000
0 <= NadE <= 1.402453e-01
-10000000000 <= NadK <= 1.377474e-01
0 <= NADR1 <= 2.313370e-01
0 <= NADR2 <= 1.034077e+00
-10000000000 <= NH3eq <= 10000000000
0 <= NOX <= 10000000000
0 <= NrdEF_ADP <= 10000000000
0 <= NrdEF_CDP <= 10000000000
0 <= NrdEF_GDP <= 10000000000
0 <= NrdEF_maintenance <= 10000000000
0 <= NrdEF_UDP <= 10000000000
0 <= Ohr_CHP <= 10000000000
0 <= Ohr_H2O2 <= 10000000000
0 <= Ohr_rx <= 10000000000
0 <= Ohr_tBHP <= 10000000000
0 <= OsmC_CHP <= 10000000000
0 <= OsmC_H2O2 <= 10000000000
0 <= OsmC_rx <= 10000000000
0 <= OsmC_tBHP <= 10000000000
-10000000000 <= PcbX <= 10000000000
-10000000000 <= Pdp1 <= 6.219393e-04
-10000000000 <= Pdp2 <= 2.379481e-04
-10000000000 <= Pdp3 <= 10000000000
0 <= PfkA1 <= 2.917641e+00
0 <= PfkA2 <= 2.917641e+00
0 <= PfkA3 <= 10000000000
0 <= PfkA4 <= 2.917641e+00
-10000000000 <= PfkA5 <= 1.144173e+00
-10000000000 <= Pgi <= 2.631544e+01
-1.176738e-01 <= Pgk1 <= 7.844921e-01
-1.098289e+00 <= Pgk2 <= 3.530214e+00
-3.146182e-03 <= Pgk3 <= 8.856504e+01
-1.019840e+00 <= Pgk4 <= 3.216418e+00
-10000000000 <= Pgk5 <= 10000000000
-10000000000 <= Pgm <= 4.790572e+00
0 <= PgpA1 <= 10000000000
0 <= PgpA2 <= 10000000000
0 <= PgpA3 <= 10000000000
-10000000000 <= PgsA1 <= 10000000000
-10000000000 <= PgsA2 <= 10000000000
-10000000000 <= PgsA3 <= 10000000000
0 <= PhnW <= 10000000000
0 <= PhnX <= 10000000000
0 <= PlsC1 <= 10000000000
0 <= PlsC2 <= 10000000000
0 <= PlsC3 <= 10000000000
0 <= PlsC4 <= 10000000000
0 <= PlsC5 <= 10000000000
0 <= PlsC6 <= 10000000000
0 <= PlsC7 <= 10000000000
0 <= PlsX1 <= 10000000000
0 <= PlsX2 <= 10000000000
0 <= PlsX3 <= 10000000000
0 <= PlsX4 <= 10000000000
0 <= PlsX5 <= 10000000000
0 <= PlsX6 <= 10000000000
0 <= PlsX7 <= 10000000000
0 <= PncB <= 6.396392e-02
0 <= Ppa <= 10000000000
0 <= Ppc <= 10000000000
-10000000000 <= PrsA <= 1.540089e+01
0 <= pSERpp <= 10000000000
-10000000000 <= Pta <= 5.749631e+01
0 <= pTHRpp <= 10000000000
0 <= pTYRpp <= 10000000000
0 <= PYDXK <= 10000000000
-10000000000 <= Pyk_ADP <= 1.563229e+01
-10000000000 <= Pyk_CDP <= 3.748560e+01
-10000000000 <= Pyk_DADP <= 10000000000
-10000000000 <= Pyk_DCDP <= 10000000000
-10000000000 <= Pyk_DGDP <= 10000000000
-10000000000 <= Pyk_DTDP <= 10000000000
-10000000000 <= Pyk_DUDP <= 10000000000
-10000000000 <= Pyk_GDP <= 3.566715e+01
-10000000000 <= Pyk_IDP <= 10000000000
-10000000000 <= Pyk_TDP <= 10000000000
-10000000000 <= Pyk_UDP <= 10000000000
0 <= RelA1 <= 10000000000
0 <= RelA2 <= 10000000000
0 <= RibF1 <= 1.132914e-02
0 <= RibF2 <= 1.132914e-02
-10000000000 <= Rpe <= 10000000000
-3.726221e+01 <= RpiB <= 10000000000
0 <= SNGLYPOr <= 10000000000
-10000000000 <= TalA <= 10000000000
0 <= TatD_dAMP <= 10000000000
0 <= TatD_dCMP <= 10000000000
0 <= TatD_dGMP <= 10000000000
0 <= TatD_dTMP <= 10000000000
0 <= Tdk1 <= 1.726734e-01
0 <= Tdk2 <= 1.726734e-01
-10000000000 <= ThiL <= 4.924975e-05
0 <= ThyA <= 7.226966e-03
-10000000000 <= Tkt1 <= 4.589589e+00
-10000000000 <= Tkt2 <= 4.589589e+00
0 <= TMDPK <= 10000000000
-10000000000 <= Tmk <= 1.413823e-02
-10000000000 <= TpiA <= 4.030404e+01
-10000000000 <= TrxB <= 7.711232e-01
0 <= TXPYDX <= 10000000000
0 <= TXTHM <= 10000000000
-10000000000 <= TX_ACAL <= 10000000000
-10000000000 <= TX_ACUP <= 10000000000
-10000000000 <= TX_AD <= 10000000000
-10000000000 <= TX_AKG <= 10000000000
-10000000000 <= TX_AROP11 <= 10000000000
-10000000000 <= TX_AROP12 <= 10000000000
-10000000000 <= TX_AROP13 <= 10000000000
-10000000000 <= TX_AROP21 <= 10000000000
-10000000000 <= TX_AROP22 <= 10000000000
-10000000000 <= TX_AROP23 <= 10000000000
-10000000000 <= TX_aza6dT <= 10000000000
-10000000000 <= TX_aza8dA <= 10000000000
-10000000000 <= TX_c5dU <= 10000000000
-10000000000 <= TX_CAP <= 10000000000
-10000000000 <= TX_carbamoylmethyl6dA <= 10000000000
0 <= TX_Cbi_Co <= 10000000000
0 <= TX_Cbi_Ni <= 10000000000
0 <= TX_CL <= 10000000000
-10000000000 <= TX_cmnm5s2U <= 10000000000
-10000000000 <= TX_CO2 <= 10000000000
-10000000000 <= TX_COA <= 10000000000
0 <= TX_ComEC_dApdAp <= 10000000000
0 <= TX_ComEC_dCpdCp <= 10000000000
0 <= TX_ComEC_dGpdGp <= 10000000000
0 <= TX_ComEC_dTpdTp <= 10000000000
-10000000000 <= TX_cyclobutane_dCdC <= 10000000000
-10000000000 <= TX_cyclobutane_dCdT <= 10000000000
-10000000000 <= TX_cyclobutane_dTdC <= 10000000000
-10000000000 <= TX_cyclobutane_dTdT <= 10000000000
-10000000000 <= TX_cyclodA <= 10000000000
-10000000000 <= TX_cyclodG <= 10000000000
-10000000000 <= TX_CYTD <= 10000000000
0 <= TX_DA <= 10000000000
-10000000000 <= TX_dC64dC <= 10000000000
-10000000000 <= TX_dC64dC_dewar <= 10000000000
-10000000000 <= TX_dC64dT <= 10000000000
-10000000000 <= TX_dC64dT_dewar <= 10000000000
-10000000000 <= TX_dC_GLYC <= 10000000000
-10000000000 <= TX_DDCA <= 10000000000
-10000000000 <= TX_DHdT <= 10000000000
-10000000000 <= TX_DHdU <= 10000000000
-10000000000 <= TX_dhpdU <= 10000000000
0 <= TX_diacylglycerolCys <= 10000000000
-10000000000 <= TX_dINS <= 10000000000
-10000000000 <= TX_dT64dC <= 10000000000
-10000000000 <= TX_dT64dC_dewar <= 10000000000
-10000000000 <= TX_dT64dT <= 10000000000
-10000000000 <= TX_dT64dT_dewar <= 10000000000
-10000000000 <= TX_dT_GLYC <= 10000000000
-10000000000 <= TX_dXAN <= 10000000000
-10000000000 <= TX_e1dA <= 10000000000
-10000000000 <= TX_e1dG <= 10000000000
-10000000000 <= TX_e3dA <= 10000000000
-10000000000 <= TX_e3dC <= 10000000000
-10000000000 <= TX_e3dG <= 10000000000
-10000000000 <= TX_e3dT <= 10000000000
-10000000000 <= TX_e4dC <= 10000000000
-10000000000 <= TX_e5dC <= 10000000000
-10000000000 <= TX_e6dA <= 10000000000
-10000000000 <= TX_e7dA <= 10000000000
-10000000000 <= TX_e7dG <= 10000000000
-10000000000 <= TX_e8dG <= 10000000000
-10000000000 <= TX_epsilon1N2dG <= 10000000000
-10000000000 <= TX_epsilon1N6dA <= 10000000000
-10000000000 <= TX_epsilon3N4dC <= 10000000000
-10000000000 <= TX_epsilonN23dG <= 10000000000
-10000000000 <= TX_f5dU <= 10000000000
-10000000000 <= TX_FALD <= 10000000000
-10000000000 <= TX_FAPydA <= 10000000000
-10000000000 <= TX_FAPydG <= 10000000000
-10000000000 <= TX_FOR <= 10000000000
0 <= TX_FRUPTS <= 10000000000
0 <= TX_GL <= 10000000000
0 <= TX_GL3P <= 10000000000
-10000000000 <= TX_glchm5dC <= 10000000000
-10000000000 <= TX_glchm5dU <= 10000000000
0 <= TX_GLCPTS <= 10000000000
0 <= TX_GLNHPQ <= 10000000000
-10000000000 <= TX_gludT <= 10000000000
0 <= TX_GLUP <= 10000000000
-10000000000 <= TX_Gm <= 10000000000
0 <= TX_GN <= 10000000000
-10000000000 <= TX_H2O <= 10000000000
-10000000000 <= TX_H2O2 <= 10000000000
-10000000000 <= TX_HDCA <= 10000000000
-10000000000 <= TX_HDCEA <= 10000000000
-10000000000 <= TX_hm3dC <= 10000000000
-10000000000 <= TX_hm5dC <= 10000000000
-10000000000 <= TX_hm5dU <= 10000000000
-10000000000 <= TX_hm5hydantoin_dRibose <= 10000000000
-10000000000 <= TX_ho5dC <= 10000000000
-10000000000 <= TX_ho5dU <= 10000000000
-10000000000 <= TX_ho5ho6dT <= 10000000000
-10000000000 <= TX_ho5ho6dU <= 10000000000
-10000000000 <= TX_ho5hydantoin_dRibose <= 10000000000
-10000000000 <= TX_k2C <= 10000000000
-10000000000 <= TX_LIPOATE <= 10000000000
-10000000000 <= TX_m1dA <= 10000000000
-10000000000 <= TX_m1dG <= 10000000000
-10000000000 <= TX_m1G <= 10000000000
-10000000000 <= TX_m2G <= 10000000000
-10000000000 <= TX_m3dA <= 10000000000
-10000000000 <= TX_m3dC <= 10000000000
-10000000000 <= TX_m3dG <= 10000000000
-10000000000 <= TX_m3dT <= 10000000000
-10000000000 <= TX_m4dC <= 10000000000
-10000000000 <= TX_m5C <= 10000000000
-10000000000 <= TX_m5dC <= 10000000000
-10000000000 <= TX_m62A <= 10000000000
-10000000000 <= TX_m6A <= 10000000000
-10000000000 <= TX_m6dA <= 10000000000
-10000000000 <= TX_m7dA <= 10000000000
-10000000000 <= TX_m7dG <= 10000000000
-10000000000 <= TX_m7G <= 10000000000
-10000000000 <= TX_m8dG <= 10000000000
0 <= TX_MdlB_doxorubicin <= 10000000000
0 <= TX_MGTA_CA <= 10000000000
0 <= TX_MGTA_CU <= 10000000000
0 <= TX_MGTA_FE2 <= 10000000000
0 <= TX_MGTA_FE3 <= 10000000000
0 <= TX_MGTA_MG <= 10000000000
0 <= TX_MGTA_MN <= 10000000000
0 <= TX_MGTA_MOBD <= 10000000000
0 <= TX_MGTA_ZN <= 10000000000
-10000000000 <= TX_n2dA <= 10000000000
0 <= TX_NAC <= 10000000000
-10000000000 <= TX_NAUP2 <= 10000000000
-10000000000 <= TX_NH4 <= 10000000000
0 <= TX_NormM_acriflavine <= 10000000000
0 <= TX_NormM_benzalkonium <= 10000000000
0 <= TX_NormM_berberine <= 10000000000
0 <= TX_NormM_ciprofloxacin <= 10000000000
0 <= TX_NormM_crystalViolet <= 10000000000
0 <= TX_NormM_doxorubicin <= 10000000000
0 <= TX_NormM_enoxacin <= 10000000000
0 <= TX_NormM_ethidium <= 10000000000
0 <= TX_NormM_lomefloxacin <= 10000000000
0 <= TX_NormM_norfloxacin <= 10000000000
0 <= TX_NormM_novobiocin <= 10000000000
0 <= TX_NormM_proflavine <= 10000000000
0 <= TX_NormM_rhodamine6G <= 10000000000
0 <= TX_NormM_tetraphenylphosphonium <= 10000000000
-10000000000 <= TX_O2 <= 10000000000
-10000000000 <= TX_O2edC <= 10000000000
-10000000000 <= TX_O2edT <= 10000000000
-10000000000 <= TX_O2mdC <= 10000000000
-10000000000 <= TX_O2mdT <= 10000000000
-10000000000 <= TX_O4edT <= 10000000000
-10000000000 <= TX_O4mdT <= 10000000000
-10000000000 <= TX_O6edG <= 10000000000
-10000000000 <= TX_O6mdG <= 10000000000
-10000000000 <= TX_OA <= 10000000000
-10000000000 <= TX_OCDCA <= 10000000000
-10000000000 <= TX_OCDCEA <= 10000000000
0 <= TX_OppB_Ala <= 10000000000
0 <= TX_OppB_Arg <= 10000000000
0 <= TX_OppB_Asn <= 10000000000
0 <= TX_OppB_Asp <= 10000000000
0 <= TX_OppB_Cys <= 10000000000
0 <= TX_OppB_Gln <= 10000000000
0 <= TX_OppB_Glu <= 10000000000
0 <= TX_OppB_Gly <= 10000000000
0 <= TX_OppB_His <= 10000000000
0 <= TX_OppB_Ile <= 10000000000
0 <= TX_OppB_Leu <= 10000000000
0 <= TX_OppB_Lys <= 10000000000
0 <= TX_OppB_Met <= 10000000000
0 <= TX_OppB_Phe <= 10000000000
0 <= TX_OppB_Pro <= 10000000000
0 <= TX_OppB_Ser <= 10000000000
0 <= TX_OppB_Thr <= 10000000000
0 <= TX_OppB_Trp <= 10000000000
0 <= TX_OppB_Tyr <= 10000000000
0 <= TX_OppB_Val <= 10000000000
-10000000000 <= TX_oxo2dA <= 10000000000
-10000000000 <= TX_oxo8dA <= 10000000000
-10000000000 <= TX_oxo8dG <= 10000000000
-10000000000 <= TX_oxoethyl7dG <= 10000000000
0 <= TX_PHNT <= 10000000000
0 <= TX_PIUP1 <= 10000000000
-10000000000 <= TX_PSIURI <= 10000000000
-10000000000 <= TX_putdT <= 10000000000
-10000000000 <= TX_RIBFLV <= 10000000000
0 <= TX_RIBUP <= 10000000000
-10000000000 <= TX_s2U <= 10000000000
-10000000000 <= TX_s4U <= 10000000000
0 <= TX_SalX_nisin <= 10000000000
0 <= TX_SO4 <= 10000000000
0 <= TX_SPMDUP1 <= 10000000000
0 <= TX_SPRDUP2 <= 10000000000
-10000000000 <= TX_STEROL <= 10000000000
0 <= TX_SunT_nisin <= 10000000000
-10000000000 <= TX_THF <= 10000000000
0 <= TX_THY <= 10000000000
-10000000000 <= TX_TRK <= 10000000000
-10000000000 <= TX_TTDCA <= 10000000000
-10000000000 <= TX_TTDCEA <= 10000000000
0 <= TX_UHPT1 <= 10000000000
0 <= TX_UHPT2 <= 10000000000
0 <= TX_UHPT3 <= 10000000000
-10000000000 <= TX_Um <= 10000000000
-10000000000 <= TX_URA <= 10000000000
0 <= Udk1 <= 6.158392e+00
0 <= Udk2 <= 6.158392e+00
0 <= Upp <= 1.862594e-03
0 <= UshA_AMP <= 10000000000
0 <= UshA_aza6dTMP <= 10000000000
0 <= UshA_aza8dAMP <= 10000000000
0 <= UshA_c5dUMP <= 10000000000
0 <= UshA_carbamoylmethyl6dAMP <= 10000000000
0 <= UshA_cmnm5s2UMP <= 10000000000
0 <= UshA_cyclobutane_dCMPdCMP <= 10000000000
0 <= UshA_cyclobutane_dCMPdTMP <= 10000000000
0 <= UshA_cyclobutane_dTMPdCMP <= 10000000000
0 <= UshA_cyclobutane_dTMPdTMP <= 10000000000
0 <= UshA_cyclodAMP <= 10000000000
0 <= UshA_cyclodGMP <= 10000000000
0 <= UshA_dCMP <= 10000000000
0 <= UshA_dCMP64dCMP <= 10000000000
0 <= UshA_dCMP64dCMP_dewar <= 10000000000
0 <= UshA_dCMP64dTMP <= 10000000000
0 <= UshA_dCMP64dTMP_dewar <= 10000000000
0 <= UshA_dCMP_GLYC <= 10000000000
0 <= UshA_DHdTMP <= 10000000000
0 <= UshA_DHdUMP <= 10000000000
0 <= UshA_dhpdUMP <= 10000000000
0 <= UshA_dINSMP <= 10000000000
0 <= UshA_dTMP <= 10000000000
0 <= UshA_dTMP64dCMP <= 10000000000
0 <= UshA_dTMP64dCMP_dewar <= 10000000000
0 <= UshA_dTMP64dTMP <= 10000000000
0 <= UshA_dTMP64dTMP_dewar <= 10000000000
0 <= UshA_dTMP_GLY <= 10000000000
0 <= UshA_dUMP <= 10000000000
0 <= UshA_dXANMP <= 10000000000
0 <= UshA_e1dAMP <= 10000000000
0 <= UshA_e1dGMP <= 10000000000
0 <= UshA_e3dAMP <= 10000000000
0 <= UshA_e3dCMP <= 10000000000
0 <= UshA_e3dGMP <= 10000000000
0 <= UshA_e3dTMP <= 10000000000
0 <= UshA_e4dCMP <= 10000000000
0 <= UshA_e5dCMP <= 10000000000
0 <= UshA_e6dAMP <= 10000000000
0 <= UshA_e7dAMP <= 10000000000
0 <= UshA_e7dGMP <= 10000000000
0 <= UshA_e8dGMP <= 10000000000
0 <= UshA_epsilon1N2dGMP <= 10000000000
0 <= UshA_epsilon1N6dAMP <= 10000000000
0 <= UshA_epsilon3N4dCMP <= 10000000000
0 <= UshA_epsilonN23dGMP <= 10000000000
0 <= UshA_f5dUMP <= 10000000000
0 <= UshA_FAPydAMP <= 10000000000
0 <= UshA_FAPydGMP <= 10000000000
0 <= UshA_glchm5dCMP <= 10000000000
0 <= UshA_glchm5dUMP <= 10000000000
0 <= UshA_gludTMP <= 10000000000
0 <= UshA_GmMP <= 10000000000
0 <= UshA_GMP <= 10000000000
0 <= UshA_hm3dCMP <= 10000000000
0 <= UshA_hm5dCMP <= 10000000000
0 <= UshA_hm5dUMP <= 10000000000
0 <= UshA_hm5hydantoin_dRiboseMP <= 10000000000
0 <= UshA_ho5dCMP <= 10000000000
0 <= UshA_ho5dUMP <= 10000000000
0 <= UshA_ho5ho6dTMP <= 10000000000
0 <= UshA_ho5ho6dUMP <= 10000000000
0 <= UshA_ho5hydantoin_dRiboseMP <= 10000000000
0 <= UshA_IMP <= 10000000000
0 <= UshA_k2CMP <= 10000000000
0 <= UshA_m1dAMP <= 10000000000
0 <= UshA_m1dGMP <= 10000000000
0 <= UshA_m1GMP <= 10000000000
0 <= UshA_m2GMP <= 10000000000
0 <= UshA_m3dAMP <= 10000000000
0 <= UshA_m3dCMP <= 10000000000
0 <= UshA_m3dGMP <= 10000000000
0 <= UshA_m3dTMP <= 10000000000
0 <= UshA_m4dCMP <= 10000000000
0 <= UshA_m5CMP <= 10000000000
0 <= UshA_m5dCMP <= 10000000000
0 <= UshA_m62AMP <= 10000000000
0 <= UshA_m6AMP <= 10000000000
0 <= UshA_m6dAMP <= 10000000000
0 <= UshA_m7dAMP <= 10000000000
0 <= UshA_m7dGMP <= 10000000000
0 <= UshA_m7GMP <= 10000000000
0 <= UshA_m8dGMP <= 10000000000
0 <= UshA_n2dAMP <= 10000000000
0 <= UshA_O2edCMP <= 10000000000
0 <= UshA_O2edTMP <= 10000000000
0 <= UshA_O2mdCMP <= 10000000000
0 <= UshA_O2mdTMP <= 10000000000
0 <= UshA_O4edTMP <= 10000000000
0 <= UshA_O4mdTMP <= 10000000000
0 <= UshA_O6edGMP <= 10000000000
0 <= UshA_O6mdGMP <= 10000000000
0 <= UshA_oxo2dAMP <= 10000000000
0 <= UshA_oxo8dAMP <= 10000000000
0 <= UshA_oxo8dGMP <= 10000000000
0 <= UshA_oxoethyl7dGMP <= 10000000000
0 <= UshA_PSIURIMP <= 10000000000
0 <= UshA_putdTMP <= 10000000000
0 <= UshA_s2UMP <= 10000000000
0 <= UshA_s4UMP <= 10000000000
0 <= UshA_UmMP <= 10000000000
0 <= vinyl_chloride_activation_1 <= 10000000000
0 <= vinyl_chloride_activation_2 <= 10000000000
0 <= YbiV_deoxyG6P <= 2.914846e-02
0 <= YbiV_F6P <= 1.995050e-02
0 <= YbiV_G6P <= 1.147046e-01
0 <= YbiV_GL2P <= 1.619359e-02
0 <= YbiV_MAN6P <= 2.590974e-02
0 <= YbiV_R5P <= 1.295487e-01
0 <= YbiV_SBT6P <= 2.539154e-02
0 <= YbiV_SNGLYP <= 2.914846e-02
0 <= YicP <= 10000000000
gamma_radiationex = 0
-12 <= A23CMPex <= 0
-12 <= A3MPex <= 0
-12 <= ACex <= 12
-12 <= ACALex <= 0
-12 <= ACCOAex <= 0
-12 <= acriflavineex <= 0
-12 <= ACTPex <= 0
-12 <= ADex <= 12
-12 <= ADNex <= 0
-12 <= ADPex <= 0
-12 <= AEPPex <= 12
-12 <= AHCYSex <= 0
-12 <= AKGex <= 0
-12 <= ALAex <= 0
-12 <= AlaAlaex <= 12
-12 <= AMETex <= 0
-12 <= aminodimethylanilineex <= 0
-12 <= AMPex <= 0
-12 <= AMP_Lysineex <= 0
-12 <= AMP_Morex <= 0
-12 <= AMP_NH2ex <= 0
-12 <= anthranilateex <= 12
-12 <= ARGex <= 0
-12 <= ArgArgex <= 12
-12 <= ASNex <= 0
-12 <= AsnAsnex <= 12
-12 <= ASPex <= 0
-12 <= AspAspex <= 12
-12 <= ATPex <= 0
-12 <= aza6dTex <= 0
-12 <= aza6dTMPex <= 0
-12 <= aza8dAex <= 0
-12 <= aza8dAMPex <= 0
-12 <= benzalkoniumex <= 0
-12 <= berberineex <= 0
-12 <= BUTex <= 0
-12 <= C23CMPex <= 0
-12 <= C3MPex <= 0
-12 <= c5dUex <= 0
-12 <= c5dUMPex <= 0
-20 <= CAex <= 20
-12 <= CAPex <= 0
-12 <= carbamoylmethyl6dAex <= 0
-12 <= carbamoylmethyl6dAMPex <= 0
-12 <= CDPex <= 0
-12 <= CDPDG120ex <= 0
-12 <= CDPDG140ex <= 0
-12 <= CDPDG141ex <= 0
-12 <= CDPDG160ex <= 0
-12 <= CDPDG161ex <= 0
-12 <= CDPDG180ex <= 0
-12 <= CDPDG181ex <= 0
-12 <= Ceramide16ex <= 0
-12 <= CH3OHex <= 0
-12 <= chloroacetaldehydeex <= 0
-12 <= chloroethylene_oxideex <= 0
-12 <= CHOLex <= 12
-12 <= CHPex <= 0
-12 <= ciprofloxacinex <= 0
-12 <= CITRex <= 0
-20 <= CLex <= 20
-12 <= CL160ex <= 0
-12 <= CL161ex <= 0
-12 <= CL181ex <= 0
-12 <= cmnm5s2Uex <= 0
-12 <= cmnm5s2UMPex <= 0
-12 <= CMPex <= 0
-20 <= COex <= 20
-12 <= CO2ex <= 12
-12 <= COAex <= 12
-12 <= crystalVioletex <= 0
-12 <= CTPex <= 0
-20 <= CUex <= 20
-12 <= cyclobutane_dCdCex <= 0
-12 <= cyclobutane_dCdTex <= 0
-12 <= cyclobutane_dCMPdCMPex <= 0
-12 <= cyclobutane_dCMPdTMPex <= 0
-12 <= cyclobutane_dTdCex <= 0
-12 <= cyclobutane_dTdTex <= 0
-12 <= cyclobutane_dTMPdCMPex <= 0
-12 <= cyclobutane_dTMPdTMPex <= 0
-12 <= cyclodAex <= 0
-12 <= cyclodAMPex <= 0
-12 <= cyclodGex <= 0
-12 <= cyclodGMPex <= 0
-12 <= CYSex <= 0
-12 <= CysCysex <= 12
-12 <= CYTDex <= 12
-12 <= c_di_AMPex <= 0
-12 <= DAex <= 12
-12 <= DADPex <= 0
-12 <= DAG160ex <= 0
-12 <= DAG161ex <= 0
-12 <= DAG181ex <= 0
-12 <= DALAex <= 0
-12 <= DAMPex <= 0
-12 <= dApdApex <= 0
-12 <= DATPex <= 0
-12 <= DCex <= 12
-12 <= dC64dCex <= 0
-12 <= dC64dC_dewarex <= 0
-12 <= dC64dTex <= 0
-12 <= dC64dT_dewarex <= 0
-12 <= DCDPex <= 0
-12 <= DCMPex <= 0
-12 <= dCMP64dCMPex <= 0
-12 <= dCMP64dCMP_dewarex <= 0
-12 <= dCMP64dTMPex <= 0
-12 <= dCMP64dTMP_dewarex <= 0
-12 <= dCMP_GLYCex <= 0
-12 <= dCpdCpex <= 0
-12 <= DCTPex <= 0
-12 <= dC_GLYCex <= 0
-12 <= DDCAex <= 12
-12 <= DDECG3Pex <= 0
-12 <= deoxyG6Pex <= 0
-12 <= deoxyGLCex <= 0
-12 <= DGex <= 12
-12 <= DGDPex <= 0
-12 <= DGMPex <= 0
-12 <= dGpdGpex <= 0
-12 <= DGTPex <= 0
-12 <= DHAex <= 0
-12 <= DHdTex <= 0
-12 <= DHdTMPex <= 0
-12 <= DHdUex <= 0
-12 <= DHdUMPex <= 0
-12 <= DHFex <= 12
-12 <= dhpdUex <= 0
-12 <= dhpdUMPex <= 0
-12 <= diacylglycerolCysex <= 0
-12 <= DIBUTYRINex <= 0
-12 <= dINSex <= 0
-12 <= dINSMPex <= 0
-12 <= DMBAex <= 0
-12 <= doxorubicinex <= 0
-12 <= DPCOAex <= 0
-12 <= DPGex <= 0
-12 <= DR1Pex <= 0
-12 <= DR5Pex <= 0
-12 <= DTex <= 12
-12 <= dT64dCex <= 0
-12 <= dT64dC_dewarex <= 0
-12 <= dT64dTex <= 0
-12 <= dT64dT_dewarex <= 0
-12 <= DTDPex <= 0
-12 <= DTMPex <= 0
-12 <= dTMP64dCMPex <= 0
-12 <= dTMP64dCMP_dewarex <= 0
-12 <= dTMP64dTMPex <= 0
-12 <= dTMP64dTMP_dewarex <= 0
-12 <= dTMP_GLYCex <= 0
-12 <= dTpdTpex <= 0
-12 <= DTTPex <= 0
-12 <= dT_GLYCex <= 0
-12 <= DUex <= 0
-12 <= DUDPex <= 0
-12 <= DUMPex <= 0
-12 <= DUTPex <= 0
-12 <= dXANex <= 0
-12 <= dXANMPex <= 0
-12 <= e1dAex <= 0
-12 <= e1dAMPex <= 0
-12 <= e1dGex <= 0
-12 <= e1dGMPex <= 0
-12 <= e3dAex <= 0
-12 <= e3dAMPex <= 0
-12 <= e3dCex <= 0
-12 <= e3dCMPex <= 0
-12 <= e3dGex <= 0
-12 <= e3dGMPex <= 0
-12 <= e3dTex <= 0
-12 <= e3dTMPex <= 0
-12 <= e4dCex <= 0
-12 <= e4dCMPex <= 0
-12 <= E4Pex <= 0
-12 <= e5dCex <= 0
-12 <= e5dCMPex <= 0
-12 <= e6dAex <= 0
-12 <= e6dAMPex <= 0
-12 <= e7dAex <= 0
-12 <= e7dAMPex <= 0
-12 <= e7dGex <= 0
-12 <= e7dGMPex <= 0
-12 <= e8dGex <= 0
-12 <= e8dGMPex <= 0
-12 <= enoxacinex <= 0
-12 <= epsilon1N2dGex <= 0
-12 <= epsilon1N2dGMPex <= 0
-12 <= epsilon1N6dAex <= 0
-12 <= epsilon1N6dAMPex <= 0
-12 <= epsilon3N4dCex <= 0
-12 <= epsilon3N4dCMPex <= 0
-12 <= epsilonN23dGex <= 0
-12 <= epsilonN23dGMPex <= 0
-12 <= ethidiumex <= 0
-12 <= ETNex <= 0
-12 <= F1Pex <= 0
-12 <= f5dUex <= 0
-12 <= f5dUMPex <= 0
-12 <= F6Pex <= 0
-12 <= FADex <= 12
-12 <= FADH2ex <= 0
-12 <= FALDex <= 0
-12 <= FAPydAex <= 0
-12 <= FAPydAMPex <= 0
-12 <= FAPydGex <= 0
-12 <= FAPydGMPex <= 0
-12 <= FDPex <= 0
-20 <= FE2ex <= 0
-20 <= FE3ex <= 20
-12 <= FMETex <= 0
-12 <= FMNex <= 0
-12 <= FMNH2ex <= 0
FORex = 0
-12 <= FRUex <= 0
-12 <= FTHF10ex <= 0
-12 <= FTHF5ex <= 0
-12 <= FUC1Pex <= 0
-12 <= G1Pex <= 0
-12 <= G23CMPex <= 0
-12 <= G2Pex <= 0
-12 <= G3MPex <= 0
-12 <= G3Pex <= 0
-12 <= G6Pex <= 0
-12 <= GalCer16ex <= 0
-12 <= GalDAG160ex <= 0
-12 <= GalDAG161ex <= 0
-12 <= GalDAG181ex <= 0
-12 <= GalGalCer16ex <= 0
-12 <= GalGalDAG160ex <= 0
-12 <= GalGalDAG161ex <= 0
-12 <= GalGalDAG181ex <= 0
-12 <= GalGalGalDAG160ex <= 0
-12 <= GalGalGalDAG161ex <= 0
-12 <= GalGalGalDAG181ex <= 0
-12 <= GalGalGlcCer16ex <= 0
-12 <= GalGalGlcDAG160ex <= 0
-12 <= GalGalGlcDAG161ex <= 0
-12 <= GalGalGlcDAG181ex <= 0
-12 <= GalGlcCer16ex <= 0
-12 <= GalGlcDAG160ex <= 0
-12 <= GalGlcDAG161ex <= 0
-12 <= GalGlcDAG181ex <= 0
-12 <= GDPex <= 0
-12 <= GDPTPex <= 0
-12 <= GLex <= 12
-12 <= GL2Pex <= 0
-12 <= GL3Pex <= 0
-12 <= GLCex <= 12
-12 <= GlcCer16ex <= 0
-12 <= GlcDAG160ex <= 0
-12 <= GlcDAG161ex <= 0
-12 <= GlcDAG181ex <= 0
-12 <= GlcGalCer16ex <= 0
-12 <= GlcGalDAG160ex <= 0
-12 <= GlcGalDAG161ex <= 0
-12 <= GlcGalDAG181ex <= 0
-12 <= GlcGalGalDAG160ex <= 0
-12 <= GlcGalGalDAG161ex <= 0
-12 <= GlcGalGalDAG181ex <= 0
-12 <= GlcGalGlcDAG160ex <= 0
-12 <= GlcGalGlcDAG161ex <= 0
-12 <= GlcGalGlcDAG181ex <= 0
-12 <= GlcGlcDAG160ex <= 0
-12 <= GlcGlcDAG161ex <= 0
-12 <= GlcGlcDAG181ex <= 0
-12 <= GlcGlcGlcDAG160ex <= 0
-12 <= GlcGlcGlcDAG161ex <= 0
-12 <= GlcGlcGlcDAG181ex <= 0
-12 <= glchm5dCex <= 0
-12 <= glchm5dCMPex <= 0
-12 <= glchm5dUex <= 0
-12 <= glchm5dUMPex <= 0
-12 <= GLNex <= 0
-12 <= GlnGlnex <= 12
-12 <= GLUex <= 0
-12 <= gludTex <= 0
-12 <= gludTMPex <= 0
-12 <= GluGluex <= 12
-12 <= GLYex <= 0
-12 <= GlyGlyex <= 12
-12 <= Gmex <= 0
-12 <= GmMPex <= 0
-12 <= GMPex <= 0
-12 <= GMP_Lysineex <= 0
-12 <= GMP_Morex <= 0
-12 <= GNex <= 12
-12 <= GPCHex <= 0
-12 <= GPEex <= 0
-12 <= GPGex <= 0
-12 <= GPIex <= 0
-12 <= GPSex <= 0
-12 <= GSNex <= 0
-12 <= GTPex <= 0
-20 <= Hex <= 20
-20 <= H2Oex <= 20
-20 <= H2O2ex <= 20
-12 <= HCYSex <= 0
-12 <= HDCAex <= 12
-12 <= HDCEAex <= 12
-12 <= HDECEG3Pex <= 0
-12 <= HDECG3Pex <= 0
-12 <= HISex <= 0
-12 <= HisHisex <= 12
-12 <= hm3dCex <= 0
-12 <= hm3dCMPex <= 0
-12 <= hm5dCex <= 0
-12 <= hm5dCMPex <= 0
-12 <= hm5dUex <= 0
-12 <= hm5dUMPex <= 0
-12 <= hm5hydantoin_dRiboseex <= 0
-12 <= hm5hydantoin_dRiboseMPex <= 0
-12 <= ho5dCex <= 0
-12 <= ho5dCMPex <= 0
-12 <= ho5dUex <= 0
-12 <= ho5dUMPex <= 0
-12 <= ho5ho6dTex <= 0
-12 <= ho5ho6dTMPex <= 0
-12 <= ho5ho6dUex <= 0
-12 <= ho5ho6dUMPex <= 0
-12 <= ho5hydantoin_dRiboseex <= 0
-12 <= ho5hydantoin_dRiboseMPex <= 0
-12 <= HPHPYRex <= 0
-20 <= hydrideex <= 0
-20 <= hydroxyl_radicalex <= 0
-12 <= HYXNex <= 12
-12 <= IDPex <= 0
-12 <= ILEex <= 0
-12 <= IleIleex <= 12
-12 <= IMPex <= 0
-12 <= INOSex <= 12
-12 <= INSex <= 0
-12 <= ITPex <= 0
-20 <= Kex <= 20
-12 <= k2Cex <= 0
-12 <= k2CMPex <= 0
-12 <= LACex <= 0
-12 <= LEUex <= 0
-12 <= LeuLeuex <= 12
-12 <= LIPOATEex <= 12
-12 <= LIPOYLAMPex <= 0
-12 <= LIPOYLLYSex <= 0
-12 <= lomefloxacinex <= 0
-12 <= LYSex <= 0
-12 <= LysLysex <= 12
-12 <= m1dAex <= 0
-12 <= m1dAMPex <= 0
-12 <= m1dGex <= 0
-12 <= m1dGMPex <= 0
-12 <= m1Gex <= 0
-12 <= m1GMPex <= 0
-12 <= m2Gex <= 0
-12 <= m2GMPex <= 0
-12 <= m3dAex <= 0
-12 <= m3dAMPex <= 0
-12 <= m3dCex <= 0
-12 <= m3dCMPex <= 0
-12 <= m3dGex <= 0
-12 <= m3dGMPex <= 0
-12 <= m3dTex <= 0
-12 <= m3dTMPex <= 0
-12 <= m4dCex <= 0
-12 <= m4dCMPex <= 0
-12 <= m5Cex <= 0
-12 <= m5CMPex <= 0
-12 <= m5dCex <= 12
-12 <= m5dCMPex <= 0
-12 <= m62Aex <= 0
-12 <= m62AMPex <= 0
-12 <= m6Aex <= 0
-12 <= m6AMPex <= 0
-12 <= m6dAex <= 0
-12 <= m6dAMPex <= 0
-12 <= m7dAex <= 0
-12 <= m7dAMPex <= 0
-12 <= m7dGex <= 0
-12 <= m7dGMPex <= 0
-12 <= m7Gex <= 0
-12 <= m7GMPex <= 0
-12 <= m8dGex <= 0
-12 <= m8dGMPex <= 0
-12 <= MALex <= 0
-12 <= MANex <= 0
-12 <= MAN1Pex <= 0
-12 <= MAN6Pex <= 0
-12 <= METex <= 0
-12 <= METHFex <= 0
-12 <= methylRedex <= 0
-12 <= METHYL_OCDCEAex <= 0
-12 <= MetMetex <= 12
-12 <= METSOXRLex <= 0
-12 <= METSOXSLex <= 0
-12 <= METTHFex <= 0
-20 <= MGex <= 20
-20 <= MNex <= 20
-20 <= MOex <= 0
-20 <= MOBDex <= 20
-12 <= MONOBUTYRINex <= 0
-12 <= MORex <= 0
-12 <= MTHFex <= 0
-12 <= n2dAex <= 0
-12 <= n2dAMPex <= 0
-20 <= NAex <= 20
-12 <= NAADex <= 0
-12 <= NACex <= 12
-12 <= NADex <= 12
-12 <= NADHex <= 0
-12 <= NADPex <= 12
-12 <= NADPHex <= 0
-12 <= NAMNex <= 0
-20 <= NH3ex <= 0
-20 <= NH4ex <= 0
-20 <= NIex <= 20
-12 <= nisinex <= 0
-12 <= NMNex <= 0
-12 <= norfloxacinex <= 0
-12 <= novobiocinex <= 0
-12 <= NRex <= 0
-20 <= O2ex <= 20
-12 <= O2edCex <= 0
-12 <= O2edCMPex <= 0
-12 <= O2edTex <= 0
-12 <= O2edTMPex <= 0
-12 <= O2mdCex <= 0
-12 <= O2mdCMPex <= 0
-12 <= O2mdTex <= 0
-12 <= O2mdTMPex <= 0
-12 <= O4edTex <= 0
-12 <= O4edTMPex <= 0
-12 <= O4mdTex <= 0
-12 <= O4mdTMPex <= 0
-12 <= O6edGex <= 0
-12 <= O6edGMPex <= 0
-12 <= O6mdGex <= 0
-12 <= O6mdGMPex <= 0
-12 <= OAex <= 0
-12 <= OCDCAex <= 12
-12 <= OCDCEAex <= 12
-12 <= ODEC11GPex <= 0
-12 <= ODECG3Pex <= 0
-20 <= OHex <= 0
-12 <= ORNex <= 0
-12 <= oxo2dAex <= 0
-12 <= oxo2dAMPex <= 0
-12 <= oxo8dAex <= 0
-12 <= oxo8dAMPex <= 0
-12 <= oxo8dGex <= 0
-12 <= oxo8dGMPex <= 0
-12 <= oxoethyl7dGex <= 0
-12 <= oxoethyl7dGMPex <= 0
-12 <= PA120ex <= 0
-12 <= PA140ex <= 0
-12 <= PA141ex <= 0
-12 <= PA160ex <= 0
-12 <= PA161ex <= 0
-12 <= PA180ex <= 0
-12 <= PA181ex <= 0
-12 <= PAPex <= 0
-12 <= pCHOLex <= 0
-12 <= PDLA5Pex <= 0
-12 <= PEPex <= 0
-12 <= PG160ex <= 0
-12 <= PG161ex <= 0
-12 <= PG181ex <= 0
-12 <= PGP160ex <= 0
-12 <= PGP161ex <= 0
-12 <= PGP181ex <= 0
-12 <= PHEex <= 0
-12 <= PhePheex <= 12
-12 <= PHPYRex <= 0
-20 <= PIex <= 20
-12 <= PL5Pex <= 12
-12 <= POE_SORBITANex <= 0
-12 <= PPAAex <= 0
-12 <= PPGPPex <= 0
-12 <= PPHGex <= 0
-20 <= PPIex <= 0
-12 <= PPIXex <= 0
-12 <= PROex <= 0
-12 <= proflavineex <= 0
-12 <= ProProex <= 12
-12 <= PRPPex <= 0
-12 <= pSERex <= 0
-12 <= PSIURIex <= 0
-12 <= PSIURIMPex <= 0
-12 <= pTHRex <= 0
-12 <= PTRCex <= 12
-12 <= pTYRex <= 0
-12 <= putdTex <= 0
-12 <= putdTMPex <= 0
-12 <= PYDXex <= 12
-12 <= PYRex <= 0
-12 <= R1Pex <= 0
-12 <= R5Pex <= 0
-12 <= RBFLVRDex <= 0
-12 <= rhodamine6Gex <= 0
-12 <= RIBex <= 0
-12 <= RIBFLVex <= 12
-12 <= RL5Pex <= 0
-12 <= s2Uex <= 0
-12 <= s2UMPex <= 0
-12 <= s4Uex <= 0
-12 <= s4UMPex <= 0
-12 <= S7Pex <= 0
-12 <= SALAex <= 0
-12 <= SBTex <= 0
-12 <= SBT6Pex <= 0
-12 <= SERex <= 0
-12 <= SerSerex <= 12
-12 <= SNGLYPex <= 0
-20 <= SO2ex <= 0
-20 <= SO4ex <= 20
-12 <= SPMDex <= 12
-12 <= STEROLex <= 12
-12 <= T3P1ex <= 0
-12 <= T3P2ex <= 0
-12 <= TAG16Pex <= 0
-12 <= TAG6Pex <= 0
-12 <= tBHPex <= 0
-12 <= tBOHex <= 0
-12 <= TDECEG3Pex <= 0
-12 <= TDECG3Pex <= 0
-12 <= TDPex <= 0
-12 <= tetraphenylphosphoniumex <= 0
0 <= THFex <= 12
-12 <= THMex <= 12
-12 <= THMPex <= 0
-12 <= THRex <= 0
-12 <= ThrThrex <= 12
-12 <= THYex <= 12
-12 <= TPPex <= 0
-12 <= TRIBUTYRINex <= 0
-12 <= TRILAURINex <= 0
-12 <= TRIMYRISTINex <= 0
-12 <= TRIOLEINex <= 0
-12 <= TRIPALMITINex <= 0
-12 <= TRISTEARINex <= 0
-12 <= TRI_HDCEA_INex <= 0
-12 <= TRI_TTDCEA_INex <= 0
-12 <= TRPex <= 0
-12 <= TrpTrpex <= 12
-12 <= TTDCAex <= 12
-12 <= TTDCEAex <= 12
-12 <= TTPex <= 0
-12 <= TWEEN80ex <= 12
-12 <= TYRex <= 0
-12 <= TyrTyrex <= 12
-12 <= U23CMPex <= 0
-12 <= U3MPex <= 0
-12 <= UDPex <= 0
-12 <= UDPGex <= 0
-12 <= UDPGALex <= 0
-12 <= UDPGALFex <= 0
-12 <= Umex <= 0
-12 <= UmMPex <= 0
-12 <= UMPex <= 0
-12 <= URAex <= 12
-12 <= URIex <= 0
-12 <= UTPex <= 12
-12 <= VALex <= 0
-12 <= ValValex <= 12
-12 <= vinyl_chlorideex <= 0
-12 <= X5Pex <= 0
-12 <= XANex <= 12
-12 <= XTSNex <= 0
-20 <= ZNex <= 20
MG_124_MONOMERex = 0
MG_287_MONOMERex = 0
MG_124_MONOMER_oxex = 0
MG_229_231_TETRAMERex = 0
MG_229_231_TETRAMER_oxex = 0
MG_287_MONOMER_ACPex = 0
MG_287_MONOMER_ddcaACPex = 0
MG_287_MONOMER_hdeACPex = 0
MG_287_MONOMER_myrsACPex = 0
MG_287_MONOMER_ocdcaACPex = 0
MG_287_MONOMER_octeACPex = 0
MG_287_MONOMER_palmACPex = 0
MG_287_MONOMER_tdeACPex = 0
MG_427_DIMERex = 0
MG_427_DIMER_oxex = 0
MG_454_DIMERex = 0
MG_454_DIMER_oxex = 0
ADix = 0
-5.495024e+00 <= ADPix <= 0
AHCYSix = 0
-2.197701e+00 <= AMPix <= 0
0 <= ATPix <= 10000000000
cmnm5s2UMPix = 0
-2.197701e+00 <= CMPix <= 0
0 <= CTPix <= 10000000000
DR5Pix = 0
FMETix = 0
FORix = 0
0 <= FTHF10ix <= 10000000000
-5.495024e+00 <= GDPix <= 0
-4.626739e-03 <= GmMPix <= 0
-2.197701e+00 <= GMPix <= 0
GNix = 0
0 <= GTPix <= 10000000000
-1.978085e+01 <= Hix <= 0
0 <= H2Oix <= 10000000000
-3.084493e-03 <= k2CMPix <= 0
LIPOYLLYSix = 0
-2.313370e-02 <= m1GMPix <= 0
m2GMPix = 0
m62AMPix = 0
m7GMPix = 0
0 <= METix <= 10000000000
NH3ix = 0
NMNix = 0
-1.099005e+01 <= PIix <= 0
-8.790804e+00 <= PPIix <= 0
-1.465134e-01 <= pSERix <= 0
PSIURIMPix = 0
-2.174567e-01 <= pTHRix <= 0
-2.313370e-02 <= pTYRix <= 0
s4UMPix = 0
SNGLYPix = 0
-4.389233e+00 <= THFix <= 0
THYix = 0
-4.626739e-03 <= UmMPix <= 0
-2.197701e+00 <= UMPix <= 0
URAix = 0
0 <= UTPix <= 10000000000
0 <= Growth <= 10000000000
-10000000000 <= Bmssex <= 0
End
