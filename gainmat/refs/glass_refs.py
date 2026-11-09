"""Host-glass properties and Judd-Ofelt (JO) inputs for rare-earth-doped FIBER glasses.

Built Nov 2026 under a 35-minute web-research budget.  The rules match
host_refs.py and phonon_refs.py:

* `verified=True` means the number was seen on the page given in `src` during
  the build.  Many of those pages are SECONDARY sources (a comparison table in a
  later paper).  In that case the primary paper is named as "via ..." or "cites
  ...", and it was NOT opened.
* `verified=False` means the value was recalled or taken over from an existing
  unverified file (ln_refs.py GAP_LAW), or it is a family proxy.  Do not
  calibrate on these values without checking them.
* `None` means nothing usable was found.  No value is invented.

Caveats (read before use)
-------------------------
1. JO parameters spread widely: +-30-50% between labs for the same nominal
   glass, and more with composition (Omega2 especially).  Treat each set as one
   draw from a broad distribution.  Where several sets were found, OMEGA gives
   their median and lists every set in `sets`.
2. Many entries are "family proxies".  Examples: "silicate glass" stands in for
   fiber silica, and sol-gel silica stands in for MCVD silica.  These are marked
   match="proxy".  Bulk-glass values are not fiber-core values.  Al/P/Ge
   codoping in silica fibers changes Omega_2 and the Er/Yb lifetimes noticeably.
3. LOMS.cz (arXiv:2507.01218) could not be downloaded.  loms.cz and
   api.figshare.com were refused by the egress proxy.  The GitHub repo
   robinkrystufek/loms-jo has only UI code (js/redb.js links to
   loms.cz/download/... files, which were unreachable).  Only the per-(ion,host)
   medians printed in the paper's Table 3 were used.  The arXiv HTML layout is
   garbled, so the host assignment of those medians is "best reading".
4. Energy-gap-law parameters come in two forms.  B*exp(-alpha*dE) is Layne
   1977.  W0*exp[-alpha*(dE-2*hw)] is the van Dijk form used in
   arXiv:1910.10609.  Where only the W0 form was seen, `nr_B_s1` is converted
   with B = W0*exp(2*alpha*hw) and labelled "derived".
5. Blocked or unanswered pages (not used): Wiley/Hindawi Zhu & Peyghambarian
   (read via scispace mirror instead), loms.cz, figshare, PMC pages (captcha),
   rp-photonics fluoride/tellurite/chalcogenide pages, en.wikipedia ZBLAN,
   pubs.aip.org (Gomes 2010 Dy:ZBLAN), inis.iaea.org.

Key sources (opened)
--------------------
ZP10   X. Zhu, N. Peyghambarian, "High-power ZBLAN glass fiber lasers: review and
       prospect", Adv. OptoElectron. 2010, 501956, Table 1 + text.
       https://scispace.com/pdf/high-power-zblan-glass-fiber-lasers-review-and-prospect-1l2wpx5qtu.pdf
HUG    M. Hughes, PhD thesis (Southampton 2007), Table 1.1, as arXiv:1412.3063.
       https://arxiv.org/pdf/1412.3063
DY20   S.D. Jackson / Dy mid-IR laser review, arXiv:2003.06030 (Tables 2, 5 + text).
       https://arxiv.org/pdf/2003.06030
KNS19  laser-cooling of silica, arXiv:1910.10609, Supplement Table S1.
NRL10  Nanoscale Res. Lett. (Er in SBA-15), Table 2 (Er Omega in hosts, refs 29-32).
       https://nanoscalereslett.springeropen.com/articles/10.1007/s11671-010-9732-9/tables/2
"""
from __future__ import annotations

from statistics import median


def _v(value, src, verified, **kw):
    d = dict(value=value, src=src, verified=verified)
    d.update(kw)
    return d


_ZP10 = ("X. Zhu, N. Peyghambarian, Adv. OptoElectron. 2010, 501956, Table 1 "
         "(https://scispace.com/pdf/high-power-zblan-glass-fiber-lasers-review-and-"
         "prospect-1l2wpx5qtu.pdf)")
_HUG = ("M. Hughes thesis, Table 1.1 (arXiv:1412.3063, https://arxiv.org/pdf/1412.3063); "
        "GLS/As2S3 values cite its refs [34,43]")
_DY20 = "Dy mid-IR laser review, arXiv:2003.06030 (https://arxiv.org/pdf/2003.06030)"
_KNS19 = ("arXiv:1910.10609 'Laser cooling of silica glass', Supplement Table S1: "
          "W_mp = W0 exp[-alpha (dE - 2 Ep)] (https://arxiv.org/pdf/1910.10609)")
_LAYNE_UNV = ("Layne, Lowdermilk, Weber, PRB 16, 10 (1977), as recalled in "
              "gainmat/refs/ln_refs.py GAP_LAW (NOT re-read; the OSTI abstract "
              "https://www.osti.gov/biblio/7100555 gives no numbers)")

# ======================================================================
# 1. GLASSES
# ======================================================================
GLASSES = {
    "silica": dict(
        composition=_v("100 SiO2", "definition", True),
        n_D=_v(1.458, "Refractive index @0.589 um: " + _ZP10, True),
        density_g_cm3=_v(2.20, _ZP10, True),
        phonon_max_cm1=_v(1100, _ZP10 + "; also arXiv:1910.10609 Table S1; "
                          "Walsh & Barnes, Appl. Phys. B 78, 325 (2004) abstract", True,
                          note="HUG lists 1150"),
        nr_B_s1=_v(2.4e12, "derived from " + _KNS19 + ": W0=7.8e7 s^-1, alpha=4.7e-3 cm, "
                   "Ep=1100 -> B=W0*exp(2*alpha*Ep)", True,
                   note="Layne 1977 silicate B=1.4e12 (unverified recall) is consistent"),
        nr_alpha_cm=_v(4.7e-3, _KNS19, True),
        uv_edge_nm=_v(160, "transmission range 0.16-4.0 um (1 mm, T>10%): " + _ZP10, True),
        ir_edge_um=_v(4.0, "transmission range 0.16-4.0 um (1 mm, T>10%): " + _ZP10, True,
                      note="practical fiber limit ~2.0-2.2 um (HUG lists window 0.16-2.0)"),
        fiber_loss_dB_per_m_min=_v(2e-4, "0.2 dB/km at 1.5 um (commercial): " + _ZP10, True,
                                   note="HUG: 0.19 dB/km at 1.5 um"),
        re_solubility_note=_v(
            "Low RE solubility in pure SiO2.  RE clusters, so Al2O3 or P2O5 codoping is "
            "standard.  Dantas et al. (Braz. J. Phys. 43, 230 (2013), "
            "https://www.redalyc.org/pdf/464/46427890004.pdf) state that silica is unsuitable "
            "for Nd3+ because of clustering.", "see text", True),
    ),
    "aluminosilicate": dict(
        composition=_v("SiO2 + few mol% Al2O3 (MCVD core); e.g. Tm fibers 0.3-11.7 mol% Al2O3; "
                       "Yb fiber 2.00 wt% Al + 2.52 wt% Yb",
                       "Fibers 6, 66 (2018) Table 1 (https://www2.mdpi.com/2079-6439/6/3/66); "
                       "arXiv:2103.02698 Supp. S2", True),
        n_D=_v(1.46, "silica value + ~1e-3 per mol% Al2O3; no tabulated n_D found", False),
        density_g_cm3=_v(None, "not found", False),
        phonon_max_cm1=_v(1100, "assumed = silica (network still SiO2)", False),
        nr_B_s1=_v(None, "use silica", False),
        nr_alpha_cm=_v(None, "use silica", False),
        uv_edge_nm=_v(None, "not found; ~silica", False),
        ir_edge_um=_v(None, "not found; ~silica", False),
        fiber_loss_dB_per_m_min=_v(
            1.8e-2, "Yb aluminosilicate fiber residual background 18 dB/km "
            "(arXiv:2103.02698, Table I); Yb-silica cooling fibers 5-10 dB/km @1200 nm "
            "(arXiv:1910.10609 Table I)", True),
        re_solubility_note=_v(
            "Al codoping disperses RE ions.  Tm 3F4 lifetime rises from 220 us (pure Si) to "
            "542 us (4.8 mol% Al2O3) and 756 us (11.7 mol%).",
            "Fibers 6, 66 (2018), Table 1 (https://www2.mdpi.com/2079-6439/6/3/66)", True),
    ),
    "germanosilicate": dict(
        composition=_v("SiO2 + GeO2 (few-20 mol%) core", "generic; no specific source opened", False),
        n_D=_v(None, "not found", False),
        density_g_cm3=_v(None, "not found", False),
        phonon_max_cm1=_v(1100, "assumed = silica (Si-O stretch dominates)", False),
        nr_B_s1=_v(None, "use silica", False),
        nr_alpha_cm=_v(None, "use silica", False),
        uv_edge_nm=_v(None, "not found", False),
        ir_edge_um=_v(None, "not found", False),
        fiber_loss_dB_per_m_min=_v(2e-4, "standard telecom (Ge-doped core) ~0.2 dB/km: " + _ZP10,
                                   True, note="ZP10 says 'silica'; telecom cores are Ge-doped"),
        re_solubility_note=_v("Ge alone does not prevent RE clustering; usually used together "
                              "with Al (e.g. Tm-Si-Al-Ge fibers in Fibers 6, 66 (2018) Table 1).",
                              "https://www2.mdpi.com/2079-6439/6/3/66", True),
    ),
    "phosphosilicate": dict(
        composition=_v("SiO2 + P2O5 core (MCVD)", "arXiv:1502.02885 (Yb P2O5-doped fiber)", True),
        n_D=_v(None, "not found", False),
        density_g_cm3=_v(None, "not found", False),
        phonon_max_cm1=_v(1300, "P=O stretch ~1300 cm^-1 locally around RE; recall, not verified",
                          False, note="MDPI Materials 12, 1238 lists phosphate 1200-1345"),
        nr_B_s1=_v(None, "not found", False),
        nr_alpha_cm=_v(None, "not found", False),
        uv_edge_nm=_v(None, "not found", False),
        ir_edge_um=_v(None, "not found", False),
        fiber_loss_dB_per_m_min=_v(None, "not found", False),
        re_solubility_note=_v("P codoping dissolves RE well.  Yb 2F5/2 lifetime is 1.45 ms "
                              "(P-silica) vs 0.83 ms (Al-silica).",
                              "arXiv:1502.02885 (https://arxiv.org/pdf/1502.02885)", True),
    ),
    "phosphate": dict(
        composition=_v("Schott LG-750/LG-770 class metaphosphate laser glass (exact mol% "
                       "proprietary)", "https://www.schott.com/products/laser-glasses-and-"
                       "components/technical-details", True),
        n_D=_v(1.516, "n at 1054 nm (not n_D) for LG-750; LG-770 1.500; APG-1 1.526: Schott "
               "technical details page", True, note="value is n(1054 nm), not n_D"),
        density_g_cm3=_v(2.83, "LG-750 (LG-770 2.59, LG-940 Er/Yb 3.04): Schott technical "
                         "details page", True),
        phonon_max_cm1=_v(1200, "MDPI Materials 12, 1238 (2019) introduction "
                          "(https://www.mdpi.com/1996-1944/12/8/1238/html), uncited there; "
                          "same paper's Sec. 3.1 lists 1345", True),
        nr_B_s1=_v(5.4e12, _LAYNE_UNV, False),
        nr_alpha_cm=_v(4.7e-3, _LAYNE_UNV, False),
        uv_edge_nm=_v(None, "not found", False),
        ir_edge_um=_v(None, "not found", False),
        fiber_loss_dB_per_m_min=_v(None, "not found (phosphate fibers typically ~1e-2-1 dB/m; "
                                   "unverified)", False),
        re_solubility_note=_v("Very high RE solubility (heavily Er/Yb-doped glasses, e.g. "
                              "Schott LG-940/950/960).  No numeric limit found.",
                              "Schott technical details page", True),
    ),
    "germanate": dict(
        composition=_v("GeO2-BaO/CaO-Na2O/Li2O (Walsh, NASA) or heavy-metal Bi2O3-GeO2-"
                       "Na2O-Ga2O3", "NASA NTRS 20080015861 "
                       "(https://ntrs.nasa.gov/api/citations/20080015861/downloads/20080015861.pdf); "
                       "MDPI Materials 12, 1238", True),
        n_D=_v(1.65, "n ~1.65 used for GeO2-BaO-Na2O glass (NASA NTRS 20080015861); "
               "Bi-germanate fiber core 2.19 (MDPI Materials 12, 1238)", True,
               note="the two values are for different compositions"),
        density_g_cm3=_v(None, "not found", False),
        phonon_max_cm1=_v(900, "NASA NTRS 20080015861 ('phonon energies extend to only "
                          "900 cm-1'); Bi-germanate 724 (MDPI Materials 12, 1238)", True),
        nr_B_s1=_v(3.4e10, _LAYNE_UNV, False),
        nr_alpha_cm=_v(4.9e-3, _LAYNE_UNV, False),
        uv_edge_nm=_v(280, "GeO2 glass transparent ~280-5000 nm (AFO Research, citing Margaryan "
                      "& Piliavin 1993; https://aforesearch.com/germanate-glasses.html); "
                      "Bi-germanate 510 nm", True),
        ir_edge_um=_v(5.0, "same AFO page; Bi-germanate 5.5 um (MDPI Materials 12, 1238)", True),
        fiber_loss_dB_per_m_min=_v(None, "not found", False),
        re_solubility_note=_v("Good (several mol% RE).  No numeric source opened.", "none", False),
    ),
    "tellurite": dict(
        composition=_v("TZN ~75TeO2-15ZnO-10Na2O (TZNS0: 75TeO2-15ZnO-9.5Na2O-0.5Er2O3)",
                       "Zheng & Zhou, Optoelectron. Lett. 10, 3 (2014) "
                       "(https://www.researching.cn/ArticlePdf/m00050/2014/10/3/OEL-2014-0209.pdf)",
                       True),
        n_D=_v(2.0, "TeO2-ZnO-BaO: n(632)=1.907-2.088; pure TeO2 2.184: Hrabovsky et al. "
               "arXiv:2308.01186 Table 2", True, note="TZN value not found; TZB range midpoint"),
        density_g_cm3=_v(None, "not found", False),
        phonon_max_cm1=_v(750, "'~750 cm-1' for tellurite (Kumi-Barimah et al., Opt. Mater. X 16, "
                          "100201 (2022), https://eprints.whiterose.ac.uk/192019/); DY20 says "
                          "~650", True),
        nr_B_s1=_v(6.3e10, "Reisfeld & Jorgensen compilation as recalled in ln_refs.py "
                   "(not re-read)", False),
        nr_alpha_cm=_v(4.7e-3, "as nr_B_s1", False),
        uv_edge_nm=_v(320, "Tauc-Lorentz gap 3.73-4.06 eV (305-332 nm) for TZB: arXiv:2308.01186", True),
        ir_edge_um=_v(6.0, "bulk IR transmission to ~6 um: arXiv:2308.01186 (its ref 21)", True),
        fiber_loss_dB_per_m_min=_v(0.26, "TZL (TeO2-ZnO-La2O3) fiber 0.26 dB/m @1.55 um: Lousteau "
                                   "et al. SPIE (https://iris.inrim.it/retrieve/756972c5-e7b2-477c-"
                                   "8f8b-40fe890cf0eb/20_Lousteau_Spie%20Photonics%20Europe.pdf); "
                                   "unpurified TZB ~10 dB/m (arXiv:2308.01186)", True),
        re_solubility_note=_v("High RE solubility (several mol%).  Not quantified on an opened page.",
                              "none", False),
    ),
    "ZBLAN": dict(
        composition=_v("53ZrF4-20BaF2-4LaF3-3AlF3-20NaF", "Le Verre Fluore, "
                       "https://leverrefluore.com/?p=112", True),
        n_D=_v(1.499, _ZP10, True),
        density_g_cm3=_v(4.33, _ZP10, True),
        phonon_max_cm1=_v(600, _ZP10 + "; 580 in arXiv:1910.10609 and Poulain et al. PNCS-ESG "
                          "2018; ~500 in Walsh & Barnes 2004 abstract", True),
        nr_B_s1=_v(1.9e5, "derived from " + _KNS19 + ": W0=1.7e4 s^-1, alpha=2.1e-3, Ep=580", True,
                   note="Disagrees strongly with Shinn et al. PRB 27, 6635 (1983) B=1.88e10, "
                        "alpha=5.77e-3 (recalled in ln_refs.py, unverified).  Prefer Shinn after "
                        "checking."),
        nr_alpha_cm=_v(2.1e-3, _KNS19, True),
        uv_edge_nm=_v(220, "0.22-8.0 um (1 mm, T>10%): " + _ZP10 + "; 0.22-7 um (3 mm): LVF", True),
        ir_edge_um=_v(8.0, _ZP10 + " (1 mm); fiber loss rises beyond ~3.8-4 um "
                      "(arXiv:1805.02857, DY20)", True),
        fiber_loss_dB_per_m_min=_v(1e-2, "commercial background 10-100 dB/km (ZP10 text); "
                                   "<10 dB/km (DY20)", True),
        re_solubility_note=_v("LaF3 site readily substituted.  RE doping up to ~10 mol% in "
                              "fluoride cores (Ho 10 mol% InF3 analog: arXiv:1804.08610).",
                              "arXiv:1804.08610", True),
    ),
    "fluoroindate": dict(
        composition=_v("40InF3-20ZnF2-20SrF2-20BaF2 (IFG)", "https://leverrefluore.com/?p=112", True),
        n_D=_v(1.496, "IZSB glass n=1.4951-1.4976 @532 nm: Manzani et al. arXiv:1210.6313", True,
               note="532 nm, not 589 nm"),
        density_g_cm3=_v(None, "not found", False),
        phonon_max_cm1=_v(509, "arXiv:1805.02857 and arXiv:1804.08610 (cite Almeida JNCS 161, 105 "
                          "(1993); Akella JNCS 213, 1 (1997)); 510 in Poulain et al. PNCS-ESG 2018 "
                          "(https://pncs-esg-2018.sciencesconf.org/215424/215424.pdf)", True),
        nr_B_s1=_v(None, "not found", False),
        nr_alpha_cm=_v(None, "not found", False),
        uv_edge_nm=_v(255, "0.255-8 um (3 mm): LVF https://leverrefluore.com/?p=112", True),
        ir_edge_um=_v(8.0, "LVF (bulk 3 mm); fibers usable to ~5.2 um (Poulain 2018)", True),
        fiber_loss_dB_per_m_min=_v(2.2e-3, "minimum 2.2 dB/km: Poulain, Poulain & Cozic, PNCS-ESG "
                                   "2018; Ho-doped core <0.2 dB/m 3.4-4.0 um (arXiv:1804.08610)", True),
        re_solubility_note=_v("Ho3+ 10 mol% core demonstrated (31InF3-30.5BaF2-19ZnF2-9.5SrF2-10HoF3).",
                              "arXiv:1804.08610", True),
    ),
    "chalcogenide": dict(
        composition=_v("Ga-La-S (GLS): 70Ga2S3-30La2S3 (+~2 wt% La2O3)", _HUG + "; Wikipedia GLS",
                       True, note="JO entries also use Ge-Ga-S / Ge-As-S / Ga-Sb-S; see OMEGA"),
        n_D=_v(2.48, "n @700 nm (not n_D): " + _HUG + "; As2S3 2.56; GLS:O 2.404 (Nistreanu 2010)",
               True),
        density_g_cm3=_v(None, "not found on an opened page (GLS ~4.1 recalled, unverified)", False),
        phonon_max_cm1=_v(425, _HUG + " (As2S3 360); Wikipedia GLS ~450; Ga-Sb-S ~306 "
                          "(Mater. Res. Bull. 70, OSTI 22475889)", True),
        nr_B_s1=_v(1.0e6, "ln_refs.py recall 'Schweizer et al. (order of magnitude)'", False,
                   note="Samson et al. PacRim3 1998 (eprints.soton.ac.uk/76653) gave revised "
                        "parameters; PDF not reachable"),
        nr_alpha_cm=_v(2.9e-3, "as nr_B_s1", False),
        uv_edge_nm=_v(500, "window 0.5-10 um: " + _HUG + "; gap 2.6 eV (475 nm) Wikipedia GLS", True),
        ir_edge_um=_v(10.0, _HUG, True),
        fiber_loss_dB_per_m_min=_v(0.5, "chalcogenide background 0.5-1 dB/m: ZP10 text; ~1 dB/m: DY20; "
                                   "GLS 1500 dB/km @1.5 um (HUG)", True,
                                   note="HUG's 0.3 dB/km @5 um is a predicted minimum, not achieved"),
        re_solubility_note=_v("GLS: 'excellent' RE solubility because La is a network constituent.  "
                              "Ge-As-S / As2S3: poor solubility (clustering).", _HUG, True),
    ),
}

# ======================================================================
# 2. OMEGA  (units 1e-20 cm^2)
#    O2/O4/O6 = the listed set, or the median when several sets are listed.
#    match: "exact" (family and composition fit) | "proxy" (related glass)
# ======================================================================


def _om(sets, src, verified, match, note=None):
    """sets: list of (O2, O4, O6, label).  Returns the OMEGA entry with medians."""
    o2 = [s[0] for s in sets]
    o4 = [s[1] for s in sets]
    o6 = [s[2] for s in sets]
    if not sets:
        o2m = o4m = o6m = None
    else:
        o2m, o4m, o6m = median(o2), median(o4), median(o6)
    d = dict(O2=o2m, O4=o4m, O6=o6m, units="1e-20 cm2", src=src, verified=verified,
             match=match, n_sets=len(sets),
             sets=[dict(O2=a, O4=b, O6=c, label=l) for a, b, c, l in sets])
    if len(sets) > 1:
        d["spread"] = dict(O2=(min(o2), max(o2)), O4=(min(o4), max(o4)), O6=(min(o6), max(o6)))
    if note:
        d["note"] = note
    return d


_NRL = ("Nanoscale Res. Lett. 'Highly efficient NIR PL of Er3+ in SBA-15', Table 2 "
        "(https://nanoscalereslett.springeropen.com/articles/10.1007/s11671-010-9732-9/tables/2); "
        "primary refs [29]-[31] not opened")
_OEL14 = ("Zheng & Zhou, Optoelectron. Lett. 10(3) (2014), Table "
          "(https://www.researching.cn/ArticlePdf/m00050/2014/10/3/OEL-2014-0209.pdf); "
          "comparison rows cite Zheng et al. OEL 9, 0461 (2013)")
_LOMS = ("LOMS.cz median, Hrabovsky, Varak, Krystufek arXiv:2507.01218 Table 3 "
         "(https://arxiv.org/abs/2507.01218); layout garbled, host assignment best reading")

OMEGA = {
    # ---------------------------------------------------------------- Er3+
    "Er:silica": _om([(4.23, 1.04, 0.61, "'silicate glass'")], _NRL + "; same row in " + _OEL14,
                     True, "proxy", note="'silicate glass', not pure fiber silica; matches "
                     "ln_refs.py Er:silica (which was unverified)"),
    "Er:aluminosilicate": _om([(5.60, 1.60, 0.61, "'aluminate glass'")], _NRL + "; same row in "
                              + _OEL14, True, "proxy",
                              note="'aluminate glass', NOT Al-codoped silica; no Al-silica set found"),
    "Er:phosphate": _om([(3.89, 1.01, 0.55, "phosphate glass [30]")], _NRL, True, "exact"),
    "Er:tellurite": _om([(5.39, 1.41, 0.93, "TZNS0 75TeO2-15ZnO-9.5Na2O"),
                         (5.05, 1.45, 1.22, "tellurite glass [31], NRL10"),
                         (7.25, 1.46, 2.25, "LOMS median TeO2-ZnO-BaO")],
                        _OEL14 + "; " + _NRL + "; " + _LOMS, True, "exact"),
    "Er:ZBLAN": _om([(2.37, 1.24, 0.73, "ZBLAN, Bullock et al. 1997"),
                     (2.91, 1.27, 1.11, "'fluoride glass' [29], NRL10"),
                     (3.12, 1.49, 1.15, "ZBLALiP, Mortier et al. 2003")],
                    "arXiv:physics/0611197 (https://arxiv.org/pdf/physics/0611197) citing Bullock "
                    "1997 and Mortier 2003; " + _NRL, True, "exact",
                    note="ln_refs.py Er:ZBLAN 2.54/1.39/0.97 (unverified) lies inside this spread"),
    "Er:germanate": _om([(5.81, 0.85, 0.28, "germanate glass [29]")], _NRL + "; same row in " + _OEL14,
                        True, "exact"),
    "Er:chalcogenide": _om([(1.287, 0.418, 1.019, "(GeS2)80(Ga2S3)20 + 0.6 mol% Er2S3"),
                            (0.560, 0.787, 0.614, "(GeS2)80(Ga2S3)20 + 1.8 mol% Er2S3")],
                           "Ivanova, Jayasimhadri, Kincl, J. Optoelectron. Adv. Mater. 11, 1269 "
                           "(2009) (https://joam.inoe.ro/articles/optical-properties-of-er-doped-"
                           "ges2-ga2s3-glasses/fulltext)", True, "exact",
                           note="Ge-Ga-S, not GLS.  Values look low for a sulfide; treat as uncertain"),
    # ---------------------------------------------------------------- Yb3+ (no JO)
    "Yb:aluminosilicate": dict(
        O2=None, O4=None, O6=None, units="1e-20 cm2",
        tau_2F52_s=_v(0.83e-3, "arXiv:1502.02885", True,
                      alt="765 us (arXiv:2103.02698 Table I); 823-832 us (arXiv:0809.0229)"),
        sigma_abs_peak_cm2=_v(2.7e-20, "976 nm, arXiv:1502.02885 (2.7 pm^2)", True),
        sigma_em_peak_cm2=_v(3.0e-20, "~976 nm, McCumber, arXiv:1502.02885 (3.0 pm^2)", True),
        peak_nm=_v(976, "arXiv:1502.02885", True),
        src="arXiv:1502.02885 (https://arxiv.org/pdf/1502.02885)", verified=True, match="exact"),
    "Yb:phosphosilicate": dict(
        O2=None, O4=None, O6=None, units="1e-20 cm2",
        tau_2F52_s=_v(1.45e-3, "arXiv:1502.02885", True),
        sigma_abs_peak_cm2=_v(1.4e-20, "974.5 nm, arXiv:1502.02885", True),
        sigma_em_peak_cm2=_v(1.5e-20, "arXiv:1502.02885", True),
        peak_nm=_v(974.5, "arXiv:1502.02885", True),
        src="arXiv:1502.02885", verified=True, match="exact"),
    # ---------------------------------------------------------------- Tm3+
    "Tm:silica": _om([], "not found (Walsh & Barnes, Appl. Phys. B 78, 325 (2004) has them; paywalled)",
                     False, "none"),
    "Tm:germanate": _om([(6.14, 1.54, 0.87, "GeO2-BaO/CaO-Na2O/Li2O, Walsh et al."),
                         (2.81, 0.29, 0.86, "GeO2-PbO, Wachtler 1998"),
                         (5.55, 2.03, 1.26, "GeO2-PbO-Nb2O5, Balda 2005"),
                         (3.53, 1.38, 0.71, "GeO2-BaO-K2O, Reisfeld & Eckstein 1975"),
                         (6.7, 1.04, 0.53, "germanate, Eyal 1987")],
                        "NASA NTRS 20080015861 (Walsh et al., Tm germanate) Table 2 "
                        "(https://ntrs.nasa.gov/api/citations/20080015861/downloads/20080015861.pdf); "
                        "LOMS median for the same glass 6.37/1.49/1.22 (" + _LOMS + ")",
                        True, "exact",
                        note="Walsh warns that some older Reisfeld values are tau_lambda, not Omega; "
                             "ion per literature row not stated (Tm context)"),
    "Tm:ZBLAN": _om([], "not found", False, "none"),
    "Tm:tellurite": _om([(4.74, 1.20, 1.31, "70TeO2-10Na2O-10ZnO-10MoO3"),
                         (6.8, 2.0, 2.2, "85TeO2-15WO3")],
                        "Zheng & Man, MEICI 2015 (https://www.atlantis-press.com/article/25198.pdf); "
                        "TeO2-WO3 paper https://core.ac.uk/download/52922062.pdf Table 2",
                        True, "proxy", note="no TZN-specific Tm set found"),
    # ---------------------------------------------------------------- Ho3+
    "Ho:ZBLAN": _om([(0.374, 2.77, 2.14, "ZBLAN fiber 0.5 mol% Ho")],
                    "Chen et al., Opt. Mater. 97, 109351 (2019) (https://core.ac.uk/download/343508826.pdf)",
                    True, "exact", note="exponents garbled in the extracted text (O2 read as "
                    "3.74e-21); check against the PDF"),
    "Ho:silica": _om([(4.480, 2.292, 2.827, "sol-gel silica, Rejikumar 2010"),
                      (2.05, 1.55, 1.68, "sol-gel silica 1000 C, Khan 2018")],
                     "J. Optoelectron. Adv. Mater. 12, 1065 (2010) (https://joam.inoe.ro/articles/"
                     "structural-and-spectroscopic-characterization-of-ho3-in-sol-gel-silica/fulltext); "
                     "Adv. Eng. Res. 178 (Atlantis 2018, https://www.atlantis-press.com/article/55909754.pdf)",
                     True, "proxy", note="sol-gel bulk, not MCVD fiber"),
    "Ho:fluoroindate": _om([], "not found (Gomes et al. Opt. Mater. 60, 618 (2016) not opened)",
                           False, "none"),
    # ---------------------------------------------------------------- Pr3+
    "Pr:ZBLAN": _om([], "not found", False, "none"),
    "Pr:silica": _om([], "not found", False, "none"),
    "Pr:chalcogenide": _om([(7.3, 6.2, 3.9, "GLS, Hewak et al. IEEE PTL 6, 5 (1994)"),
                            (0.614, 2.808, 0.738, "70Ga2S3-30La2S3 + 2.95 wt% O, Nistreanu 2010")],
                           "Nistreanu et al., J. Optoelectron. Adv. Mater. 12, 800 (2010) "
                           "(https://joam.inoe.ro/articles/the-spontaneous-emission-probabilities-"
                           "for-ga-la-spr3-glass/fulltext), which also quotes Hewak 1994",
                           True, "exact", note="Ga-La-S.  Factor ~10 disagreement in O2; Pr JO "
                           "fits are notoriously unstable"),
    "Pr:fluoroindate": _om([(0.7, 5.3, 5.0, "InF3-ZnF2-SrF2-BaF2 (IZSB)")],
                           "Manzani et al. arXiv:1210.6313 Table 3", True, "exact",
                           note="Florez et al. JNCS 213-214, 315 (1997) say Pr:fluoroindate needs "
                                "odd-rank terms for a consistent fit"),
    # ---------------------------------------------------------------- Nd3+
    "Nd:silica": _om([], "not found; nearest SiO2-Na2O-Al2O3-B2O3 1.2/2.3/0.8 (Dantas 2013 Table 1)",
                     False, "none"),
    "Nd:phosphate": _om([(4.30, 5.00, 5.60, "Schott LG-770"),
                         (3.30, 5.10, 5.60, "Kigre Q88"),
                         (3.60, 5.00, 5.50, "Hoya LHG-80"),
                         (3.42, 4.09, 4.35, "KBAP phosphate")],
                        "KNU slides https://indico.knu.ac.kr/event/711/contributions/1231/attachments/"
                        "1025/1695/Eakgapon.pdf (citing A. Jose et al. 2022); Dantas et al. Braz. J. "
                        "Phys. 43, 230 (2013) Table 1 (https://www.redalyc.org/pdf/464/46427890004.pdf)",
                        True, "exact"),
    "Nd:ZBLAN": _om([(2.66, 3.05, 4.08, "ZBLAN [42]")],
                    "Dantas et al. Braz. J. Phys. 43, 230 (2013) Table 1, ref [42] not opened",
                    True, "exact"),
    # ---------------------------------------------------------------- Dy3+
    "Dy:ZBLAN": _om([(1.86, 1.42, 2.73, "ZBLAN [40]"),
                     (3.03, 1.32, 2.06, "ZBLAN [41]"),
                     (3.12, 0.95, 2.3, "ZBLAN, Piramidowicz et al. Opt. Mater. 30, 707 (2008)")],
                    _DY20 + " Table 5", True, "exact"),
    "Dy:chalcogenide": _om([(10.53, 3.17, 1.17, "Ge30As10S60 [47]"),
                            (11.9, 3.58, 2.17, "Ge25Ga5S70 [48]"),
                            (8.51, 2.09, 1.60, "Ga-Sb-S, Mater. Res. Bull. 70")],
                           _DY20 + " Table 5; OSTI 22475889 (https://www.osti.gov/biblio/22475889)",
                           True, "exact", note="Ge-As-S / Ge-Ga-S / Ga-Sb-S; no GLS set found"),
    "Dy:tellurite": _om([(4.28, 1.32, 2.53, "TeO2 [46]")], _DY20 + " Table 5", True, "exact"),
    "Dy:fluoroindate": _om([], "not found", False, "none"),
}


def omega(ion: str, family: str):
    """Return (O2, O4, O6) in 1e-20 cm^2 or None."""
    e = OMEGA.get(f"{ion}:{family}")
    if e is None or e["O2"] is None:
        return None
    return e["O2"], e["O4"], e["O6"]
