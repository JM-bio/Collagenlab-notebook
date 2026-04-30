"""
CollagenLab Notebook v2.0
Detailed buffer calculator + SOP library + Lab notes + Voice input
University of Sydney | Corneal Bioengineering Lab
"""

import streamlit as st
import pandas as pd
import numpy as np
import json, os, datetime
from pathlib import Path
import io

st.set_page_config(
    page_title="CollagenLab Notebook v2",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .main { background-color: #f8fafc; }
    .stButton > button {
        border-radius: 8px; border: 1.5px solid #1a6b8a;
        color: #1a6b8a; background: white; font-weight: 600;
        padding: 0.4rem 1.2rem; transition: all 0.2s;
    }
    .stButton > button:hover { background: #1a6b8a; color: white; }
    .note-card {
        background: white; border-left: 4px solid #1a6b8a;
        border-radius: 8px; padding: 1rem 1.2rem; margin: 0.6rem 0;
        box-shadow: 0 1px 2px rgba(0,0,0,0.05);
        border-top: 0.5px solid #e2e8f0;
        border-right: 0.5px solid #e2e8f0;
        border-bottom: 0.5px solid #e2e8f0;
    }
    .result-box {
        background: #f0fdf4; border: 1px solid #86efac;
        border-radius: 8px; padding: 1.2rem; margin-top: 1rem;
    }
    .warn-box {
        background: #fffbeb; border: 1px solid #fcd34d;
        border-radius: 8px; padding: 1rem; margin-top: 0.8rem;
    }
    .info-box {
        background: #eff6ff; border: 1px solid #bfdbfe;
        border-radius: 8px; padding: 1rem; margin-top: 0.8rem;
    }
    .sop-step {
        background: white; border-radius: 8px; padding: 1rem 1.2rem;
        margin: 0.5rem 0; border: 0.5px solid #e2e8f0;
        border-left: 4px solid #7c3aed;
    }
    .sop-critical {
        background: #fef2f2; border-left: 4px solid #dc2626 !important;
    }
    .section-hdr {
        font-size: 1.1rem; font-weight: 700; color: #1e293b;
        border-bottom: 2px solid #e2e8f0; padding-bottom: 0.3rem;
        margin: 1.2rem 0 0.8rem;
    }
    .tag {
        display: inline-block; background: #e0f2fe; color: #0369a1;
        border-radius: 20px; padding: 2px 10px; font-size: 0.72rem;
        font-weight: 600; margin-right: 4px;
    }
    .metric-row {
        display: flex; gap: 12px; margin: 0.8rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ── Data persistence ───────────────────────────────────────────────────────────
DATA_DIR = Path("lab_data")
DATA_DIR.mkdir(exist_ok=True)
NOTES_FILE = DATA_DIR / "notes.json"

def load_notes():
    if NOTES_FILE.exists():
        with open(NOTES_FILE) as f:
            return json.load(f)
    return []

def save_notes(notes):
    with open(NOTES_FILE, "w") as f:
        json.dump(notes, f, indent=2, ensure_ascii=False)

if "notes" not in st.session_state:
    st.session_state.notes = load_notes()

# ══════════════════════════════════════════════════════════════════════════════
# BUFFER CALCULATION FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def calc_acetic_acid(vol_ml, conc_mM):
    glacial_M = 17.4
    target_M = conc_mM / 1000
    acetic_ul = (target_M / glacial_M) * vol_ml * 1000
    water_ml = vol_ml - acetic_ul / 1000
    return {
        "Glacial acetic acid (17.4 M)": f"{acetic_ul:.2f} µL",
        "Milli-Q water (to volume)":    f"{water_ml:.2f} mL",
        "Final volume":                  f"{vol_ml} mL",
        "Final concentration":           f"{conc_mM} mM  =  {target_M:.4f} M",
        "Expected pH":                   "~2.8–3.2",
        "Storage":                       "4°C, up to 1 month",
    }

def calc_pbs(vol_ml, conc):
    f = 10 if "10" in conc else 1
    return {
        "NaCl":                         f"{8.0*f*vol_ml/1000:.4f} g",
        "KCl":                          f"{0.2*f*vol_ml/1000:.4f} g",
        "Na₂HPO₄ · 7H₂O":             f"{1.44*f*vol_ml/1000:.4f} g",
        "KH₂PO₄":                      f"{0.24*f*vol_ml/1000:.4f} g",
        "Milli-Q water (to volume)":    f"{vol_ml:.0f} mL",
        "pH target":                    "7.4 ± 0.05",
        "Adjust pH with":               "1M HCl or 1M NaOH",
        "Autoclave / filter sterilise": "121°C 15 min  OR  0.22 µm filter",
        "Storage":                      "RT up to 3 months  |  4°C preferred",
    }

def calc_hepes(vol_ml, hepes_mM, nacl_mM, kcl_mM, glucose_mM, pH):
    return {
        "HEPES (MW 238.30 g/mol)":      f"{(hepes_mM/1000)*(vol_ml/1000)*238.30:.4f} g",
        "NaCl (MW 58.44 g/mol)":        f"{(nacl_mM/1000)*(vol_ml/1000)*58.44:.4f} g",
        "KCl (MW 74.55 g/mol)":         f"{(kcl_mM/1000)*(vol_ml/1000)*74.55:.4f} g",
        "D-Glucose (MW 180.16 g/mol)":  f"{(glucose_mM/1000)*(vol_ml/1000)*180.16:.4f} g" if glucose_mM>0 else "0 g (not added)",
        "Milli-Q water (to volume)":    f"{vol_ml:.0f} mL",
        "pH adjustment":                f"Adjust to pH {pH} with 1M NaOH",
        "Filter sterilise":             "0.22 µm filter",
        "Storage":                      "4°C up to 1 month",
    }

def calc_dialysis(vol_ml, acetic_mM, buffer_changes):
    glacial_M = 17.4
    ul = (acetic_mM/1000/glacial_M)*vol_ml*1000
    total_vol = vol_ml * buffer_changes
    total_ul  = ul * buffer_changes
    return {
        "── Per change ──":             "",
        "Glacial acetic acid":          f"{ul:.2f} µL",
        "Milli-Q water (to volume)":    f"{vol_ml - ul/1000:.2f} mL",
        "── Total (all changes) ──":    "",
        f"Total acetic acid ({buffer_changes}×)": f"{total_ul:.2f} µL",
        f"Total water ({buffer_changes}×)":        f"{total_vol - total_ul/1000:.2f} mL",
        "Dialysis membrane MWCO":       "12–14 kDa  (collagen MW ~300 kDa)",
        "Duration per change":          "8–12h at 4°C",
        "Total dialysis time":          f"{buffer_changes * 10}h (~{buffer_changes} days)",
    }

def calc_neutralization(vol_ml, col_conc, ratio_col, ratio_pbs, ratio_naoh):
    total = ratio_col + ratio_pbs + ratio_naoh
    col_vol  = vol_ml * (ratio_col  / total)
    pbs_vol  = vol_ml * (ratio_pbs  / total)
    naoh_vol = vol_ml * (ratio_naoh / total)
    final_conc = col_conc * (ratio_col / total)
    return {
        f"Collagen stock ({col_conc} mg/mL)": f"{col_vol:.3f} mL",
        "10× PBS":                             f"{pbs_vol:.3f} mL",
        "0.1 M NaOH":                          f"{naoh_vol:.3f} mL",
        "Total volume":                        f"{vol_ml:.2f} mL",
        "Final collagen concentration":        f"~{final_conc:.2f} mg/mL",
        "Target pH after mixing":              "7.2–7.4",
        "Gelation":                            "37°C, 30 min (keep on ice while mixing)",
        "Keep on ice":                         "Yes — collagen gels at RT",
    }

def calc_trypsin(vol_ml, conc_pct, stock_pct):
    dilution = conc_pct / stock_pct
    trypsin_ml = dilution * vol_ml
    pbs_ml = vol_ml - trypsin_ml
    return {
        f"Trypsin stock ({stock_pct}%)":  f"{trypsin_ml:.3f} mL",
        "1× PBS":                          f"{pbs_ml:.3f} mL",
        "Final volume":                    f"{vol_ml:.1f} mL",
        "Final concentration":             f"{conc_pct}% trypsin",
        "Incubation (CEC)":               "3–5 min at 37°C",
        "Storage":                         "-20°C (stock)  |  4°C up to 2 weeks (working)",
    }

def calc_sds_page_gel(vol_ml, acrylamide_pct, gel_type):
    """Resolving or stacking gel recipe."""
    if gel_type == "Resolving":
        acrylamide_vol = (acrylamide_pct / 30) * vol_ml * (30/100)
        # Standard resolving gel (30% acrylamide stock)
        acrylamide_ml = (acrylamide_pct / 30) * vol_ml
        buffer_ml     = vol_ml / 4           # 1.5M Tris pH 8.8, 1/4 volume
        sds_ul        = vol_ml * 10           # 10% SDS, 1/100 volume → µL
        water_ml      = vol_ml - acrylamide_ml - buffer_ml - sds_ul/1000 - 0.05 - 0.1
        return {
            "30% Acrylamide/Bis":           f"{acrylamide_ml:.3f} mL",
            "1.5M Tris-HCl pH 8.8":        f"{buffer_ml:.3f} mL",
            "10% SDS":                      f"{sds_ul:.1f} µL",
            "Milli-Q water":                f"{max(water_ml,0):.3f} mL",
            "10% APS":                      "50 µL (add just before pouring)",
            "TEMED":                        "5 µL (add just before pouring)",
            "Expected separation range":    f"~{int(200/acrylamide_pct*10)}–{int(1000/acrylamide_pct*10)} kDa",
            "Note":                         "For collagen α chains (~180 kDa): use 6–8% gel",
        }
    else:  # Stacking
        return {
            "30% Acrylamide/Bis":           f"{(5/30)*vol_ml:.3f} mL  (5% stacking)",
            "0.5M Tris-HCl pH 6.8":        f"{vol_ml/4:.3f} mL",
            "10% SDS":                      f"{vol_ml*10:.1f} µL",
            "Milli-Q water":                f"{vol_ml - (5/30)*vol_ml - vol_ml/4 - vol_ml*10/1000 - 0.04 - 0.008:.3f} mL",
            "10% APS":                      "40 µL",
            "TEMED":                        "8 µL",
        }

def calc_collagenase(vol_ml, conc_units_ml, stock_units_ml):
    enzyme_vol = (conc_units_ml / stock_units_ml) * vol_ml
    buffer_vol = vol_ml - enzyme_vol
    return {
        f"Collagenase stock ({stock_units_ml} U/mL)": f"{enzyme_vol:.3f} mL",
        "HBSS or PBS":                                  f"{buffer_vol:.3f} mL",
        "Final volume":                                 f"{vol_ml:.1f} mL",
        "Final activity":                               f"{conc_units_ml} U/mL",
        "Incubation":                                   "37°C, 30–60 min with gentle agitation",
        "Filter sterilise":                             "0.22 µm after preparation",
        "Storage":                                      "Prepare fresh; freeze aliquots at -80°C",
    }

def calc_bca_assay(num_standards, sample_volume_ul, sample_dilution):
    bsa_concs = [0, 25, 125, 250, 500, 750, 1000, 1500, 2000][:num_standards]
    reagent_a_ml  = (num_standards + sample_volume_ul/25) * 0.2
    reagent_b_ml  = reagent_a_ml / 50
    working_ml    = reagent_a_ml + reagent_b_ml
    return {
        "BSA standard concentrations":   f"{bsa_concs} µg/mL",
        "BCA Reagent A (approx.)":       f"{reagent_a_ml:.2f} mL",
        "BCA Reagent B (approx.)":       f"{reagent_b_ml:.3f} mL",
        "Working reagent (A:B = 50:1)":  f"{working_ml:.2f} mL",
        "Sample dilution factor":         f"{sample_dilution}× (adjust to fall in standard range)",
        "Incubation":                     "37°C, 30 min  (or RT 2h)",
        "Read at":                        "562 nm",
        "Detection range":                "20–2000 µg/mL",
    }

def calc_elisa_dilution(stock_conc, final_conc, final_vol_ml):
    dilution = stock_conc / final_conc
    sample_vol_ul = (final_vol_ml / dilution) * 1000
    diluent_ul    = final_vol_ml * 1000 - sample_vol_ul
    return {
        "Stock concentration":            f"{stock_conc} µg/mL",
        "Target concentration":           f"{final_conc} µg/mL",
        "Dilution factor":                f"{dilution:.1f}×",
        "Sample/antibody volume":         f"{sample_vol_ul:.2f} µL",
        "Diluent volume":                 f"{diluent_ul:.2f} µL",
        "Final volume":                   f"{final_vol_ml*1000:.0f} µL",
    }

# ══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.markdown("### 🔬 CollagenLab v2")
    st.caption("University of Sydney | Corneal Bioengineering")
    st.divider()

    page = st.radio("Navigation", [
        "📔 Lab Book Entry",
        "📝 Lab Notes",
        "🧪 Buffer Calculator",
        "🛠️ Custom Buffer Builder",
        "📋 SOP Library",
        "📊 Data Dashboard",
        "⚙️ Settings",
    ], label_visibility="collapsed")

    st.divider()
    today = datetime.date.today().isoformat()
    today_notes = [n for n in st.session_state.notes if n["timestamp"].startswith(today)]
    c1, c2 = st.columns(2)
    c1.metric("Today", len(today_notes))
    c2.metric("Total", len(st.session_state.notes))

    openai_key = st.text_input("OpenAI API Key (voice)", type="password",
                                placeholder="sk-...", help="For Whisper voice transcription")

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 2 — LAB NOTES
# ══════════════════════════════════════════════════════════════════════════════
if page == "📝 Lab Notes":
    st.title("📝 Lab Notes")
    st.caption(datetime.datetime.now().strftime("%A, %d %B %Y  |  %H:%M"))

    st.markdown('<div class="section-hdr">New Note</div>', unsafe_allow_html=True)
    c1, c2 = st.columns([2,1])
    with c1:
        note_title = st.text_input("Experiment title",
            placeholder="e.g., Type IV Collagen Extraction Run 4 – Placenta P012")
        note_text  = st.text_area("Note content", height=180,
            placeholder="Protocol step, observation, result, issue...\n\ne.g., Dissolved 50mg tissue in 0.1M acetic acid 4°C 24h.\nSDS-PAGE: clear band at ~180kDa. No fibrinogen contamination detected.")
    with c2:
        note_type     = st.selectbox("Type", ["Observation","Protocol","Result",
                                               "Issue/Troubleshoot","Analysis","Meeting","Other"])
        note_tags_raw = st.text_input("Tags", placeholder="collagen, type-iv, SDS-PAGE")
        note_priority = st.select_slider("Priority",
                                          ["Low","Medium","High","Critical"], value="Medium")
        uploaded_file = st.file_uploader("Attach file",
                                          type=["png","jpg","jpeg","csv","pdf","xlsx"])

    # Voice input
    with st.expander("🎙️ Voice Input (Whisper AI — Korean + English)"):
        st.info("Whisper recognises mixed Korean/English and scientific terms: collagen, PBS, HEPES, SDS-PAGE, fibrinogen, iPSC, bioink, chromatography...")
        try:
            from audio_recorder_streamlit import audio_recorder
            audio_bytes = audio_recorder(
                text="Click to record",
                recording_color="#1a6b8a", neutral_color="#64748b",
                icon_name="microphone", icon_size="2x",
            )
            if audio_bytes and openai_key:
                import openai
                client = openai.OpenAI(api_key=openai_key)
                buf = io.BytesIO(audio_bytes); buf.name = "rec.wav"
                with st.spinner("Transcribing..."):
                    result = client.audio.transcriptions.create(
                        model="whisper-1", file=buf, language="ko",
                        prompt="collagen, Type IV collagen, acetic acid, PBS, HEPES, placenta, corneal endothelial cells, fibrinogen, SDS-PAGE"
                    )
                st.text_area("Transcript (edit before saving)", result.text, height=100, key="tx")
            elif audio_bytes:
                st.warning("Enter OpenAI API key in sidebar to enable transcription.")
        except ImportError:
            st.warning("Install: `pip install audio-recorder-streamlit`")
            st.markdown("**Alternative**: Use device voice keyboard (iOS/Android mic button, or Mac Fn+Fn)")

    sc1, sc2, _ = st.columns([1,1,4])
    with sc1:
        if st.button("💾 Save Note", type="primary"):
            if note_text.strip() or note_title.strip():
                tags = [t.strip() for t in note_tags_raw.split(",") if t.strip()]
                note = {
                    "id": len(st.session_state.notes)+1,
                    "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
                    "title": note_title or "(Untitled)",
                    "text": note_text, "type": note_type,
                    "tags": tags, "priority": note_priority,
                    "has_file": uploaded_file is not None,
                    "filename": uploaded_file.name if uploaded_file else None,
                }
                st.session_state.notes.insert(0, note)
                save_notes(st.session_state.notes)
                if uploaded_file:
                    fd = DATA_DIR/"uploads"; fd.mkdir(exist_ok=True)
                    with open(fd/uploaded_file.name,"wb") as f: f.write(uploaded_file.getbuffer())
                st.success(f"Saved! {note['timestamp']}"); st.rerun()
            else:
                st.warning("Enter title or content.")
    with sc2:
        if st.button("🗑️ Clear"): st.rerun()

    st.markdown('<div class="section-hdr">Previous Notes</div>', unsafe_allow_html=True)
    fc1, fc2, fc3 = st.columns(3)
    with fc1: ft = st.selectbox("Type", ["All","Observation","Protocol","Result","Issue/Troubleshoot","Analysis","Meeting","Other"])
    with fc2: fp = st.selectbox("Priority", ["All","Critical","High","Medium","Low"])
    with fc3: fs = st.text_input("Search", placeholder="keyword...")

    shown = st.session_state.notes
    if ft != "All": shown = [n for n in shown if n["type"]==ft]
    if fp != "All": shown = [n for n in shown if n["priority"]==fp]
    if fs: shown = [n for n in shown if fs.lower() in n["title"].lower() or fs.lower() in n["text"].lower()]

    pc = {"Critical":"#dc2626","High":"#f59e0b","Medium":"#3b82f6","Low":"#6b7280"}
    if not shown:
        st.info("No notes found.")
    for n in shown[:30]:
        col = pc.get(n.get("priority","Medium"),"#64748b")
        tags_html = "".join(f'<span class="tag">{t}</span>' for t in n.get("tags",[]))
        st.markdown(f"""
        <div class="note-card">
          <div style="display:flex;justify-content:space-between;">
            <strong>{n['title']}</strong>
            <span style="font-size:0.72rem;font-weight:700;color:{col};background:{col}18;padding:2px 8px;border-radius:12px;">{n.get('priority','')}</span>
          </div>
          <div style="font-size:0.75rem;color:#94a3b8;">🕐 {n['timestamp'].replace('T',' ')} · {n['type']}</div>
          <div style="font-size:0.9rem;color:#374151;margin:0.4rem 0;white-space:pre-wrap;">{n['text'][:400]}{'...' if len(n['text'])>400 else ''}</div>
          <div>{tags_html}</div>
        </div>""", unsafe_allow_html=True)

    if st.session_state.notes:
        df_exp = pd.DataFrame(st.session_state.notes)
        st.download_button("⬇️ Export all notes (CSV)",
                            df_exp.to_csv(index=False).encode(),
                            file_name=f"lab_notes_{today}.csv", mime="text/csv")

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 2 — BUFFER CALCULATOR (DETAILED)
# ══════════════════════════════════════════════════════════════════════════════
elif page == "🧪 Buffer Calculator":
    st.title("🧪 Buffer Calculator")
    st.caption("Detailed recipes for collagen extraction, purification & corneal tissue engineering")

    buf_category = st.selectbox("Category", [
        "── Collagen Extraction ──",
        "Acetic Acid (collagen solubilization)",
        "Dialysis Buffer (purification)",
        "Collagen Neutralization Mix (gelation / bioink)",
        "Collagenase Solution (tissue digestion)",
        "── Cell Culture ──",
        "PBS — Phosphate-Buffered Saline",
        "HEPES Buffer (CEC / iPSC culture)",
        "Trypsin Working Solution",
        "── Protein Analysis ──",
        "SDS-PAGE Gel Recipe",
        "BCA Protein Assay",
        "ELISA Antibody Dilution",
    ])

    st.divider()
    result = None
    tip = ""
    warning = ""

    # ── Acetic Acid ─────────────────────────────────────────────────────────
    if buf_category == "Acetic Acid (collagen solubilization)":
        st.markdown("**Acetic acid** — primary solvent for collagen extraction from placenta / tissue")
        c1, c2 = st.columns(2)
        with c1: vol  = st.number_input("Target volume (mL)", 10, 5000, 500, 50)
        with c2: conc = st.number_input("Concentration (mM)", 1, 1000, 500,
                                         help="Type I collagen: 500mM (0.5M) | Type IV: 50–100mM")
        result = calc_acetic_acid(vol, conc)
        tip = ("💡 **Type IV collagen** (basement membrane): use **50–100 mM** acetic acid to "
               "preserve the 7S domain and NC1 domain quaternary structure. "
               "Higher concentrations may denature the triple-helix.\n\n"
               "💡 **Type I collagen** (placenta fibrillar): use **0.5 M** (500 mM) for efficient extraction.")
        warning = "⚠️ Use glacial acetic acid in fume hood. Corrosive — wear gloves and eye protection."

    # ── Dialysis ────────────────────────────────────────────────────────────
    elif buf_category == "Dialysis Buffer (purification)":
        st.markdown("**Dialysis buffer** — remove salts and small molecule impurities after collagen extraction")
        c1, c2, c3 = st.columns(3)
        with c1: vol     = st.number_input("Volume per change (mL)", 500, 20000, 5000, 500)
        with c2: conc    = st.number_input("Acetic acid (mM)", 1, 100, 10,
                                            help="Type IV: 10mM | Type I: 50mM")
        with c3: changes = st.number_input("Number of buffer changes", 2, 6, 3)
        result  = calc_dialysis(vol, conc, changes)
        tip = ("💡 **Standard protocol**: 3 changes × 10h = ~30h total at 4°C\n\n"
               "💡 **Membrane MWCO**: Use 12–14 kDa cut-off — collagen MW ~300 kDa (trimer) "
               "so it is fully retained.\n\n"
               "💡 **Final dialysis**: Can dialyse against Milli-Q water (no acetic acid) "
               "for lyophilisation / freeze-drying.")
        warning = "⚠️ Always dialyse at 4°C to prevent collagen degradation and bacterial growth."

    # ── Neutralization ────────────────────────────────────────────────────────
    elif buf_category == "Collagen Neutralization Mix (gelation / bioink)":
        st.markdown("**Neutralization mix** — for collagen hydrogel formation and bioink preparation")
        c1, c2 = st.columns(2)
        with c1:
            vol      = st.number_input("Total final volume (mL)", 0.1, 50.0, 2.0, 0.5)
            col_conc = st.number_input("Collagen stock concentration (mg/mL)", 0.5, 20.0, 3.0, 0.5)
        with c2:
            ratio_col  = st.number_input("Collagen ratio", 1, 10, 8)
            ratio_pbs  = st.number_input("10× PBS ratio", 1, 5, 1)
            ratio_naoh = st.number_input("0.1M NaOH ratio", 1, 5, 1)
        result = calc_neutralization(vol, col_conc, ratio_col, ratio_pbs, ratio_naoh)
        tip = ("💡 **Standard ratio**: 8:1:1 (collagen : 10×PBS : 0.1M NaOH) → final pH ~7.4\n\n"
               "💡 **Bioink**: For 3D printing, target **2–4 mg/mL** for printability. "
               "Too low (<1mg/mL) → no gelation. Too high (>6mg/mL) → nozzle clogging.\n\n"
               "💡 **Corneal construct**: Mix with Type IV collagen (1:1) for basement membrane mimicry.")
        warning = "⚠️ Keep ALL components on ice. Collagen gels at 37°C — work quickly!"

    # ── Collagenase ─────────────────────────────────────────────────────────
    elif buf_category == "Collagenase Solution (tissue digestion)":
        st.markdown("**Collagenase** — enzymatic tissue digestion for cell isolation")
        c1, c2, c3 = st.columns(3)
        with c1: vol        = st.number_input("Working volume (mL)", 1.0, 100.0, 10.0, 1.0)
        with c2: target_u   = st.number_input("Target activity (U/mL)", 10, 2000, 200,
                                               help="Placenta digestion: 200–500 U/mL")
        with c3: stock_u    = st.number_input("Stock activity (U/mL)", 100, 10000, 2000)
        result = calc_collagenase(vol, target_u, stock_u)
        tip = ("💡 **Placenta digestion**: 200–500 U/mL collagenase Type I or IV, 37°C, 60 min\n\n"
               "💡 **Corneal stroma**: 1 mg/mL collagenase A, 37°C, 2–4h with agitation\n\n"
               "💡 **Cell viability**: Monitor by trypan blue every 30 min. "
               "Over-digestion reduces CEC yield.")
        warning = "⚠️ Prepare fresh on day of use. Filter sterilise (0.22 µm) before use on cells."

    # ── PBS ──────────────────────────────────────────────────────────────────
    elif buf_category == "PBS — Phosphate-Buffered Saline":
        st.markdown("**PBS** — wash buffer, cell suspension, dilution")
        c1, c2 = st.columns(2)
        with c1: vol  = st.number_input("Volume (mL)", 100, 20000, 1000, 100)
        with c2: conc = st.selectbox("Concentration", ["1× (working solution)", "10× (stock solution)"])
        result = calc_pbs(vol, conc)
        tip = ("💡 **For cell washing**: Use 1× PBS at 37°C (pre-warm before washing cells)\n\n"
               "💡 **Ca²⁺ / Mg²⁺ free**: This recipe is Ca²⁺/Mg²⁺ free — "
               "suitable for cell detachment. For cell adhesion assays, use PBS with Ca²⁺/Mg²⁺.\n\n"
               "💡 **Verification**: pH 7.4 ± 0.05 at 25°C. Check osmolality: 280–300 mOsm/kg.")

    # ── HEPES ────────────────────────────────────────────────────────────────
    elif buf_category == "HEPES Buffer (CEC / iPSC culture)":
        st.markdown("**HEPES buffer** — for corneal endothelial cell (CEC) and iPSC-CEC culture")
        c1, c2, c3 = st.columns(3)
        with c1:
            vol      = st.number_input("Volume (mL)", 50, 2000, 500, 50)
            hepes_mM = st.number_input("HEPES (mM)", 5, 50, 20,
                                        help="Typical: 10–25 mM for cell culture")
        with c2:
            nacl_mM = st.number_input("NaCl (mM)", 0, 200, 140)
            kcl_mM  = st.number_input("KCl (mM)", 0, 20, 5)
        with c3:
            glucose_mM = st.number_input("D-Glucose (mM)", 0, 50, 10, help="0 = do not add")
            pH         = st.number_input("Target pH", 6.8, 7.8, 7.4, 0.05)
        result = calc_hepes(vol, hepes_mM, nacl_mM, kcl_mM, glucose_mM, pH)
        tip = ("💡 **CEC culture**: Use 20 mM HEPES in OptiMEM-based medium for CO₂-independent buffering\n\n"
               "💡 **iPSC-CEC differentiation**: HEPES-buffered B8 medium with ROCK inhibitor (Y-27632) "
               "at 10 µM for first 24h post-seeding on collagen IV substrate\n\n"
               "💡 **pH drift**: HEPES is stable pH 6.8–8.2; use when CO₂ incubator is unavailable.")
        warning = "⚠️ Do not autoclave HEPES — filter sterilise only (0.22 µm)."

    # ── Trypsin ──────────────────────────────────────────────────────────────
    elif buf_category == "Trypsin Working Solution":
        st.markdown("**Trypsin** — cell detachment for CEC and iPSC-derived cells")
        c1, c2, c3 = st.columns(3)
        with c1: vol       = st.number_input("Working volume (mL)", 1.0, 100.0, 10.0, 1.0)
        with c2: final_pct = st.number_input("Final concentration (%)", 0.01, 0.5, 0.05, 0.01,
                                              help="CEC: 0.02–0.05% | Standard: 0.25%")
        with c3: stock_pct = st.selectbox("Stock concentration (%)", [2.5, 0.5, 0.25])
        result = calc_trypsin(vol, final_pct, stock_pct)
        tip = ("💡 **CEC detachment**: Use **0.02–0.05%** trypsin (3–5 min at 37°C) — "
               "CECs are fragile. Higher concentrations damage tight junctions.\n\n"
               "💡 **Collagen substrate**: Pre-coat with collagen IV (2 µg/cm²) "
               "to improve CEC attachment after re-seeding.\n\n"
               "💡 **ROCK inhibitor**: Add Y-27632 (10 µM) after passaging to improve survival.")
        warning = "⚠️ Neutralise trypsin promptly with complete medium (FBS inhibits trypsin). Do not over-trypsinise CECs."

    # ── SDS-PAGE ─────────────────────────────────────────────────────────────
    elif buf_category == "SDS-PAGE Gel Recipe":
        st.markdown("**SDS-PAGE** — collagen chain analysis (α1, α2 ~130 kDa; pro-α ~180 kDa)")
        c1, c2, c3 = st.columns(3)
        with c1: vol      = st.number_input("Gel volume (mL)", 5, 50, 10, 5,
                                             help="Mini gel: 10mL | Standard: 20mL")
        with c2: gel_type = st.selectbox("Gel type", ["Resolving", "Stacking"])
        with c3: pct      = st.number_input("Acrylamide %", 4, 18, 7,
                                             help="For collagen (130–200 kDa): 6–8%")
        result = calc_sds_page_gel(vol, pct, gel_type)
        tip = ("💡 **Collagen SDS-PAGE**: Use **6–8%** resolving gel for collagen α chains (~130–200 kDa)\n\n"
               "💡 **Sample prep**: Heat at 95°C for 5 min in Laemmli buffer with β-ME "
               "(reduces disulphide bonds in Type IV collagen)\n\n"
               "💡 **Loading**: 5–15 µg collagen per lane. Stain with Coomassie R-250 (1h) "
               "or silver stain for low-abundance samples.\n\n"
               "💡 **Expected bands**: Type I: α1 (~138 kDa) + α2 (~129 kDa) + β dimers (~250 kDa)\n"
               "Type IV: α1/α2 chains (~185 kDa), NC1 domain (~26 kDa)")
        warning = "⚠️ Add APS and TEMED last — gel polymerises immediately."

    # ── BCA ──────────────────────────────────────────────────────────────────
    elif buf_category == "BCA Protein Assay":
        st.markdown("**BCA assay** — quantify total collagen protein concentration")
        c1, c2, c3 = st.columns(3)
        with c1: n_std      = st.number_input("Number of standards", 4, 9, 8)
        with c2: sample_ul  = st.number_input("Sample volume (µL)", 5, 25, 25)
        with c3: dilution   = st.number_input("Sample dilution factor", 1, 100, 10,
                                               help="Dilute to fall within 20–2000 µg/mL range")
        result = calc_bca_assay(n_std, sample_ul, dilution)
        tip = ("💡 **Collagen samples**: Acidic samples (acetic acid) may affect BCA reading. "
               "Dialyse or neutralise before assay.\n\n"
               "💡 **Alternative**: Hydroxyproline assay is collagen-specific "
               "(detects Hyp residue unique to collagen) — more accurate for collagen quantification.\n\n"
               "💡 **Expected yield**: Good collagen extraction: 1–5 mg/g tissue (wet weight)")

    # ── ELISA ────────────────────────────────────────────────────────────────
    elif buf_category == "ELISA Antibody Dilution":
        st.markdown("**ELISA** — Type IV collagen, fibronectin, laminin quantification")
        c1, c2, c3 = st.columns(3)
        with c1: stock_c = st.number_input("Stock concentration (µg/mL)", 0.1, 10000.0, 1000.0)
        with c2: final_c = st.number_input("Target concentration (µg/mL)", 0.001, 100.0, 1.0, 0.1)
        with c3: final_v = st.number_input("Final volume (mL)", 0.1, 50.0, 1.0, 0.1)
        result = calc_elisa_dilution(stock_c, final_c, final_v)
        tip = ("💡 **Type IV collagen ELISA**: Coat plate with 2 µg/mL antibody overnight at 4°C. "
               "Block with 1% BSA/PBS 1h. Sample dilution: 1:10–1:100.\n\n"
               "💡 **Capture antibody**: Anti-Col IV clone CIV22 (Abcam ab6586) "
               "or equivalent monoclonal.\n\n"
               "💡 **Detection**: HRP-conjugated secondary, TMB substrate, read 450nm (ref 570nm).")

    # Display result
    if result and not buf_category.startswith("──"):
        st.markdown('<div class="section-hdr">📋 Recipe</div>', unsafe_allow_html=True)
        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        for k, v in result.items():
            if v == "":
                st.markdown(f"**{k}**")
            else:
                col1, col2 = st.columns([2, 3])
                col1.markdown(f"**{k}**")
                col2.markdown(v)
        st.markdown('</div>', unsafe_allow_html=True)

        if tip:
            st.markdown(f'<div class="info-box">{tip}</div>', unsafe_allow_html=True)
        if warning:
            st.markdown(f'<div class="warn-box">{warning}</div>', unsafe_allow_html=True)

        if st.button("📝 Save this calculation to Lab Notes"):
            note_body = f"Buffer: {buf_category}\nDate: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}\n\nRecipe:\n"
            for k, v in result.items():
                if v: note_body += f"  {k}: {v}\n"
            new_note = {
                "id": len(st.session_state.notes)+1,
                "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
                "title": f"Buffer — {buf_category.split('(')[0].strip()}",
                "text": note_body, "type": "Protocol",
                "tags": ["buffer","calculator"], "priority": "Low",
                "has_file": False, "filename": None,
            }
            st.session_state.notes.insert(0, new_note)
            save_notes(st.session_state.notes)
            st.success("Saved to Lab Notes!")

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 3 — CUSTOM BUFFER BUILDER
# ══════════════════════════════════════════════════════════════════════════════
elif page == "🛠️ Custom Buffer Builder":
    st.title("🛠️ Custom Buffer Builder")
    st.caption("Create, save, and recalculate any buffer recipe — fully editable")

    # ── Persistent custom buffer storage ──────────────────────────────────────
    CUSTOM_BUFFERS_FILE = DATA_DIR / "custom_buffers.json"

    def load_custom_buffers():
        if CUSTOM_BUFFERS_FILE.exists():
            with open(CUSTOM_BUFFERS_FILE) as f:
                return json.load(f)
        # Pre-loaded starter templates (3 originals + 6 requested additions)
        now = datetime.datetime.now().isoformat(timespec="seconds")
        return [
            # ════════════════════════════════════════
            # ORIGINAL 3 TEMPLATES
            # ════════════════════════════════════════

            # ── 1. 0.5M Acetic Acid ──────────────────────────────────────────
            {
                "id": 1,
                "name": "0.5M Acetic Acid — Type IV Collagen Extraction",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 500,
                "notes": (
                    "Acetic acid solubilisation buffer for Type IV collagen extraction from placenta.\n"
                    "• 0.5 M used for initial tissue homogenisation step\n"
                    "• Stir tissue at 4°C for 24–48 h with protease inhibitors (PMSF 1mM, EDTA 5mM)\n"
                    "• For Type IV: some protocols use 0.05–0.1 M to preserve 7S / NC1 domains\n"
                    "• Use in fume hood — glacial acetic acid is corrosive"
                ),
                "components": [
                    {"name": "Glacial acetic acid (17.4 M, MW 60.05 g/mol)", "amount": 14.368, "unit": "mL", "mw": 60.05, "role": "Solvent"},
                    {"name": "Milli-Q water (to final volume)", "amount": 485.632, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "~2.8 (no adjustment needed)", "storage": "4°C, up to 1 month", "created": now,
            },

            # ── 2. 1× PBS pH 7.4 ─────────────────────────────────────────────
            {
                "id": 2,
                "name": "1× PBS pH 7.4",
                "category": "Cell Culture",
                "total_volume_ml": 1000,
                "notes": (
                    "Standard phosphate-buffered saline.\n"
                    "• Adjust pH to 7.4 ± 0.05 with 1M HCl or 1M NaOH\n"
                    "• Autoclave (121°C, 15 min) OR filter sterilise (0.22 µm)\n"
                    "• Ca²⁺/Mg²⁺-free — suitable for cell detachment"
                ),
                "components": [
                    {"name": "NaCl (MW 58.44 g/mol)", "amount": 8.000, "unit": "g", "mw": 58.44, "role": "Salt"},
                    {"name": "KCl (MW 74.55 g/mol)", "amount": 0.200, "unit": "g", "mw": 74.55, "role": "Salt"},
                    {"name": "Na₂HPO₄ (MW 141.96 g/mol)", "amount": 1.440, "unit": "g", "mw": 141.96, "role": "Buffer"},
                    {"name": "KH₂PO₄ (MW 136.09 g/mol)", "amount": 0.240, "unit": "g", "mw": 136.09, "role": "Buffer"},
                    {"name": "Milli-Q water (to final volume)", "amount": 1000.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "7.4", "storage": "RT up to 3 months | 4°C preferred", "created": now,
            },

            # ── 3. HEPES Buffer ───────────────────────────────────────────────
            {
                "id": 3,
                "name": "HEPES Buffer 20mM — Type IV Collagen Extraction",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 500,
                "notes": (
                    "HEPES-based washing/equilibration buffer used during Type IV collagen extraction.\n"
                    "• Used to wash decellularised placenta tissue before acid solubilisation\n"
                    "• Maintains physiological ionic strength without disrupting basement membrane structure\n"
                    "• Adjust pH to 7.4 with 1M NaOH\n"
                    "• Filter sterilise only (0.22 µm) — do NOT autoclave"
                ),
                "components": [
                    {"name": "HEPES (MW 238.30 g/mol)", "amount": 2.383, "unit": "g", "mw": 238.30, "role": "Buffer"},
                    {"name": "NaCl (MW 58.44 g/mol)", "amount": 4.088, "unit": "g", "mw": 58.44, "role": "Salt"},
                    {"name": "KCl (MW 74.55 g/mol)", "amount": 0.186, "unit": "g", "mw": 74.55, "role": "Salt"},
                    {"name": "D-Glucose (MW 180.16 g/mol)", "amount": 0.900, "unit": "g", "mw": 180.16, "role": "Carbon source"},
                    {"name": "Milli-Q water (to final volume)", "amount": 500.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "7.4", "storage": "4°C, up to 1 month", "created": now,
            },

            # ════════════════════════════════════════
            # 6 NEW TEMPLATES (requested additions)
            # ════════════════════════════════════════

            # ── 4. 0.5M Acetic Acid (for dialysis / washing) ─────────────────
            {
                "id": 4,
                "name": "0.5M Acetic Acid — Dialysis Buffer (Type IV Collagen)",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 1000,
                "notes": (
                    "Larger-volume acetic acid for dialysis and washing steps in Type IV collagen purification.\n"
                    "• Use as dialysis buffer: 3× changes of 5L over 48h at 4°C\n"
                    "• Dialysis membrane: 12–14 kDa MWCO\n"
                    "• For final dialysis change: reduce to 0.01–0.05M to lower acidity before lyophilisation\n"
                    "• Scale using the 'Scale to volume' field →"
                ),
                "components": [
                    {"name": "Glacial acetic acid (17.4 M, MW 60.05 g/mol)", "amount": 28.736, "unit": "mL", "mw": 60.05, "role": "Solvent"},
                    {"name": "Milli-Q water (to final volume)", "amount": 971.264, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "~2.8", "storage": "4°C, up to 1 month", "created": now,
            },

            # ── 5. 1.5M Tris-HCl pH 7.6 ─────────────────────────────────────
            {
                "id": 5,
                "name": "1.5M Tris-HCl pH 7.6 — Type IV Collagen Extraction",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 100,
                "notes": (
                    "Tris-HCl stock buffer (pH 7.6) for Type IV collagen extraction workflow.\n"
                    "• Used as base for collagenase digestion buffer (e.g. 50mM Tris + 0.4M NaCl, pH 7.6)\n"
                    "• Also used in ion-exchange chromatography equilibration (e.g. DEAE-cellulose, CM-52)\n"
                    "• Used for pepsin neutralisation step: raise pH to ≥7.5 to inactivate pepsin\n"
                    "⚠️ Dissolve Tris base in ~80 mL water → adjust pH to 7.6 with conc. HCl → bring to 100 mL\n"
                    "⚠️ pH of Tris changes ~0.03 units per °C — calibrate at room temperature"
                ),
                "components": [
                    {"name": "Tris base (MW 121.14 g/mol)", "amount": 18.171, "unit": "g", "mw": 121.14, "role": "Buffer"},
                    {"name": "Milli-Q water (dissolve, ~80 mL first)", "amount": 80.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                    {"name": "Concentrated HCl — titrate to pH 7.6 (approx. 8–9 mL)", "amount": 8.5, "unit": "mL", "mw": 36.46, "role": "pH adjuster"},
                    {"name": "Milli-Q water (to final volume 100 mL)", "amount": 100.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "7.6", "storage": "4°C, up to 3 months", "created": now,
            },

            # ── 6. 1.5M Tris-HCl pH 6.5 ─────────────────────────────────────
            {
                "id": 6,
                "name": "1.5M Tris-HCl pH 6.5 — Type IV Collagen Extraction",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 100,
                "notes": (
                    "Tris-HCl stock buffer (pH 6.5) for SDS-PAGE stacking gel in Type IV collagen purity analysis.\n"
                    "• Used as stacking gel buffer — collagen SDS-PAGE (6–8% resolving gel)\n"
                    "• Expected bands: α1/α2 chains ~185 kDa, NC1 domain ~26 kDa\n"
                    "• Also used in some Type IV collagen chromatography elution buffers\n"
                    "⚠️ Dissolve Tris base in ~80 mL water → adjust pH to 6.5 with conc. HCl → bring to 100 mL\n"
                    "⚠️ More HCl required than pH 7.6 — add dropwise, check pH frequently"
                ),
                "components": [
                    {"name": "Tris base (MW 121.14 g/mol)", "amount": 18.171, "unit": "g", "mw": 121.14, "role": "Buffer"},
                    {"name": "Milli-Q water (dissolve, ~80 mL first)", "amount": 80.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                    {"name": "Concentrated HCl — titrate to pH 6.5 (approx. 10–12 mL)", "amount": 11.0, "unit": "mL", "mw": 36.46, "role": "pH adjuster"},
                    {"name": "Milli-Q water (to final volume 100 mL)", "amount": 100.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "6.5", "storage": "4°C, up to 3 months", "created": now,
            },

            # ── 7. 0.4M Sodium Acetate Buffer ────────────────────────────────
            {
                "id": 7,
                "name": "0.4M Sodium Acetate Buffer pH 4.8 — Type IV Collagen Extraction",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 500,
                "notes": (
                    "Sodium acetate / acetic acid buffer (pH 4.8) used in Type IV collagen extraction workflow.\n"
                    "• Used for hydroxyproline assay (collagen-specific quantification)\n"
                    "• Used as mild acid extraction buffer to preserve NC1 / 7S domains\n"
                    "• pKa of acetic acid = 4.76 → pH 4.8 is within optimal buffering range\n"
                    "• Henderson-Hasselbalch: pH 4.8 ≈ 55% NaOAc + 45% AcOH (molar)\n"
                    "⚠️ Verify pH with calibrated pH meter after mixing"
                ),
                "components": [
                    {"name": "Sodium acetate trihydrate (MW 136.08 g/mol)", "amount": 22.165, "unit": "g", "mw": 136.08, "role": "Buffer"},
                    {"name": "Glacial acetic acid (MW 60.05 g/mol) — titrate to pH 4.8", "amount": 1.15, "unit": "mL", "mw": 60.05, "role": "pH adjuster"},
                    {"name": "Milli-Q water (to final volume)", "amount": 500.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "4.8", "storage": "RT up to 3 months | 4°C preferred", "created": now,
            },

            # ── 8. 5M NaCl Stock ──────────────────────────────────────────────
            {
                "id": 8,
                "name": "5M NaCl Stock Solution — Type IV Collagen Extraction",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 500,
                "notes": (
                    "High-concentration NaCl stock for Type IV collagen salting-out purification.\n"
                    "• Add to crude extract to 0.7 M final → fibrinogen precipitates (DISCARD pellet)\n"
                    "• Increase to 2.5 M final → Type IV collagen precipitates (KEEP pellet)\n"
                    "• Repeat precipitation 2–3× for higher purity\n"
                    "• Stir at RT until fully dissolved — may take 20–30 min\n"
                    "• 5M stock: add volume = (target M × final vol) / 5 to your extract"
                ),
                "components": [
                    {"name": "NaCl (MW 58.44 g/mol)", "amount": 146.100, "unit": "g", "mw": 58.44, "role": "Salt"},
                    {"name": "Milli-Q water (to final volume)", "amount": 500.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "~6.0–7.0 (neutral, no adjustment)", "storage": "RT, up to 6 months", "created": now,
            },

            # ── 9. 6M Urea ────────────────────────────────────────────────────
            {
                "id": 9,
                "name": "6M Urea — Type IV Collagen Extraction / Chromatography",
                "category": "Collagen IV Extraction",
                "total_volume_ml": 500,
                "notes": (
                    "Denaturing urea buffer for Type IV collagen solubilisation and purification chromatography.\n"
                    "• Used in DEAE-cellulose or CM-Sepharose column chromatography running buffer\n"
                    "• Dissolves insoluble Type IV collagen aggregates before column loading\n"
                    "• Combine with 50 mM Tris-HCl pH 7.6 + 0.4 M NaCl for chromatography buffer\n"
                    "⚠️ Warm gently ≤37°C to dissolve — higher temperature increases cyanate formation\n"
                    "⚠️ Urea → cyanate (RT degradation) → carbamylates Lys/N-terminus on collagen\n"
                    "⚠️ ALWAYS prepare fresh on day of use — discard after 24 h\n"
                    "• Optional: deionise with AG 501-X8 mixed-bed resin before use"
                ),
                "components": [
                    {"name": "Urea (MW 60.06 g/mol)", "amount": 180.18, "unit": "g", "mw": 60.06, "role": "Denaturant"},
                    {"name": "Milli-Q water (to final volume)", "amount": 500.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "Neutral — combine with 50 mM Tris-HCl pH 7.6 for pH control",
                "storage": "⚠️ Prepare FRESH — use within 24 h", "created": now,
            },
        ]

    def save_custom_buffers(buffers):
        with open(CUSTOM_BUFFERS_FILE, "w") as f:
            json.dump(buffers, f, indent=2, ensure_ascii=False)

    if "custom_buffers" not in st.session_state:
        st.session_state.custom_buffers = load_custom_buffers()

    UNITS    = ["g", "mg", "µg", "mL", "µL", "mM", "M", "U/mL", "%", "units"]
    ROLES    = ["Buffer", "Salt", "Solvent", "Diluent", "Chelator", "Detergent",
                "Enzyme", "Substrate", "Inhibitor", "Supplement", "Other"]
    CATEGORIES = ["Collagen IV Extraction", "Collagen Extraction", "Purification",
                  "Cell Culture", "Protein Analysis", "Gelation / Bioink", "Custom"]

    # ── Layout: left = buffer list, right = editor ─────────────────────────
    left, right = st.columns([1, 2])

    with left:
        st.markdown('<div class="section-hdr">Saved Buffers</div>', unsafe_allow_html=True)

        # Category filter
        cat_filter = st.selectbox("Filter by category",
                                   ["All"] + CATEGORIES, key="cat_filter")
        shown_bufs = st.session_state.custom_buffers
        if cat_filter != "All":
            shown_bufs = [b for b in shown_bufs if b.get("category") == cat_filter]

        if not shown_bufs:
            st.info("No buffers saved yet.")

        selected_id = st.session_state.get("selected_buffer_id", None)

        for buf in shown_bufs:
            is_sel = buf["id"] == selected_id
            border = "border: 2px solid #1a6b8a;" if is_sel else "border: 0.5px solid #e2e8f0;"
            bg     = "background:#eff6ff;" if is_sel else "background:white;"
            if st.button(f"{'▶ ' if is_sel else ''}{buf['name']}", key=f"sel_{buf['id']}"):
                st.session_state.selected_buffer_id = buf["id"]
                st.rerun()

        st.divider()
        if st.button("➕ New blank buffer", use_container_width=True):
            new_id = max((b["id"] for b in st.session_state.custom_buffers), default=0) + 1
            new_buf = {
                "id": new_id,
                "name": "New Buffer",
                "category": "Custom",
                "total_volume_ml": 100,
                "notes": "",
                "components": [
                    {"name": "Component 1", "amount": 0.0, "unit": "g", "mw": 0.0, "role": "Buffer"},
                    {"name": "Milli-Q water", "amount": 100.0, "unit": "mL", "mw": 18.02, "role": "Diluent"},
                ],
                "pH": "", "storage": "4°C",
                "created": datetime.datetime.now().isoformat(timespec="seconds"),
            }
            st.session_state.custom_buffers.append(new_buf)
            save_custom_buffers(st.session_state.custom_buffers)
            st.session_state.selected_buffer_id = new_id
            st.rerun()

    # ── Right panel: editor ────────────────────────────────────────────────
    with right:
        sel_id  = st.session_state.get("selected_buffer_id")
        buf_idx = next((i for i, b in enumerate(st.session_state.custom_buffers)
                        if b["id"] == sel_id), None)

        if buf_idx is None:
            st.info("👈 Select a buffer from the list, or click **New blank buffer** to create one.")
        else:
            buf = st.session_state.custom_buffers[buf_idx]

            st.markdown('<div class="section-hdr">Edit Buffer</div>', unsafe_allow_html=True)

            # ── Header fields ───────────────────────────────────────────────
            hc1, hc2, hc3 = st.columns([3, 2, 1])
            with hc1:
                buf["name"] = st.text_input("Buffer name", value=buf["name"], key="buf_name")
            with hc2:
                buf["category"] = st.selectbox("Category", CATEGORIES,
                                                index=CATEGORIES.index(buf.get("category","Custom")),
                                                key="buf_cat")
            with hc3:
                buf["total_volume_ml"] = st.number_input("Total vol. (mL)",
                                                           min_value=0.001, value=float(buf["total_volume_ml"]),
                                                           format="%.3f", key="buf_vol")

            mc1, mc2, mc3 = st.columns(3)
            with mc1: buf["pH"]      = st.text_input("pH", value=buf.get("pH",""), placeholder="e.g. 7.4", key="buf_ph")
            with mc2: buf["storage"] = st.text_input("Storage", value=buf.get("storage","4°C"), key="buf_storage")
            with mc3:
                scale_vol = st.number_input("Scale to volume (mL)",
                                             min_value=0.001,
                                             value=float(buf["total_volume_ml"]),
                                             format="%.3f", key="scale_vol",
                                             help="Change this to auto-scale all component amounts")

            scale_factor = scale_vol / buf["total_volume_ml"] if buf["total_volume_ml"] > 0 else 1.0

            buf["notes"] = st.text_area("Notes / protocol", value=buf.get("notes",""),
                                         height=60, key="buf_notes",
                                         placeholder="Add protocol notes, references, tips...")

            # ── Component table ─────────────────────────────────────────────
            st.markdown('<div class="section-hdr">Components</div>', unsafe_allow_html=True)

            # Column headers
            hdr = st.columns([3, 2, 1, 2, 2, 1])
            for h, label in zip(hdr, ["Reagent name", "Amount (base)", "Unit", "MW (g/mol)", "Role", "Del"]):
                h.markdown(f"<small style='color:#64748b;font-weight:600;'>{label}</small>",
                            unsafe_allow_html=True)

            components_to_delete = []
            for ci, comp in enumerate(buf["components"]):
                cols = st.columns([3, 2, 1, 2, 2, 1])
                with cols[0]:
                    comp["name"]   = st.text_input("", value=comp["name"], key=f"cname_{buf_idx}_{ci}",
                                                    label_visibility="collapsed")
                with cols[1]:
                    comp["amount"] = st.number_input("", value=float(comp["amount"]),
                                                      min_value=0.0, format="%.4f",
                                                      key=f"camt_{buf_idx}_{ci}",
                                                      label_visibility="collapsed")
                with cols[2]:
                    comp["unit"]   = st.selectbox("", UNITS,
                                                   index=UNITS.index(comp["unit"]) if comp["unit"] in UNITS else 0,
                                                   key=f"cunit_{buf_idx}_{ci}",
                                                   label_visibility="collapsed")
                with cols[3]:
                    comp["mw"]     = st.number_input("", value=float(comp.get("mw", 0.0)),
                                                      min_value=0.0, format="%.2f",
                                                      key=f"cmw_{buf_idx}_{ci}",
                                                      label_visibility="collapsed")
                with cols[4]:
                    comp["role"]   = st.selectbox("", ROLES,
                                                   index=ROLES.index(comp["role"]) if comp["role"] in ROLES else 0,
                                                   key=f"crole_{buf_idx}_{ci}",
                                                   label_visibility="collapsed")
                with cols[5]:
                    if st.button("🗑", key=f"cdel_{buf_idx}_{ci}"):
                        components_to_delete.append(ci)

            for ci in reversed(components_to_delete):
                buf["components"].pop(ci)
                st.rerun()

            ac1, ac2 = st.columns([1, 5])
            with ac1:
                if st.button("➕ Add row"):
                    buf["components"].append(
                        {"name": "", "amount": 0.0, "unit": "g", "mw": 0.0, "role": "Buffer"})
                    st.rerun()

            # ── Scaled recipe output ─────────────────────────────────────────
            st.markdown('<div class="section-hdr">📋 Scaled Recipe</div>', unsafe_allow_html=True)
            st.markdown(
                f"Scaling **{buf['total_volume_ml']} mL → {scale_vol:.3f} mL** "
                f"(factor: {scale_factor:.4f}×)" if scale_factor != 1.0
                else f"Showing base recipe ({buf['total_volume_ml']} mL)",
                unsafe_allow_html=True
            )

            st.markdown('<div class="result-box">', unsafe_allow_html=True)
            has_molar = False
            for comp in buf["components"]:
                scaled = comp["amount"] * scale_factor
                unit   = comp["unit"]
                name   = comp["name"] or "(unnamed)"
                role   = comp.get("role","")

                # Compute molarity if MW provided
                molar_str = ""
                if comp.get("mw", 0) > 0 and unit == "g" and scale_vol > 0:
                    moles  = scaled / comp["mw"]
                    conc_M = moles / (scale_vol / 1000)
                    molar_str = f" &nbsp;→&nbsp; <span style='color:#0369a1;font-size:0.82rem;'>{conc_M*1000:.2f} mM</span>"
                    has_molar = True
                elif comp.get("mw", 0) > 0 and unit == "mg" and scale_vol > 0:
                    moles  = (scaled / 1000) / comp["mw"]
                    conc_M = moles / (scale_vol / 1000)
                    molar_str = f" &nbsp;→&nbsp; <span style='color:#0369a1;font-size:0.82rem;'>{conc_M*1000:.4f} mM</span>"
                    has_molar = True

                # Format number nicely
                if scaled >= 100:    disp = f"{scaled:.2f}"
                elif scaled >= 1:    disp = f"{scaled:.3f}"
                elif scaled >= 0.01: disp = f"{scaled:.4f}"
                else:                disp = f"{scaled:.6f}"

                role_badge = f"<span style='font-size:0.72rem;background:#f1f5f9;color:#475569;padding:1px 7px;border-radius:10px;margin-left:6px;'>{role}</span>" if role else ""
                st.markdown(
                    f"<div style='display:flex;justify-content:space-between;align-items:center;"
                    f"padding:5px 0;border-bottom:0.5px solid #dcfce7;font-size:0.9rem;'>"
                    f"<span><strong>{name}</strong>{role_badge}</span>"
                    f"<span style='color:#15803d;font-weight:600;'>{disp} {unit}{molar_str}</span>"
                    f"</div>",
                    unsafe_allow_html=True
                )

            # Totals
            st.markdown(
                f"<div style='padding:6px 0;font-size:0.85rem;color:#64748b;margin-top:4px;'>"
                f"Total volume: <strong>{scale_vol:.3f} mL</strong> &nbsp;|&nbsp; "
                f"pH: <strong>{buf['pH'] or 'not set'}</strong> &nbsp;|&nbsp; "
                f"Storage: <strong>{buf['storage']}</strong>"
                f"</div>",
                unsafe_allow_html=True
            )
            st.markdown('</div>', unsafe_allow_html=True)

            if has_molar:
                st.markdown('<div class="info-box" style="font-size:0.82rem;">💡 Blue values show molar concentration calculated from MW and volume. '
                            'Enter MW (g/mol) for each reagent to enable this.</div>',
                            unsafe_allow_html=True)

            # ── Action buttons ───────────────────────────────────────────────
            st.divider()
            bc1, bc2, bc3, bc4 = st.columns(4)

            with bc1:
                if st.button("💾 Save buffer", type="primary"):
                    save_custom_buffers(st.session_state.custom_buffers)
                    st.success("Buffer saved!")

            with bc2:
                if st.button("📝 Save to Lab Notes"):
                    body = f"Buffer: {buf['name']}\nVolume: {scale_vol:.2f} mL\npH: {buf['pH']}\nStorage: {buf['storage']}\n\nRecipe:\n"
                    for comp in buf["components"]:
                        scaled = comp["amount"] * scale_factor
                        body  += f"  {comp['name']}: {scaled:.4f} {comp['unit']}\n"
                    if buf["notes"]:
                        body += f"\nNotes: {buf['notes']}"
                    new_note = {
                        "id": len(st.session_state.notes)+1,
                        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
                        "title": f"Buffer — {buf['name']} ({scale_vol:.1f} mL)",
                        "text": body, "type": "Protocol",
                        "tags": ["buffer", buf.get("category","").lower().replace(" ","-")],
                        "priority": "Low", "has_file": False, "filename": None,
                    }
                    st.session_state.notes.insert(0, new_note)
                    save_notes(st.session_state.notes)
                    st.success("Saved to Lab Notes!")

            with bc3:
                # Export as CSV
                rows = []
                for comp in buf["components"]:
                    scaled = comp["amount"] * scale_factor
                    rows.append({
                        "Reagent": comp["name"],
                        f"Amount ({scale_vol:.1f}mL)": round(scaled, 6),
                        "Unit": comp["unit"],
                        "MW (g/mol)": comp.get("mw",""),
                        "Role": comp.get("role",""),
                    })
                df_exp = pd.DataFrame(rows)
                st.download_button(
                    "⬇️ Export CSV",
                    df_exp.to_csv(index=False).encode(),
                    file_name=f"{buf['name'].replace(' ','_')}_{int(scale_vol)}mL.csv",
                    mime="text/csv",
                )

            with bc4:
                if st.button("🗑️ Delete buffer"):
                    st.session_state.custom_buffers.pop(buf_idx)
                    save_custom_buffers(st.session_state.custom_buffers)
                    st.session_state.selected_buffer_id = None
                    st.rerun()

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 1 — LAB BOOK ENTRY (BIENCO template) — DEFAULT LANDING PAGE
# ══════════════════════════════════════════════════════════════════════════════
elif page == "📔 Lab Book Entry":
    st.markdown("""
    <div style='background:#1a6b8a;color:white;padding:10px 18px;border-radius:8px;
                display:flex;justify-content:space-between;align-items:center;margin-bottom:16px;'>
        <div>
            <div style='font-size:1.1rem;font-weight:700;'>🔬 BIENCO Lab Book</div>
            <div style='font-size:0.78rem;opacity:0.85;'>University of Sydney — Corneal Bioengineering</div>
        </div>
        <div style='font-size:0.75rem;opacity:0.8;font-style:italic;'>Addressing the global challenge of corneal blindness</div>
    </div>
    """, unsafe_allow_html=True)

    # ── Template selector ──────────────────────────────────────────────────────
    template = st.selectbox("📄 Select Lab Book Template", [
        "Day 6–8: Ultracentrifugation & Salt Precipitation",
        "Day 9–10: Pepsin Digestion & Salt Precipitation",
        "Day 13: DEAE Chromatography",
        "Day 15: CM Chromatography",
        "Custom / Blank Template",
    ])

    st.divider()

    # ── Helper to render editable table ───────────────────────────────────────
    def editable_table(key, columns, default_rows=2):
        """Render an editable table via st.data_editor."""
        import pandas as pd
        empty = {c: [""] * default_rows for c in columns}
        df = pd.DataFrame(empty)
        return st.data_editor(df, key=key, num_rows="dynamic",
                              use_container_width=True, hide_index=True)

    # ══════════════════════════════════════════════════════════════════════════
    # COMMON HEADER SECTIONS (all templates)
    # ══════════════════════════════════════════════════════════════════════════

    # ── Document info ──────────────────────────────────────────────────────────
    st.markdown('<div class="section-hdr">📋 Document Information</div>', unsafe_allow_html=True)
    ci1, ci2 = st.columns(2)
    with ci1:
        doc_name  = st.text_input("Document Name", value=template)
        exp_date  = st.date_input("Date of Experiment", value=datetime.date.today())
    with ci2:
        researcher = st.text_input("Researcher(s)", placeholder="e.g., J. Moon")
        start_time = st.text_input("Procedure Start Date and Time",
                                    placeholder="e.g., 30 Apr 2026, 09:00")

    # ── Purpose ───────────────────────────────────────────────────────────────
    st.markdown('<div class="section-hdr">🎯 Purpose</div>', unsafe_allow_html=True)
    aim = st.text_input("Aim", value="Purification of collagen IV from human placenta extract")

    # Scope pre-fills by template
    scope_defaults = {
        "Day 6–8: Ultracentrifugation & Salt Precipitation":
            "Purification and repurification of collagen extracts from Day 6 to Day 8. "
            "This phase utilizes ultracentrifugation and differential salt precipitation to produce intact collagen IV.",
        "Day 9–10: Pepsin Digestion & Salt Precipitation":
            "Extraction of collagen from human placenta tissue from Day 9 to Day 10. "
            "This phase includes enzymatic digestion (pepsin) and salt precipitation.",
        "Day 13: DEAE Chromatography":
            "Purification and repurification of collagen extracts from Day 10 to Day 23. "
            "This phase utilizes Diethylaminoethyl (DEAE) chromatography to produce high-purity collagen IV.",
        "Day 15: CM Chromatography":
            "Purification and repurification of collagen extracts from Day 10 to Day 21. "
            "This phase utilizes ion-exchange chromatography (Carboxymethyl, CM) and lyophilisation to produce high-purity collagen IV.",
        "Custom / Blank Template": "",
    }
    scope = st.text_area("Scope of Experiment",
                          value=scope_defaults.get(template, ""),
                          height=80)

    # ── Sample information ─────────────────────────────────────────────────────
    st.markdown('<div class="section-hdr">🧫 Sample Information</div>', unsafe_allow_html=True)
    sample_defaults = {
        "Day 6–8: Ultracentrifugation & Salt Precipitation":
            "0.03M Tris-HCl (pH 7.6) and 0.02M NaCl dissolved and dialysed collagen extracts",
        "Day 9–10: Pepsin Digestion & Salt Precipitation": "",
        "Day 13: DEAE Chromatography":
            "0.5M Tris-HCl (pH 6.5), 2M Urea dialysed sample",
        "Day 15: CM Chromatography":
            "0.04M Sodium Acetate (pH 4.8), 2M Urea dialysed sample",
        "Custom / Blank Template": "",
    }
    sc1, sc2 = st.columns(2)
    with sc1:
        sample_desc = st.text_input("Sample Description",
                                     value=sample_defaults.get(template, ""))
    with sc2:
        sample_id = st.text_input("Sample ID", placeholder="e.g., PLX-2026-004")

    # ══════════════════════════════════════════════════════════════════════════
    # MATERIALS SECTION
    # ══════════════════════════════════════════════════════════════════════════
    st.markdown('<div class="section-hdr">🧰 Materials Used</div>', unsafe_allow_html=True)

    with st.expander("📦 Consumables", expanded=True):
        consumable_defaults = {
            "Day 6–8: Ultracentrifugation & Salt Precipitation": [
                ["Microcentrifuge Tubes (1.5 mL)", "Eppendorf 000702", "Room 122 C1A", "N/A", "1 bag"],
                ["50 mL Centrifuge Tubes", "Corning CLS430828", "Room 122 C1A / Room 128", "N/A", "2 bags"],
                ["50 mL Tubes", "Corning CLS430790", "Room 122 C1A / Room 128", "N/A", "2 bags"],
                ["Stir Magnetic Bar 5×50 mm", "", "Room 122 C1C", "N/A", "2 ea"],
                ["Plastic Dropper 3 mL", "", "Room 122 C1C / Room 128", "N/A", "10 ea"],
                ["Serological Pipette (5/10/25 mL)", "Corning Costar", "Each bench / Room 128", "N/A", "5 ea"],
            ],
            "Day 9–10: Pepsin Digestion & Salt Precipitation": [
                ["Microcentrifuge Tubes (1.5 mL)", "Eppendorf 000702", "Room 122 C1A", "N/A", "1 bag"],
                ["50 mL Centrifuge Tubes", "Corning CLS430828", "Room 122 C1A / Room 128", "N/A", "2 bags"],
                ["Stir Magnetic Bar 5×50 mm", "", "Room 122 C1C", "N/A", "2 ea"],
                ["Serological Pipette (5/10/25 mL)", "Corning Costar", "Each bench / Room 128", "N/A", "5 ea"],
            ],
            "Day 13: DEAE Chromatography": [
                ["Microcentrifuge Tubes (1.5 mL)", "Eppendorf 000702", "Room 122 C1A", "N/A", "1 bag"],
                ["50 mL Centrifuge Tubes", "Corning CLS430828", "Room 122 C1A / Room 128", "N/A", "2 bags"],
                ["15 mL Tubes", "Corning CLS430790", "Room 122 C1A / Room 128", "N/A", "2 bags"],
                ["Stir Magnetic Bar 5×50 mm", "", "Room 122 C1C", "N/A", "2 ea"],
                ["Plastic Dropper 3 mL", "", "Room 122 C1C / Room 128", "N/A", "10 ea"],
                ["Serological Pipette (5/10/25 mL)", "Corning Costar", "Each bench / Room 128", "N/A", "5 ea"],
            ],
            "Day 15: CM Chromatography": [
                ["Microcentrifuge Tubes (1.5 mL)", "Eppendorf 000702", "Room 122 C1A", "N/A", "1 bag"],
                ["50 mL Centrifuge Tubes", "Corning CLS430828", "Room 122 C1A / Room 128", "N/A", "2 bags"],
                ["15 mL Tubes", "Corning CLS430790", "Room 122 C1A / Room 128", "N/A", "2 bags"],
                ["Stir Magnetic Bar 5×50 mm", "", "Room 122 C1C", "N/A", "2 ea"],
                ["Plastic Dropper 3 mL", "", "Room 122 C1C / Room 128", "N/A", "10 ea"],
                ["Serological Pipette (5/10/25 mL)", "Corning Costar", "Each bench / Room 128", "N/A", "5 ea"],
            ],
        }
        rows = consumable_defaults.get(template, [["", "", "", "", ""]])
        df_cons = pd.DataFrame(rows, columns=["Material Name", "Material Code",
                                               "Location", "Use By Date", "Notes / Amount"])
        consumables_df = st.data_editor(df_cons, key="consumables", num_rows="dynamic",
                                         use_container_width=True, hide_index=True)

    with st.expander("🧪 Reagents and Buffers", expanded=True):
        reagent_defaults = {
            "Day 6–8: Ultracentrifugation & Salt Precipitation": [
                ["NaCl", "AJA465-2.5KG", "", ""],
                ["0.03M Tris-HCl (pH 7.6) 0.2M NaCl", "Working solution", "", "Cold room"],
                ["0.5M Acetic Acid", "Working solution", "1 month (4°C)", "Cold room"],
                ["32% Hydrochloric Acid", "Chemsuppl HA020-500M", "", "Fume hood cabinet"],
            ],
            "Day 9–10: Pepsin Digestion & Salt Precipitation": [
                ["0.5M Acetic Acid pH 2.0", "Ajax Finechem AJA1 2.5L", "1 month (4°C)", "Cold room 4°C (10L)"],
                ["Pepsin", "Merck P7012-1G, −30°C", "Prepare fresh", "Room 122, FZ 122 C2 Door bottom"],
                ["10mM HCl", "Working solution", "RT", "Room 122"],
                ["NaCl", "Ajax Finechem 465-500G", "Manufacturer exp.", "Room 122 Chemical cabinet"],
                ["0.03M Tris (pH 6.5) 2M Urea", "Working solution", "1 month (4°C)", "Cold room"],
            ],
            "Day 13: DEAE Chromatography": [
                ["1M HCl", "Working solution from 32% HCl", "6 months RT", "Room 122 C1C"],
                ["DEAE beads", "Cytiva GE17-0710-01 500mL", "Manufacture date", "Room 122 FR 122 C2 top shelf"],
                ["0.03M Tris-HCl (pH 6.5), 2M Urea", "Working solution", "1 month (4°C)", "Cold room"],
                ["0.04M Sodium Acetate (pH 4.8), 2M Urea", "Working solution", "1 month (4°C)", "Cold room"],
            ],
            "Day 15: CM Chromatography": [
                ["1M HCl", "Working solution from 32% HCl", "", ""],
                ["CM beads (Carboxymethyl)", "Merck", "N/A", ""],
                ["0.04M Sodium Acetate (pH 4.8), 2M Urea", "Working solution", "1 month (4°C)", "Cold room"],
                ["Dialysis flask 10 kDa MWCO 250 mL", "ThermoFisher 87762", "N/A", "Room 122 C1B"],
            ],
        }
        rows_r = reagent_defaults.get(template, [["", "", "", ""]])
        df_reag = pd.DataFrame(rows_r, columns=["Reagent Name", "Material Code",
                                                  "Use By Date", "Notes / Location"])
        reagents_df = st.data_editor(df_reag, key="reagents", num_rows="dynamic",
                                      use_container_width=True, hide_index=True)

    # ── Equipment / Instruments ────────────────────────────────────────────────
    with st.expander("🔧 Equipment & Instruments Used"):
        equip_defaults = {
            "Day 6–8: Ultracentrifugation & Salt Precipitation": [
                ["Ultracentrifuge", "Beckman Coulter Optima XPN-100", "Annual — TSS Managed (K25-MFB Room 124)", "Booking required"],
                ["Multifuge (centrifuge)", "Thermo Scientific 75009915", "Annual — Room 122", ""],
                ["Stir Plate", "IKA C-MAG HS 7", "Condition-based", "Cold room"],
                ["pH Meter", "Mettler-Toledo SevenDirect SD20", "Annual (tech) / Weekly (operator)", "Room 122"],
                ["Scale", "Mettler-Toledo MX205DU", "Autocalibration", "Room 122"],
                ["Auto Dialysis Equipment", "GR-001", "Condition-based", "Cold room"],
            ],
            "Day 9–10: Pepsin Digestion & Salt Precipitation": [
                ["pH Probe", "Mettler-Toledo", "N/A", "Room 122 C1A"],
                ["Stir Plate (digital)", "IKA C-MAG HS 7", "N/A", "Cold room"],
                ["Snakeskin Dialysis Tubing 10 kDa", "ThermoFisher 87762", "N/A", "Room 122 C1B"],
                ["Beakers 500/1000/2000 mL", "Pyrex", "N/A", "Room 122 buffer storage"],
                ["Auto Dialysis Equipment", "GR-001", "N/A", "Cold room"],
            ],
            "Day 13: DEAE Chromatography": [
                ["pH Probe", "Mettler-Toledo SevenDirect SD20", "Daily–weekly (operator)", "Room 122"],
                ["Centrifuge Multifuge X4R", "Thermo Scientific 75009915", "Annual — Room 122", ""],
                ["NanoDrop One (UV-Vis)", "Thermo Scientific", "TSS Managed — K25-MFB G03", "Booking required"],
                ["Chromatography Column 5×30 cm", "Biorad glass", "N/A", "Room 122 D1B"],
                ["Chromatography Column 5×50 cm", "Biorad glass", "N/A", "Room 122 D1B"],
                ["Microsart Mini Vacuum Pump", "Sartorius 16694-2-50-06", "N/A", "Room 122 D1A"],
                ["Microcentrifuge", "Sigma 150021", "N/A", "Room 122 C2A"],
                ["Auto Dialysis System", "GR-001", "N/A", "Cold room"],
            ],
            "Day 15: CM Chromatography": [
                ["pH Probe", "Mettler-Toledo SevenDirect SD20", "Daily–weekly (operator)", "Room 122"],
                ["Centrifuge Multifuge X4R", "Thermo Scientific 75009915", "Annual — Room 122", ""],
                ["NanoDrop One (UV-Vis)", "Thermo Scientific", "TSS Managed — K25-MFB G03", "Booking required"],
                ["Heating Block", "Thermoline TC0401003", "N/A", "Room 122 C2A"],
                ["Microcentrifuge", "Sigma 150021", "N/A", "Room 122 C2A"],
            ],
        }
        rows_eq = equip_defaults.get(template, [["", "", "", ""]])
        df_eq = pd.DataFrame(rows_eq, columns=["Equipment / Instrument", "Model / Serial No.",
                                                 "Next Calibration", "Notes / Location"])
        st.data_editor(df_eq, key="equipment", num_rows="dynamic",
                        use_container_width=True, hide_index=True)

    # ── Test Groups ────────────────────────────────────────────────────────────
    st.markdown('<div class="section-hdr">🔬 Test Groups</div>', unsafe_allow_html=True)
    tg_defaults = {
        "Day 15: CM Chromatography": [
            ["1", "Negative control", "5 µg human fibrinogen", "", ""],
            ["2", "Positive control", "5 µg human collagen IV (Sigma)", "", ""],
            ["3", "Samples", "20 µL of CM chromatography fractions", "", ""],
        ],
        "Day 13: DEAE Chromatography": [
            ["1", "Negative control", "", "", ""],
            ["2", "Positive control", "", "", ""],
        ],
    }
    rows_tg = tg_defaults.get(template, [
        ["1", "Negative control", "", "", ""],
        ["2", "Positive control", "", "", ""],
        ["3", "", "", "", ""],
    ])
    df_tg = pd.DataFrame(rows_tg, columns=["#", "Name", "Description / Justification",
                                             "Colour", "Result"])
    st.data_editor(df_tg, key="test_groups", num_rows="dynamic",
                    use_container_width=True, hide_index=True)

    # ══════════════════════════════════════════════════════════════════════════
    # PRE-CHECK LIST (template-specific)
    # ══════════════════════════════════════════════════════════════════════════
    st.markdown('<div class="section-hdr">✅ Pre-Check List</div>', unsafe_allow_html=True)

    prechecks = {
        "Day 6–8: Ultracentrifugation & Salt Precipitation": [
            "Verify pH meter has been recently calibrated",
            "Confirm ultracentrifuge reservation (Room 124, K25-MFB)",
            "Check pH of each bottle of allocated buffer",
            "Check required number of centrifuge bottles — no visible damage",
            "Pre-label all 50 mL tubes and microcentrifuge tubes",
        ],
        "Day 9–10: Pepsin Digestion & Salt Precipitation": [
            "Check pH of each bottle of allocated buffer",
            "Check required number of centrifuge bottles — no visible damage",
            "Check pH probe accuracy (2-point calibration)",
        ],
        "Day 13: DEAE Chromatography": [
            "Record initial sample pH",
            "Confirm quantity of 0.45 µm filter units available",
            "Verify DEAE bead volume",
            "Confirm Urea buffers are <7 days old and stored at 4°C",
            "Confirm all buffers are at pH 7.0 (0.03M Tris-HCl pH 6.5, 2M Urea)",
            "Confirm NanoDrop, Centrifuge, and Freeze-dryer reservations",
            "Verify pH probe and NanoDrop have been recently calibrated",
            "Pre-label all 50 mL tubes and eppendorf tubes for MS / SDS-PAGE aliquots",
        ],
        "Day 15: CM Chromatography": [
            "Record initial sample pH",
            "Verify CM bead volume",
            "Confirm Urea buffers are <7 days old and stored at 4°C",
            "Confirm all buffers are at pH 4.8–5.2",
            "Confirm NanoDrop, Centrifuge, and Freeze-dryer reservations",
            "Verify pH probe and NanoDrop have been recently calibrated",
            "Pre-label all 15 mL tubes for fractions and eppendorf tubes for MS / SDS-PAGE aliquots",
        ],
        "Custom / Blank Template": [
            "Verify all equipment calibrations",
            "Check buffer pH values",
            "Confirm instrument reservations",
        ],
    }

    checks = prechecks.get(template, [])
    check_states = {}
    for i, step in enumerate(checks):
        col_cb, col_step, col_note = st.columns([0.5, 4, 3])
        with col_cb:
            check_states[i] = st.checkbox("", key=f"chk_{i}")
        with col_step:
            st.markdown(f"<div style='padding:6px 0;font-size:0.9rem;'>{step}</div>",
                         unsafe_allow_html=True)
        with col_note:
            st.text_input("", placeholder="Note / change recorded here",
                           key=f"chknote_{i}", label_visibility="collapsed")

    pc_initials = st.text_input("Completed by (initials and date)",
                                  placeholder="e.g., JM — 30 Apr 2026")

    # ══════════════════════════════════════════════════════════════════════════
    # SOP DATA RECORDING TABLES (template-specific)
    # ══════════════════════════════════════════════════════════════════════════
    st.markdown('<div class="section-hdr">📊 SOP Data Recording</div>', unsafe_allow_html=True)

    sop_notes = st.text_area("Protocol Notes / SOP deviations",
                               height=80,
                               placeholder="Record any deviations from SOP here. "
                               "If following SOP exactly, write 'No deviations'.")

    # ── Day 6–8 tables ────────────────────────────────────────────────────────
    if template == "Day 6–8: Ultracentrifugation & Salt Precipitation":
        st.markdown("#### SOP 8.5.2 — Ultracentrifugation")
        uc1, uc2 = st.columns(2)
        with uc1:
            st.markdown("**Tube weights**")
            df_uc = pd.DataFrame({
                "No. of tubes": [""],
                "Total empty tube weight (g)": [""],
                "Total filled tube weight (g)": [""],
                "Average tube weight (g)": [""],
            })
            st.data_editor(df_uc, key="uc_weights", hide_index=True, use_container_width=True)
        with uc2:
            st.markdown("**Pellet / supernatant**")
            df_pel = pd.DataFrame({
                "Supernatant volume (mL)": [""],
                "Total tube + pellet (g)": [""],
                "Total empty tube (g)": [""],
                "Pellet mass (g)": [""],
            })
            st.data_editor(df_pel, key="uc_pellet", hide_index=True, use_container_width=True)

        st.markdown("#### SOP 8.5.3 — Salt Precipitation")
        df_salt = pd.DataFrame({
            "Step": ["0.7M NaCl (fibrinogen removal)", "2.5M NaCl (collagen precipitation)", "Re-precipitation 1", "Re-precipitation 2"],
            "NaCl added (mL from 5M stock)": ["", "", "", ""],
            "Final NaCl conc. (M)": ["0.7", "2.5", "", ""],
            "Pellet kept? (Y/N)": ["N — discard", "Y — keep", "", ""],
            "Notes": ["", "", "", ""],
        })
        st.data_editor(df_salt, key="salt_precip", hide_index=True, use_container_width=True)

    # ── Day 9–10 tables ───────────────────────────────────────────────────────
    elif template == "Day 9–10: Pepsin Digestion & Salt Precipitation":
        st.markdown("#### SOP 8.2.13.1 — Sample Volume & Tube Weight")
        d9c1, d9c2 = st.columns(2)
        with d9c1:
            st.text_input("Volume of Acetic Acid (mL)", key="d9_aavol")
            st.text_input("Start Time", key="d9_start")
        with d9c2:
            st.text_input("Sample mass / tissue pellet mass (g)", key="d9_mass")

        st.markdown("#### SOP 8.2.14.3 — pH Adjustment to pH 2.0")
        df_ph = pd.DataFrame({
            "Initial pH": [""],
            "Final pH": [""],
            "Volume of HCl used (mL)": [""],
        })
        st.data_editor(df_ph, key="d9_ph", hide_index=True, use_container_width=True)

        st.markdown("#### SOP 8.2.14.5 — Pepsin Digest (1 mg/mL pepsin solution)")
        d9p1, d9p2 = st.columns(2)
        with d9p1:
            st.text_input("Mass of pepsin weighed (mg)", key="d9_pepsin_mass")
            st.text_input("Volume of 10mM HCl to dissolve pepsin (mL)", key="d9_hcl_vol")
        with d9p2:
            st.text_input("Volume of pepsin solution added to sample (mL)", key="d9_pepsin_vol")
            st.text_input("Time pepsin solution added", key="d9_pepsin_time")

        st.markdown("#### SOP 8.2.14.6 — Cold Room Rotator")
        st.text_input("Time placed in cold room", key="d9_coldroom")

        st.markdown("#### Day 10 — SOP 8.4.1.1")
        st.text_input("Time sample collected from cold room", key="d10_collect")

        st.markdown("#### Salt Precipitation Summary")
        df_salt9 = pd.DataFrame({
            "Step": ["1.8M NaCl precipitation", "Reprecipitation"],
            "NaCl added (mL)": ["", ""],
            "Pellet mass (g)": ["", ""],
            "Notes": ["", ""],
        })
        st.data_editor(df_salt9, key="d9_salt", hide_index=True, use_container_width=True)

    # ── Day 13 DEAE tables ────────────────────────────────────────────────────
    elif template == "Day 13: DEAE Chromatography":
        st.markdown("#### SOP 8.9.2.4 — Sample pH Measurement & Adjustment")
        df_deae_ph = pd.DataFrame({
            "": ["pH"],
            "Initial": [""],
            "Final": [""],
            "Volume used (mL)": [""],
        })
        st.data_editor(df_deae_ph, key="deae_ph", hide_index=True, use_container_width=True)

        st.markdown("#### SOP 8.9.2.7–8.9.2.10 — Filtration (0.45 µm)")
        df_filt = pd.DataFrame({
            "Filter unit # (500 mL)": [""],
            "Sample vol. before filtration (mL)": [""],
            "Filtered total volume (mL)": [""],
            "Notes on filtration": [""],
        })
        st.data_editor(df_filt, key="deae_filt", hide_index=True, use_container_width=True)

        st.markdown("#### SOP 8.9.4 — DEAE Bead Mixing")
        df_deae = pd.DataFrame({
            "": ["DEAE-1", "DEAE-2"],
            "DEAE bead volume (mL)": ["", ""],
            "Start time": ["", ""],
        })
        st.data_editor(df_deae, key="deae_beads", hide_index=True, use_container_width=True)

        st.markdown("#### SOP 8.9.4.14 — NanoDrop Absorbance after DEAE Chromatography")
        df_nano = pd.DataFrame({
            "Fraction": [f"F{i}" for i in range(1, 9)],
            "A280 (mg/mL)": [""] * 8,
            "Volume (mL)": [""] * 8,
            "Notes": [""] * 8,
        })
        st.data_editor(df_nano, key="deae_nano", hide_index=True, use_container_width=True)

    # ── Day 15 CM tables ──────────────────────────────────────────────────────
    elif template == "Day 15: CM Chromatography":
        st.markdown("#### SOP 8.10.2.4 — Sample pH")
        st.text_input("Initial sample pH", key="cm_ph")

        st.markdown("#### SOP 8.10.4.2 — Swelled CM Bead Volume")
        st.text_input("Swelled CM bead volume (mL)", key="cm_bead_vol")

        st.markdown("#### SOP 8.10.4.3 — Mixing Sample with CM Beads")
        st.text_input("Start time", key="cm_mix_time")

        st.markdown("#### SOP 8.10.4.4 — Elution Buffers (NaCl gradient)")
        df_elut = pd.DataFrame({
            "Buffer": ["Buffer 1", "Buffer 2", "Buffer 3", "Buffer 4"],
            "Composition": [
                "0.04M Sodium Acetate (pH 4.8), 2M Urea",
                "0.04M Sodium Acetate (pH 4.8), 2M Urea + 0.05M NaCl",
                "0.04M Sodium Acetate (pH 4.8), 2M Urea + 0.10M NaCl",
                "0.04M Sodium Acetate (pH 4.8), 2M Urea + 0.15M NaCl",
            ],
            "Volume prepared (mL)": ["500", "500", "500", "500"],
            "pH verified": ["", "", "", ""],
            "Notes": ["", "", "", ""],
        })
        st.data_editor(df_elut, key="cm_buffers", hide_index=True, use_container_width=True)

        st.markdown("#### Fraction Collection & NanoDrop")
        df_cm_nano = pd.DataFrame({
            "Fraction": [f"F{i}" for i in range(1, 13)],
            "Buffer used": [""] * 12,
            "A280 (mg/mL)": [""] * 12,
            "Volume (mL)": [""] * 12,
            "Notes": [""] * 12,
        })
        st.data_editor(df_cm_nano, key="cm_fractions", hide_index=True,
                        use_container_width=True)

    # ── Custom template ───────────────────────────────────────────────────────
    else:
        st.markdown("#### Protocol Steps / Data Recording")
        df_custom = pd.DataFrame({
            "SOP Step": [""] * 5,
            "Parameter": [""] * 5,
            "Value / Result": [""] * 5,
            "Notes": [""] * 5,
        })
        st.data_editor(df_custom, key="custom_data", num_rows="dynamic",
                        hide_index=True, use_container_width=True)

    # ── Additional observations ────────────────────────────────────────────────
    st.markdown('<div class="section-hdr">📝 Additional Observations / Results</div>',
                unsafe_allow_html=True)
    observations = st.text_area("", height=120,
                                  placeholder="Record additional observations, unexpected results, "
                                  "issues encountered, or next steps...")

    # ══════════════════════════════════════════════════════════════════════════
    # SAVE TO LAB NOTES
    # ══════════════════════════════════════════════════════════════════════════
    st.divider()
    sc1, sc2 = st.columns([1, 4])
    with sc1:
        if st.button("💾 Save Lab Book Entry to Notes", type="primary"):
            # Build structured summary text
            completed_checks = [checks[i] for i, done in check_states.items() if done]
            body = (
                f"BIENCO Lab Book Entry\n"
                f"Template: {template}\n"
                f"Date: {exp_date}\n"
                f"Researcher: {researcher}\n"
                f"Start: {start_time}\n\n"
                f"AIM: {aim}\n"
                f"SCOPE: {scope}\n\n"
                f"SAMPLE: {sample_desc} | ID: {sample_id}\n\n"
                f"PRE-CHECKS COMPLETED ({len(completed_checks)}/{len(checks)}):\n"
                + "\n".join(f"  ✅ {c}" for c in completed_checks)
                + (f"\n\nSOP NOTES:\n{sop_notes}" if sop_notes else "")
                + (f"\n\nOBSERVATIONS:\n{observations}" if observations else "")
            )
            new_note = {
                "id": len(st.session_state.notes) + 1,
                "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
                "title": f"{template} — {exp_date.strftime('%d %b %Y')}",
                "text": body,
                "type": "Protocol",
                "tags": ["lab-book", "BIENCO", "collagen-IV",
                          template.split(":")[0].lower().replace(" ", "-").replace("–", "")],
                "priority": "High",
                "has_file": False,
                "filename": None,
            }
            st.session_state.notes.insert(0, new_note)
            save_notes(st.session_state.notes)
            st.success(f"Lab book entry saved to Lab Notes! ({new_note['timestamp']})")

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 5 — SOP LIBRARY
# ══════════════════════════════════════════════════════════════════════════════
elif page == "📋 SOP Library":
    st.title("📋 SOP Library")
    st.caption("Standard Operating Procedures — Collagen extraction, purification & corneal cell work")

    sop_choice = st.selectbox("Select SOP", [
        "SOP-01: Type IV Collagen Extraction from Human Placenta",
        "SOP-02: Type I Collagen Extraction from Placenta",
        "SOP-03: Collagen Purification by Dialysis",
        "SOP-04: SDS-PAGE for Collagen Purity Analysis",
        "SOP-05: Collagen Quantification (BCA + Hydroxyproline)",
        "SOP-06: Corneal Endothelial Cell (CEC) Culture on Collagen Substrate",
        "SOP-07: Collagen Bioink Preparation for 3D Printing",
        "SOP-08: Fibrinogen Removal from Collagen Preparation",
    ])

    def sop_step(number, title, content, critical=False):
        css_class = "sop-step sop-critical" if critical else "sop-step"
        prefix = "🚨 CRITICAL — " if critical else ""
        st.markdown(f"""
        <div class="{css_class}">
            <strong>Step {number}: {prefix}{title}</strong><br/>
            <span style="font-size:0.9rem;color:#374151;line-height:1.6;">{content}</span>
        </div>""", unsafe_allow_html=True)

    # ── SOP 01 ─────────────────────────────────────────────────────────────
    if sop_choice.startswith("SOP-01"):
        st.subheader("SOP-01: Type IV Collagen Extraction from Human Placenta")
        col1, col2, col3 = st.columns(3)
        col1.metric("Duration", "3–4 days")
        col2.metric("Temperature", "4°C throughout")
        col3.metric("Expected yield", "0.5–2 mg/g tissue")

        st.markdown("**Materials required**")
        st.markdown("""
- Human placenta tissue (fresh or -80°C stored)
- 0.05–0.1 M acetic acid (see Buffer Calculator)
- Pepsin (1 mg/mL in 0.5M acetic acid) — optional for pepsin-soluble collagen
- Dialysis membrane (12–14 kDa MWCO)
- Ultracentrifuge (100,000×g) or high-speed centrifuge (10,000×g)
- Protease inhibitor cocktail (PMSF 1mM, EDTA 5mM)
- NaCl for salting-out (optional purification step)
""")
        st.divider()
        sop_step(1, "Tissue preparation",
                 "Wash fresh placenta 3× with ice-cold PBS. Remove blood vessels and decidua. "
                 "Cut into small pieces (~1cm³). Weigh wet weight. Flash-freeze unused tissue at -80°C.", critical=True)
        sop_step(2, "Decellularisation (optional but recommended)",
                 "Incubate tissue in 1% Triton X-100 / PBS for 24h at 4°C with gentle agitation. "
                 "Wash 5× with PBS. DNase I (50 µg/mL) treatment 2h at 37°C. Verify acellularity by H&E.")
        sop_step(3, "Acid solubilisation",
                 "Homogenise tissue in 0.05 M acetic acid (20 mL per gram wet weight) using a blender or Polytron. "
                 "Stir at 4°C for 24–48h. Add PMSF (1 mM) and EDTA (5 mM) to inhibit proteases.", critical=True)
        sop_step(4, "Centrifugation",
                 "Centrifuge at 10,000×g, 4°C, 60 min. Collect supernatant. "
                 "Re-extract pellet once more with 0.05M acetic acid. Pool supernatants.")
        sop_step(5, "NaCl precipitation (optional — removes fibrinogen)",
                 "Add NaCl to 0.7 M final concentration. Stir 4°C, 1h. "
                 "Centrifuge 10,000×g, 30 min — discard pellet (contains fibrinogen). "
                 "Increase NaCl to 2.5 M to precipitate collagen. Centrifuge, collect pellet.")
        sop_step(6, "Dialysis",
                 "Resuspend pellet in 0.05 M acetic acid. Transfer to dialysis tubing (12–14 kDa MWCO). "
                 "Dialyse against 5L of 0.01 M acetic acid, 3 changes over 48h at 4°C.")
        sop_step(7, "Quantification and storage",
                 "Measure protein concentration by BCA or hydroxyproline assay. "
                 "Verify purity by SDS-PAGE (6–8% gel). Aliquot and store at -80°C or lyophilise.")

        st.markdown('<div class="info-box">💡 <strong>Key QC checks</strong>: SDS-PAGE should show α1 (~185 kDa) and α2 bands. '
                    'No fibrinogen band at ~48 kDa. No β-actin contamination. '
                    'Hydroxyproline:protein ratio should be ~13% for Type IV collagen.</div>',
                    unsafe_allow_html=True)

    # ── SOP 02 ─────────────────────────────────────────────────────────────
    elif sop_choice.startswith("SOP-02"):
        st.subheader("SOP-02: Type I Collagen Extraction from Placenta")
        col1, col2, col3 = st.columns(3)
        col1.metric("Duration", "2–3 days")
        col2.metric("Temperature", "4°C throughout")
        col3.metric("Expected yield", "2–8 mg/g tissue")

        st.markdown("**Materials required**")
        st.markdown("""
- Human placenta tissue
- 0.5 M acetic acid (see Buffer Calculator)
- Pepsin from porcine stomach, 2,500–3,000 U/mg (Sigma P7012)
- 4 M NaCl solution
- Dialysis tubing 12–14 kDa MWCO
- 1 M NaOH (for pepsin inactivation)
""")
        st.divider()
        sop_step(1, "Tissue preparation",
                 "Mince placenta tissue on ice. Remove Wharton's jelly and blood. "
                 "Wash 3× PBS. Dry blot. Weigh (~5–10g per extraction).")
        sop_step(2, "Pepsin digestion",
                 "Prepare pepsin solution: 1 mg/mL pepsin in 0.5M acetic acid. "
                 "Add tissue at 1:10 ratio (w/v). Stir at 4°C for 24–48h. "
                 "Pepsin solubilises non-covalently cross-linked collagen.", critical=True)
        sop_step(3, "Pepsin inactivation",
                 "Raise pH to 7.5–8.0 using 1M NaOH (neutralise acid slowly on ice). "
                 "Stir 30 min at 4°C. This inactivates pepsin (active only at pH <4).", critical=True)
        sop_step(4, "Salting out",
                 "Add NaCl to 2.5 M. Stir 1h at 4°C. Centrifuge 10,000×g, 30 min, 4°C. "
                 "Discard supernatant. Resuspend pellet in 0.5M acetic acid.")
        sop_step(5, "Re-precipitation",
                 "Repeat step 4 twice for higher purity. "
                 "Final pellet = crude Type I collagen.")
        sop_step(6, "Dialysis and storage",
                 "Dialyse against 0.05M acetic acid (3× 5L changes, 48h, 4°C). "
                 "Store at -80°C or lyophilise at -50°C.")

    # ── SOP 03 ─────────────────────────────────────────────────────────────
    elif sop_choice.startswith("SOP-03"):
        st.subheader("SOP-03: Collagen Purification by Dialysis")
        col1, col2 = st.columns(2)
        col1.metric("Duration", "48–72h")
        col2.metric("Temperature", "4°C")

        sop_step(1, "Membrane preparation",
                 "Cut dialysis tubing (12–14 kDa MWCO) to appropriate length. "
                 "Boil in 2% NaHCO₃ + 1mM EDTA for 10 min. Rinse thoroughly with Milli-Q water. "
                 "Store in 0.1% NaHCO₃ at 4°C until use.")
        sop_step(2, "Sample loading",
                 "Clamp one end of tubing. Load collagen solution (max 60% of tubing volume). "
                 "Remove air bubbles. Clamp other end securely.", critical=True)
        sop_step(3, "First dialysis",
                 "Dialyse against 5L of 0.01M acetic acid at 4°C for 8–12h with gentle stirring. "
                 "Use stir bar in the beaker — not in contact with dialysis tube.")
        sop_step(4, "Buffer changes",
                 "Change buffer × 3 total over 48h. Final change: 0.01M acetic acid or Milli-Q water.")
        sop_step(5, "Recovery",
                 "Carefully remove collagen from tubing using a wide-bore pipette. "
                 "Avoid foaming. Transfer to pre-cooled tubes.", critical=True)
        sop_step(6, "Quality check",
                 "SDS-PAGE: check for purity. Conductivity measurement: should drop to <50 µS/cm after dialysis. "
                 "BCA: measure final protein concentration.")

    # ── SOP 04 ─────────────────────────────────────────────────────────────
    elif sop_choice.startswith("SOP-04"):
        st.subheader("SOP-04: SDS-PAGE for Collagen Purity Analysis")
        col1, col2, col3 = st.columns(3)
        col1.metric("Gel %", "6–8%")
        col2.metric("Run time", "~90 min at 100V")
        col3.metric("Detection", "Coomassie / Silver")

        sop_step(1, "Sample preparation",
                 "Mix sample with 4× Laemmli buffer + β-mercaptoethanol (5%). "
                 "Heat 95°C, 5 min. Cool on ice. Spin briefly.", critical=True)
        sop_step(2, "Gel casting",
                 "Cast 7% resolving gel (see Buffer Calculator → SDS-PAGE). Allow to polymerise 30 min. "
                 "Layer 5% stacking gel. Insert comb. Wait 20 min.")
        sop_step(3, "Loading",
                 "Load 5–15 µg protein per lane. Include molecular weight marker (10–250 kDa range). "
                 "Load same volume per lane.")
        sop_step(4, "Running",
                 "Run at 80V through stacking gel, then 100–120V through resolving gel. "
                 "Stop when dye front reaches bottom (≈ 90 min for mini gel).")
        sop_step(5, "Staining",
                 "Fix 30 min in 40% methanol / 10% acetic acid. "
                 "Stain with 0.1% Coomassie R-250 for 1h. "
                 "Destain until bands are clear (multiple changes of 40% MeOH / 10% AcOH).")
        sop_step(6, "Analysis",
                 "Expected bands: α1 (~138 kDa), α2 (~129 kDa) for Type I. "
                 "Type IV: ~185 kDa. β dimers: ~250 kDa. γ trimers: top of gel. "
                 "No band at ~48 kDa (fibrinogen α) or ~94 kDa (fibronectin).")

        st.markdown('<div class="info-box">💡 <strong>Purity criterion</strong>: α1/α2 band intensity ratio ~2:1 for Type I collagen. '
                    'Total collagen bands should be >90% of total protein staining.</div>',
                    unsafe_allow_html=True)

    # ── SOP 05 ─────────────────────────────────────────────────────────────
    elif sop_choice.startswith("SOP-05"):
        st.subheader("SOP-05: Collagen Quantification (BCA + Hydroxyproline)")

        st.markdown("**Method A: BCA Protein Assay (total protein)**")
        sop_step(1, "Standard preparation",
                 "Prepare BSA standards: 0, 25, 125, 250, 500, 750, 1000, 1500, 2000 µg/mL in same buffer as samples.")
        sop_step(2, "Working reagent",
                 "Mix BCA Reagent A : Reagent B = 50:1 (v/v). Prepare fresh. "
                 "Volume: 200 µL per well for 96-well plate.")
        sop_step(3, "Assay",
                 "Add 25 µL sample/standard + 200 µL working reagent per well. "
                 "Incubate 37°C, 30 min. Read absorbance at 562 nm.", critical=True)

        st.markdown("**Method B: Hydroxyproline Assay (collagen-specific)** ⭐ More accurate")
        sop_step(4, "Hydrolysis",
                 "Add 6M HCl (1:1 v/v) to sample. Heat 110°C for 18h in sealed vials. "
                 "Neutralise with 6M NaOH. Adjust to pH 6–7.", critical=True)
        sop_step(5, "Chloramine-T oxidation",
                 "Add chloramine-T reagent (56 mM in acetate-citrate buffer). "
                 "Incubate RT, 20 min.")
        sop_step(6, "Color development",
                 "Add Ehrlich's reagent (DAB in perchloric acid/isopropanol). "
                 "Heat 60°C, 15 min. Read at 550 nm.")
        sop_step(7, "Calculation",
                 "Use trans-4-hydroxyproline standard curve (0–100 µg/mL). "
                 "Collagen content = hydroxyproline × 7.69 (conversion factor, as Hyp is ~13% of collagen by mass).")

        st.markdown('<div class="info-box">💡 <strong>Recommendation</strong>: Use hydroxyproline assay for collagen-specific quantification. '
                    'BCA overestimates if other proteins are present. '
                    'Use both for cross-validation of purity.</div>',
                    unsafe_allow_html=True)

    # ── SOP 06 ─────────────────────────────────────────────────────────────
    elif sop_choice.startswith("SOP-06"):
        st.subheader("SOP-06: Corneal Endothelial Cell (CEC) Culture on Collagen Substrate")
        col1, col2 = st.columns(2)
        col1.metric("Substrate", "Type IV Collagen 2 µg/cm²")
        col2.metric("Medium", "OptiMEM + supplements")

        sop_step(1, "Substrate coating",
                 "Dilute Type IV collagen to 2 µg/cm² in 0.05M acetic acid. "
                 "Add to tissue culture surface. Incubate 1h at 37°C. "
                 "Aspirate excess. Allow to air-dry briefly. Do NOT wash off.")
        sop_step(2, "CEC seeding",
                 "Seed CECs at 2,500–5,000 cells/cm² in complete OptiMEM medium "
                 "(10% FBS, 1% Pen/Strep, 2 mM L-glutamine, 20 µg/mL ascorbic acid). "
                 "Add Y-27632 (10 µM) for first 24h.", critical=True)
        sop_step(3, "Medium change",
                 "Change medium every 2–3 days. CECs form hexagonal monolayer at confluence (~14 days). "
                 "Assess morphology: regular hexagonal shape (similar to in vivo).")
        sop_step(4, "Characterisation",
                 "Stain for ZO-1 (tight junctions), Na⁺/K⁺-ATPase (pump function), CD166, N-cadherin. "
                 "Measure transendothelial electrical resistance (TEER) if applicable.")
        sop_step(5, "Passaging",
                 "Use 0.025–0.05% trypsin + 0.53 mM EDTA, 37°C, 3–5 min. "
                 "Neutralise with medium. Centrifuge 200×g, 5 min. Re-seed on fresh collagen substrate.", critical=True)

        st.markdown('<div class="warn-box">⚠️ CECs are contact-inhibited post-confluence. '
                    'Do not over-grow. Avoid sub-culturing beyond passage 3–4 (phenotype drift risk). '
                    'iPSC-derived CECs may require adjusted conditions.</div>', unsafe_allow_html=True)

    # ── SOP 07 ─────────────────────────────────────────────────────────────
    elif sop_choice.startswith("SOP-07"):
        st.subheader("SOP-07: Collagen Bioink Preparation for 3D Printing")
        col1, col2, col3 = st.columns(3)
        col1.metric("Collagen concentration", "2–4 mg/mL")
        col2.metric("Temperature", "On ice throughout")
        col3.metric("Print temperature", "4°C nozzle / 37°C bed")

        sop_step(1, "Component preparation",
                 "Pre-cool all components and pipettes on ice. "
                 "Prepare: collagen stock solution, 10× PBS, 0.1M NaOH. "
                 "All components must be at 4°C.", critical=True)
        sop_step(2, "Mixing (8:1:1 ratio)",
                 "In order: add collagen stock → 10× PBS → 0.1M NaOH. "
                 "Mix gently by pipetting (avoid bubbles). Do NOT vortex. "
                 "Check pH: target 7.2–7.4 with pH strip or meter.")
        sop_step(3, "Optional: add cells",
                 "For cell-laden bioink: resuspend cells in small volume of medium, "
                 "add to collagen mix as final step. Final density: 1–5 × 10⁶ cells/mL. "
                 "Work quickly — total time on ice <15 min.", critical=True)
        sop_step(4, "Load into cartridge",
                 "Transfer bioink to pre-cooled printing cartridge on ice. "
                 "Centrifuge briefly (200×g, 1 min) to remove bubbles if needed.")
        sop_step(5, "Printing parameters",
                 "Nozzle temperature: 4°C (for extrusion-based). "
                 "Print bed: 37°C (triggers gelation). "
                 "Pressure: 5–15 kPa. Speed: 5–10 mm/s (optimise for your printer).")
        sop_step(6, "Post-printing gelation",
                 "After printing, incubate construct at 37°C for 30–60 min to complete gelation. "
                 "Add culture medium for cell-laden constructs. "
                 "Assess print fidelity by measuring filament diameter vs design.")

    # ── SOP 08 ─────────────────────────────────────────────────────────────
    elif sop_choice.startswith("SOP-08"):
        st.subheader("SOP-08: Fibrinogen Removal from Collagen Preparation")
        st.markdown("Fibrinogen (MW ~340 kDa) is a major contaminant in placenta-derived collagen. "
                    "It appears at ~48 kDa (α chain) on SDS-PAGE.")

        sop_step(1, "NaCl differential precipitation",
                 "Add NaCl to **0.7 M** to crude collagen solution. Stir 1h at 4°C. "
                 "Centrifuge 10,000×g, 30 min, 4°C. "
                 "**Discard PELLET** (fibrinogen precipitates at 0.7M NaCl). "
                 "Retain supernatant.", critical=True)
        sop_step(2, "Collagen precipitation",
                 "Increase NaCl in supernatant to **2.5 M**. Stir 1h at 4°C. "
                 "Centrifuge 10,000×g, 30 min. "
                 "**Retain PELLET** (collagen). Discard supernatant.")
        sop_step(3, "Wash",
                 "Resuspend pellet in 0.5M acetic acid. Repeat NaCl precipitation steps 1–2. "
                 "Repeat 2–3× for high purity.")
        sop_step(4, "Ammonium sulphate fractionation (alternative)",
                 "Alternatively: add ammonium sulphate to 30% saturation → spin → "
                 "fibrinogen precipitates. Collect supernatant. "
                 "Increase to 60% saturation → collagen precipitates.")
        sop_step(5, "Verification by SDS-PAGE",
                 "Run SDS-PAGE (7% gel). No band at ~48 kDa (fibrinogen α chain) = successful removal. "
                 "Western blot with anti-fibrinogen antibody for confirmation.", critical=True)

        st.markdown('<div class="info-box">💡 <strong>Other contaminants to check</strong>:<br/>'
                    '• Fibronectin (~220 kDa): removed by 0.5M NaCl precipitation<br/>'
                    '• Laminin (~400 kDa): requires immunoaffinity chromatography<br/>'
                    '• Metal ions (Zn²⁺, Cu²⁺): remove by EDTA dialysis (5mM EDTA in first dialysis)<br/>'
                    '• DNA: DNase I treatment (50 µg/mL, 37°C, 2h) before acid extraction</div>',
                    unsafe_allow_html=True)

    # Print / save SOP button
    st.divider()
    if st.button("📝 Save this SOP reference to Lab Notes"):
        note = {
            "id": len(st.session_state.notes)+1,
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
            "title": f"SOP Reference — {sop_choice}",
            "text": f"SOP accessed: {sop_choice}\nDate: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
            "type": "Protocol", "tags": ["SOP","reference"], "priority": "Medium",
            "has_file": False, "filename": None,
        }
        st.session_state.notes.insert(0, note)
        save_notes(st.session_state.notes)
        st.success("Saved to Lab Notes!")

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 5 — DATA DASHBOARD
# ══════════════════════════════════════════════════════════════════════════════
elif page == "📊 Data Dashboard":
    st.title("📊 Data Dashboard")

    uploaded_data = st.file_uploader("Upload data (CSV or Excel)",
                                      type=["csv","xlsx"])
    if uploaded_data:
        df = pd.read_csv(uploaded_data) if uploaded_data.name.endswith(".csv") else pd.read_excel(uploaded_data)
        st.success(f"Loaded: {uploaded_data.name} — {df.shape[0]} rows × {df.shape[1]} cols")
        with st.expander("Data preview", expanded=True):
            st.dataframe(df, use_container_width=True)

        num_cols = df.select_dtypes(include=np.number).columns.tolist()
        if len(num_cols) >= 1:
            c1, c2, c3 = st.columns(3)
            with c1: x_col = st.selectbox("X axis", df.columns.tolist())
            with c2: y_col = st.selectbox("Y axis", num_cols)
            with c3: chart = st.selectbox("Chart type", ["Bar","Line","Area","Scatter"])

            if chart == "Bar":    st.bar_chart(df.set_index(x_col)[y_col])
            elif chart == "Line": st.line_chart(df.set_index(x_col)[y_col])
            elif chart == "Area": st.area_chart(df.set_index(x_col)[y_col])
            else:                 st.scatter_chart(df, x=x_col, y=y_col)

            st.subheader("Descriptive statistics")
            st.dataframe(df[num_cols].describe().round(3), use_container_width=True)
    else:
        st.info("No data uploaded. Showing demo — collagen extraction yield across 8 runs.")
        demo = pd.DataFrame({
            "Run":       [f"Run {i}" for i in range(1,9)],
            "Yield_mg":  [42.1,38.7,51.3,48.9,55.2,47.6,60.1,53.8],
            "Purity_%":  [78,82,85,80,88,84,91,87],
            "pH_final":  [3.1,3.0,2.9,3.2,3.0,3.1,2.9,3.0],
        })
        st.dataframe(demo, use_container_width=True)
        c1,c2 = st.columns(2)
        with c1:
            st.markdown("**Yield (mg)**")
            st.bar_chart(demo.set_index("Run")["Yield_mg"])
        with c2:
            st.markdown("**Purity (%)**")
            st.line_chart(demo.set_index("Run")["Purity_%"])

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 6 — SETTINGS
# ══════════════════════════════════════════════════════════════════════════════
elif page == "⚙️ Settings":
    st.title("⚙️ Settings")
    st.subheader("Researcher Profile")
    c1,c2 = st.columns(2)
    with c1:
        st.text_input("Name", value="Jiyoung Moon")
        st.text_input("Institution", value="University of Sydney")
    with c2:
        st.text_input("Lab", value="Corneal Bioengineering Lab")
        st.text_input("Project", value="Type IV Collagen – Corneal Blindness Treatment")

    st.subheader("Data Management")
    c1,c2 = st.columns(2)
    with c1:
        if st.button("🗑️ Clear all notes"):
            st.session_state.notes = []
            save_notes([])
            st.success("Cleared.")
    with c2:
        if st.session_state.notes:
            st.download_button("⬇️ Export notes",
                                pd.DataFrame(st.session_state.notes).to_csv(index=False).encode(),
                                file_name="all_lab_notes.csv", mime="text/csv")

    st.subheader("Version")
    st.markdown("""
**CollagenLab Notebook v2.0**

New in v2:
- 10 detailed buffer calculators (Acetic acid, PBS, HEPES, Dialysis, Neutralization, Collagenase, Trypsin, SDS-PAGE, BCA, ELISA)
- 8 full SOPs (Type I/IV extraction, purification, SDS-PAGE, quantification, CEC culture, bioink, fibrinogen removal)
- Protocol tips and critical step warnings for each SOP
- Save buffer calculations and SOP references to Lab Notes

Roadmap v3:
- ML yield prediction
- Chromatography peak analysis
- Cloud sync
- Mobile PWA
""")
