"""
AI structure advisor: identifies what a receptor structure actually IS
(protein family, function, bound cofactors/ions, relevant literature) and asks
an LLM to turn that into concrete, protein-aware receptor-prep guidance -
which HETATM groups to keep vs strip, expected protonation-sensitive
residues, metal-site handling, known gotchas from the literature, etc.

Two backends are supported, selectable with --backend:
  - "anthropic" (default): Claude via the Anthropic API. Needs ANTHROPIC_API_KEY.
    New Console accounts get a small one-time free trial credit; ongoing use
    is pay-as-you-go. See https://console.anthropic.com
  - "gemini": Google's Gemini API, which has an ongoing free tier (rate-limited,
    no credit card) for the Flash/Flash-Lite models. Needs GEMINI_API_KEY, get
    one free at https://aistudio.google.com/apikey. Note: on Gemini's free
    tier, Google's terms allow using your inputs/outputs to improve their
    models - fine for public PDB structures, worth knowing regardless.

DESIGN PRINCIPLE - grounded, not vibes-based:
The LLM is never asked to "recall" a PDB structure from memory alone (that
invites hallucination on anything less famous than lysozyme). Instead this
module first pulls hard facts:
  - RCSB entry metadata (title, keywords, polymer descriptions, organism,
    experimental method/resolution, primary citation + PubMed ID)
  - The primary citation's PubMed abstract, if available
  - The ACTUAL HETATM groups present in the user's own structure file,
    parsed directly off the coordinates (not guessed)
...and only then asks the LLM to reason over those specific facts and its
general biochemistry knowledge. The report always shows what was actually
fetched, so you can sanity-check the AI's reasoning against the same sources
it used.

THIS IS ADVISORY ONLY. See LIMITATIONS at the bottom of this file (also
printed with --help) before trusting it unsupervised. Nothing in this module
mutates receptor.pdb or receptor.pdbqt automatically - it only produces a
report for a human to read before preparing the structure.

Requires (only when --ai-advise is actually used):
  - Internet access
  - ANTHROPIC_API_KEY (default backend) or GEMINI_API_KEY (--backend gemini)

Standalone use:
    python structure_advisor.py --pdb-id 8FK4 --output-dir ./prepared
    python structure_advisor.py --pdb-file receptor.pdb --pdb-id 8FK4 --output-dir ./prepared --backend gemini

Library use (see prepare_structures.py for the integration point):
    from structure_advisor import run_advisor
    result = run_advisor(pdb_id="8FK4", pdb_path="receptor.pdb", backend="gemini")
"""

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request

RCSB_ENTRY_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{id}"
PUBMED_ESUMMARY_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    "?db=pubmed&retmode=json&id={pmid}"
)
PUBMED_EFETCH_ABSTRACT_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=pubmed&rettype=abstract&retmode=text&id={pmid}"
)

ANTHROPIC_API_URL = "https://api.anthropic.com/v1/messages"
ANTHROPIC_VERSION = "2023-06-01"
# Model string as of this writing - Anthropic updates these over time, so if
# this starts failing with a "model not found"-style error, check
# https://docs.claude.com/en/docs/about-claude/models/overview and update it
# (or just pass --model explicitly).
DEFAULT_ANTHROPIC_MODEL = "claude-sonnet-5"

GEMINI_API_URL_TEMPLATE = "https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent"
# gemini-3.5-flash is (as of this writing, July 2026) Google's current
# general-availability Flash model and is inside Gemini's free tier - a
# rate-limited but ongoing no-credit-card allowance, unlike Anthropic's
# one-time trial credit. The Gemini 2.x line (including the older
# gemini-2.5-flash this used to default to) has been retired for new users.
# Free-tier model availability changes over time - check
# https://ai.google.dev/gemini-api/docs/pricing if this stops working,
# or pass --model explicitly (e.g. a Flash-Lite variant for higher rate limits).
DEFAULT_GEMINI_MODEL = "gemini-3.5-flash"

DEFAULT_MODEL_BY_BACKEND = {"anthropic": DEFAULT_ANTHROPIC_MODEL, "gemini": DEFAULT_GEMINI_MODEL}
API_KEY_ENV_BY_BACKEND = {"anthropic": "ANTHROPIC_API_KEY", "gemini": "GEMINI_API_KEY"}

REQUEST_TIMEOUT = 30

WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD"}
ION_RESNAMES = {
    "NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO",
    "CD", "HG", "LI", "CS", "BA", "SR", "AG", "AU", "PT", "AL", "GA", "IN", "PB", "SN", "TL",
}
CRYO_BUFFER_RESNAMES = {
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "DMS", "ACT", "FMT",
    "TRS", "BME", "MPD", "IPA", "MOH", "EOH", "CO3", "NO3", "UNK", "UNX", "UNL", "IOD",
    "BR", "AZI", "ACY", "CIT", "TAR", "MRD", "BOG", "LDA", "P6G", "P4G", "OGA", "BU3",
    "SIN", "MES", "IMD", "EPE", "BEZ", "PLM", "OLA", "MYR", "PLC", "GSH", "NH4", "OXL",
}


class AdvisorError(Exception):
    """Raised for advisor-specific failures (network, API, parsing) so callers
    can catch this one type and fall back to plain mechanical prep."""


# ---------------------------------------------------------------------------
# Step 1: hard facts from RCSB / PubMed
# ---------------------------------------------------------------------------

def _http_get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "structure-advisor/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
            return json.loads(resp.read().decode("utf-8", errors="replace"))
    except urllib.error.HTTPError as e:
        raise AdvisorError(f"HTTP {e.code} fetching {url}")
    except urllib.error.URLError as e:
        raise AdvisorError(f"network error fetching {url}: {e.reason}")


def _http_get_text(url):
    req = urllib.request.Request(url, headers={"User-Agent": "structure-advisor/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        raise AdvisorError(f"HTTP {e.code} fetching {url}")
    except urllib.error.URLError as e:
        raise AdvisorError(f"network error fetching {url}: {e.reason}")


def fetch_rcsb_entry(pdb_id):
    """Pull title/keywords/method/citation metadata for a PDB ID from RCSB's
    public data API. No API key needed. Returns a plain dict of the fields we
    care about, or raises AdvisorError."""
    pdb_id = pdb_id.strip().upper()
    if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
        raise AdvisorError(f"'{pdb_id}' doesn't look like a 4-character PDB ID")

    data = _http_get_json(RCSB_ENTRY_URL_TEMPLATE.format(id=pdb_id))

    struct = data.get("struct", {}) or {}
    exptl = data.get("exptl", []) or []
    citations = data.get("rcsb_primary_citation") or {}
    keywords = (data.get("struct_keywords") or {}).get("pdbx_keywords", "")
    entity_descriptions = []
    for poly in (data.get("rcsb_entry_container_identifiers", {}) or {}).get(
        "polymer_entity_ids", []
    ) or []:
        entity_descriptions.append(poly)  # placeholder ids; full names need a second call, skipped for brevity

    pubmed_id = citations.get("pdbx_database_id_pub_med")

    return {
        "pdb_id": pdb_id,
        "title": struct.get("title", ""),
        "keywords": keywords,
        "experimental_method": exptl[0].get("method") if exptl else None,
        "resolution": (data.get("rcsb_entry_info", {}) or {}).get(
            "resolution_combined", [None]
        )[0],
        "citation_title": citations.get("title"),
        "citation_journal": citations.get("rcsb_journal_abbrev"),
        "citation_year": citations.get("year"),
        "pubmed_id": str(pubmed_id) if pubmed_id else None,
    }


def fetch_pubmed_abstract(pubmed_id, max_chars=2500):
    """Best-effort abstract fetch. Returns None (never raises) if anything
    goes wrong - literature context is a nice-to-have, not load-bearing."""
    try:
        text = _http_get_text(PUBMED_EFETCH_ABSTRACT_URL.format(pmid=pubmed_id))
        text = text.strip()
        if not text or "error" in text.lower()[:200]:
            return None
        return text[:max_chars]
    except AdvisorError:
        return None


# ---------------------------------------------------------------------------
# Step 2: hard facts from the user's OWN structure file (not guessed)
# ---------------------------------------------------------------------------

def summarize_hetero_groups(pdb_path):
    """Parse HETATM groups directly out of the local file so the LLM reasons
    about what is ACTUALLY in this structure, not a generic guess. Mirrors the
    classification logic in redock_and_validate.py."""
    if not os.path.isfile(pdb_path):
        raise AdvisorError(f"no such file: {pdb_path}")

    groups = {}
    with open(pdb_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith("HETATM") or len(line) < 20:
                continue
            res_name = line[17:20].strip()
            chain_id = line[21].strip() or " "
            res_seq = line[22:26].strip()
            elem = line[76:78].strip().upper() if len(line) >= 78 else ""
            if elem == "H":
                continue
            key = (res_name, chain_id, res_seq)
            groups.setdefault(key, 0)
            groups[key] += 1

    summary = {"waters": 0, "ions": [], "buffers": [], "other_hetero": []}
    for (res_name, chain_id, res_seq), n_atoms in groups.items():
        if res_name in WATER_RESNAMES:
            summary["waters"] += 1
        elif res_name in ION_RESNAMES:
            summary["ions"].append(f"{res_name} (chain {chain_id}, resSeq {res_seq})")
        elif res_name in CRYO_BUFFER_RESNAMES:
            summary["buffers"].append(f"{res_name} (chain {chain_id}, resSeq {res_seq})")
        else:
            summary["other_hetero"].append(
                f"{res_name} (chain {chain_id}, resSeq {res_seq}, {n_atoms} heavy atoms)"
            )
    return summary


# ---------------------------------------------------------------------------
# Step 3: build a grounded prompt and call Claude
# ---------------------------------------------------------------------------

_SYSTEM_PROMPT = """You are assisting a computational chemist in preparing a protein \
structure for AutoDock Vina docking. You will be given hard facts pulled directly from \
RCSB PDB, PubMed, and the user's own structure file - not your own recollection of the \
structure. Reason ONLY from the facts provided plus your general biochemistry/structural \
biology knowledge of this protein family. If you are not confident about something \
specific to THIS structure, say so explicitly rather than inventing detail - a wrong \
recommendation here can quietly corrupt a docking study.

Respond with ONLY a single JSON object (no markdown fences, no commentary before or \
after) matching exactly this schema:

{
  "protein_family": "short description, e.g. 'Serine/threonine protein kinase, CMGC group'",
  "function_summary": "2-3 sentences on biological function relevant to structure prep",
  "identification_confidence": "high" | "medium" | "low",
  "hetero_group_guidance": [
    {"group": "e.g. ZN chain A resSeq 301", "recommendation": "keep as structural cofactor" | "keep as catalytic cofactor" | "strip as crystallization additive" | "review manually", "rationale": "one sentence"}
  ],
  "protonation_notes": "notable residues whose protonation state matters for this protein family/mechanism and why, or 'nothing unusual' if so",
  "water_notes": "guidance on structural/catalytic waters worth keeping vs bulk solvent, framed generally since exact water IDs aren't in the facts given",
  "flexibility_warnings": ["short warnings about known flexible loops, induced-fit behavior, alternate conformations etc. for this protein family, if relevant - empty list if none"],
  "literature_notes": "1-3 sentences summarizing anything prep-relevant from the citation/abstract given, or 'no relevant literature details available' if none was provided",
  "caveats": ["explicit list of things this advisory could NOT verify and that the user should check manually"],
  "overall_confidence": "high" | "medium" | "low"
}
"""


def build_user_prompt(rcsb_meta, abstract_text, hetero_summary, pdb_id):
    lines = []
    lines.append(f"PDB ID: {pdb_id or 'unknown / locally-supplied file, no RCSB lookup'}")
    if rcsb_meta:
        lines.append(f"RCSB title: {rcsb_meta.get('title') or 'n/a'}")
        lines.append(f"RCSB keywords: {rcsb_meta.get('keywords') or 'n/a'}")
        lines.append(f"Experimental method: {rcsb_meta.get('experimental_method') or 'n/a'}")
        lines.append(f"Resolution: {rcsb_meta.get('resolution') or 'n/a'}")
        if rcsb_meta.get("citation_title"):
            lines.append(
                f"Primary citation: \"{rcsb_meta['citation_title']}\" "
                f"({rcsb_meta.get('citation_journal') or 'n/a'}, {rcsb_meta.get('citation_year') or 'n/a'}), "
                f"PubMed ID {rcsb_meta.get('pubmed_id') or 'n/a'}"
            )
    else:
        lines.append("(No RCSB metadata available - proceeding on structure-file facts only.)")

    if abstract_text:
        lines.append("\nPrimary citation abstract (from PubMed):")
        lines.append(abstract_text)
    else:
        lines.append("\n(No abstract available.)")

    lines.append("\nHETATM groups actually present in the user's structure file:")
    lines.append(f"  Water molecules: {hetero_summary['waters']} (resnames HOH/WAT/etc.)")
    lines.append(f"  Ions: {', '.join(hetero_summary['ions']) or 'none detected'}")
    lines.append(
        f"  Likely crystallization buffer/cryoprotectant components: "
        f"{', '.join(hetero_summary['buffers']) or 'none detected'}"
    )
    lines.append(
        f"  Other HETATM groups (candidate ligands/cofactors): "
        f"{', '.join(hetero_summary['other_hetero']) or 'none detected'}"
    )
    lines.append(
        "\nGive prep guidance for docking with AutoDock Vina, per the JSON schema in "
        "the system prompt."
    )
    return "\n".join(lines)


def call_claude(system_prompt, user_prompt, api_key, model):
    """Anthropic API backend."""
    body = json.dumps({
        "model": model,
        "max_tokens": 1500,
        "system": system_prompt,
        "messages": [{"role": "user", "content": user_prompt}],
    }).encode("utf-8")

    req = urllib.request.Request(
        ANTHROPIC_API_URL,
        data=body,
        method="POST",
        headers={
            "Content-Type": "application/json",
            "x-api-key": api_key,
            "anthropic-version": ANTHROPIC_VERSION,
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
    except urllib.error.HTTPError as e:
        detail = e.read().decode("utf-8", errors="replace")[:500]
        raise AdvisorError(f"Anthropic API returned HTTP {e.code}: {detail}")
    except urllib.error.URLError as e:
        raise AdvisorError(f"network error calling Anthropic API: {e.reason}")

    text_parts = [block["text"] for block in data.get("content", []) if block.get("type") == "text"]
    raw_text = "".join(text_parts).strip()
    if not raw_text:
        raise AdvisorError(f"empty response from Anthropic API: {data}")
    return _parse_json_response(raw_text)


def call_gemini(system_prompt, user_prompt, api_key, model):
    """Gemini API backend (generativelanguage.googleapis.com) - free tier for
    Flash/Flash-Lite models, rate-limited, no credit card required."""
    body = json.dumps({
        "contents": [{"role": "user", "parts": [{"text": user_prompt}]}],
        "systemInstruction": {"parts": [{"text": system_prompt}]},
        "generationConfig": {"responseMimeType": "application/json"},
    }).encode("utf-8")

    url = GEMINI_API_URL_TEMPLATE.format(model=model)
    req = urllib.request.Request(
        url,
        data=body,
        method="POST",
        headers={"Content-Type": "application/json", "x-goog-api-key": api_key},
    )
    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
    except urllib.error.HTTPError as e:
        detail = e.read().decode("utf-8", errors="replace")[:500]
        raise AdvisorError(f"Gemini API returned HTTP {e.code}: {detail}")
    except urllib.error.URLError as e:
        raise AdvisorError(f"network error calling Gemini API: {e.reason}")

    candidates = data.get("candidates") or []
    if not candidates:
        feedback = data.get("promptFeedback")
        raise AdvisorError(f"no candidates in Gemini response (blocked or empty). promptFeedback={feedback}")
    finish_reason = candidates[0].get("finishReason")
    parts = (candidates[0].get("content") or {}).get("parts") or []
    raw_text = "".join(p.get("text", "") for p in parts).strip()
    if not raw_text:
        raise AdvisorError(f"empty response from Gemini API (finishReason={finish_reason}): {data}")
    return _parse_json_response(raw_text)


def _parse_json_response(raw_text):
    """Shared JSON extraction: strips markdown fences defensively even though
    both backends are asked for raw JSON, since models sometimes add them anyway."""
    cleaned = re.sub(r"^```(?:json)?\s*|\s*```$", "", raw_text.strip(), flags=re.MULTILINE).strip()
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError as e:
        raise AdvisorError(f"could not parse JSON from model response ({e}). Raw response:\n{raw_text}")


_CALL_FN_BY_BACKEND = {"anthropic": call_claude, "gemini": call_gemini}


def call_llm(system_prompt, user_prompt, api_key, model, backend):
    if backend not in _CALL_FN_BY_BACKEND:
        raise AdvisorError(f"unknown backend '{backend}' - must be one of {list(_CALL_FN_BY_BACKEND)}")
    return _CALL_FN_BY_BACKEND[backend](system_prompt, user_prompt, api_key, model)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def run_advisor(pdb_id=None, pdb_path=None, api_key=None, model=None, backend="anthropic"):
    """Runs the full advisory pipeline. Returns a dict:
        {"status": "ok", "recommendation": {...}, "sources": {...}}
        {"status": "skipped", "reason": "..."}
        {"status": "error", "reason": "..."}
    Never raises - callers (e.g. prepare_structures.py) can treat this as a
    pure side-channel report and keep going with mechanical prep regardless.

    backend: "anthropic" (default, needs ANTHROPIC_API_KEY) or "gemini"
    (needs GEMINI_API_KEY - has an ongoing free tier, see module docstring).
    """
    if backend not in DEFAULT_MODEL_BY_BACKEND:
        return {"status": "error", "reason": f"unknown backend '{backend}' - use 'anthropic' or 'gemini'"}

    model = model or DEFAULT_MODEL_BY_BACKEND[backend]
    api_key = api_key or os.environ.get(API_KEY_ENV_BY_BACKEND[backend])
    if not api_key:
        return {
            "status": "skipped",
            "reason": f"no {API_KEY_ENV_BY_BACKEND[backend]} set and no --api-key given "
                      f"(backend={backend})",
        }
    if not pdb_path and not pdb_id:
        return {"status": "skipped", "reason": "no --pdb-id or --pdb-file given"}

    rcsb_meta = None
    abstract_text = None
    if pdb_id:
        try:
            rcsb_meta = fetch_rcsb_entry(pdb_id)
            if rcsb_meta.get("pubmed_id"):
                abstract_text = fetch_pubmed_abstract(rcsb_meta["pubmed_id"])
        except AdvisorError as e:
            print(f"  NOTE: RCSB metadata lookup failed ({e}); continuing with structure-file facts only.")

    hetero_summary = {"waters": 0, "ions": [], "buffers": [], "other_hetero": []}
    if pdb_path:
        try:
            hetero_summary = summarize_hetero_groups(pdb_path)
        except AdvisorError as e:
            print(f"  NOTE: could not parse HETATM groups from {pdb_path} ({e}).")

    user_prompt = build_user_prompt(rcsb_meta, abstract_text, hetero_summary, pdb_id)
    try:
        recommendation = call_llm(_SYSTEM_PROMPT, user_prompt, api_key, model, backend)
    except AdvisorError as e:
        return {"status": "error", "reason": str(e)}

    return {
        "status": "ok",
        "recommendation": recommendation,
        "sources": {
            "rcsb_meta": rcsb_meta,
            "abstract_used": bool(abstract_text),
            "hetero_summary": hetero_summary,
            "model": model,
            "backend": backend,
        },
    }


def format_report(result):
    if result["status"] == "skipped":
        return f"AI advisory skipped: {result['reason']}"
    if result["status"] == "error":
        return f"AI advisory FAILED: {result['reason']}"

    rec = result["recommendation"]
    src = result["sources"]
    lines = []
    lines.append("=" * 60)
    lines.append("AI STRUCTURE ADVISORY (informational - not auto-applied)")
    lines.append("=" * 60)
    lines.append(f"Backend: {src.get('backend')}   |   Model: {src.get('model')}   |   Abstract used: {src.get('abstract_used')}")
    lines.append("")
    lines.append(f"Protein family:  {rec.get('protein_family', 'n/a')}")
    lines.append(f"Function:        {rec.get('function_summary', 'n/a')}")
    lines.append(f"ID confidence:   {rec.get('identification_confidence', 'n/a')}")
    lines.append("")
    lines.append("Hetero group guidance (compare against what's actually in your file):")
    for g in rec.get("hetero_group_guidance", []):
        lines.append(f"  - {g.get('group', '?')}: {g.get('recommendation', '?')} - {g.get('rationale', '')}")
    if not rec.get("hetero_group_guidance"):
        lines.append("  (none)")
    lines.append("")
    lines.append(f"Protonation notes: {rec.get('protonation_notes', 'n/a')}")
    lines.append(f"Water notes:       {rec.get('water_notes', 'n/a')}")
    lines.append("")
    warnings = rec.get("flexibility_warnings") or []
    lines.append("Flexibility warnings:")
    for w in warnings:
        lines.append(f"  - {w}")
    if not warnings:
        lines.append("  (none)")
    lines.append("")
    lines.append(f"Literature notes: {rec.get('literature_notes', 'n/a')}")
    lines.append("")
    caveats = rec.get("caveats") or []
    lines.append("Caveats / things NOT verified by this advisory:")
    for c in caveats:
        lines.append(f"  - {c}")
    if not caveats:
        lines.append("  (none listed by the model - treat that as a gap, not a clean bill of health)")
    lines.append("")
    lines.append(f"Overall confidence: {rec.get('overall_confidence', 'n/a')}")
    lines.append("=" * 60)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="AI structure identification + literature-informed receptor prep advisory (report-only)."
    )
    parser.add_argument("--pdb-id", default=None, help="4-character PDB ID to fetch RCSB/PubMed context for")
    parser.add_argument("--pdb-file", default=None, help="Local .pdb file to scan for actual HETATM groups")
    parser.add_argument("--output-dir", default=".", help="Where to save the advisory report")
    parser.add_argument("--backend", choices=["anthropic", "gemini"], default="anthropic",
                         help="Which API to use. 'anthropic' (default) needs ANTHROPIC_API_KEY. "
                              "'gemini' needs GEMINI_API_KEY and has an ongoing free tier - "
                              "get a key free at https://aistudio.google.com/apikey")
    parser.add_argument("--model", default=None,
                         help="Model string override (default depends on --backend: "
                              f"{DEFAULT_ANTHROPIC_MODEL} for anthropic, {DEFAULT_GEMINI_MODEL} for gemini)")
    parser.add_argument("--api-key", default=None,
                         help="API key (default: ANTHROPIC_API_KEY or GEMINI_API_KEY env var, per --backend)")
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.pdb_id and not args.pdb_file:
        print("ERROR: provide at least one of --pdb-id or --pdb-file.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    result = run_advisor(
        pdb_id=args.pdb_id, pdb_path=args.pdb_file, api_key=args.api_key,
        model=args.model, backend=args.backend,
    )
    print(format_report(result))

    if result["status"] == "ok":
        stem = args.pdb_id or os.path.splitext(os.path.basename(args.pdb_file))[0]
        json_path = os.path.join(args.output_dir, f"{stem}_ai_advisory.json")
        md_path = os.path.join(args.output_dir, f"{stem}_ai_advisory.md")
        with open(json_path, "w") as f:
            json.dump(result, f, indent=2)
        with open(md_path, "w") as f:
            f.write(format_report(result))
        print(f"\nSaved: {json_path}")
        print(f"Saved: {md_path}")
    elif result["status"] == "error":
        sys.exit(1)


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# LIMITATIONS (see also the chat message this shipped with)
# ---------------------------------------------------------------------------
# 1. Report-only. Nothing here writes to receptor.pdb/.pdbqt automatically -
#    by design, so a bad AI call can't silently corrupt a prep run. Wiring
#    its recommendations into actual prepare_receptor/obabel flags is left
#    for you to do deliberately once you trust specific recommendation types.
# 2. Hetero-group guidance is matched to your file's ACTUAL resnames, but the
#    model never sees exact 3D coordinates/geometry - e.g. it can't tell a
#    genuinely catalytic Zn from a crystallization artifact Zn purely by
#    resname; it's reasoning from typical biology for that protein family.
# 3. Literature grounding is limited to the single RCSB primary citation's
#    PubMed abstract (if present) - not a real literature search. No
#    PubMed/Google Scholar search across multiple papers about the target.
# 4. RCSB polymer entity descriptions (chain-level names) are not fully
#    pulled in this version - only entry-level title/keywords - so multi-
#    chain hetero-complexes may get a shallower identification than a single
#    well-annotated monomer.
# 5. No sequence-based identification (BLAST/UniProt) - relies entirely on
#    PDB metadata, so a stripped/renamed local file with no matching PDB ID
#    gets identified from HETATM content alone, which is much weaker.
# 6. No caching/rate-limit handling for RCSB, PubMed, or the Anthropic API -
#    fine for one-off interactive use, not built for batch-processing a
#    large structure library without adding throttling yourself.
# 7. LLM output can still be wrong or overconfident despite the grounding
#    effort - always treat this as a second opinion from someone who has
#    read a lot of structural biology, not a validated computational result.
# 8. Requires internet access + an API key on the machine actually running
#    this script - ANTHROPIC_API_KEY (default backend) or GEMINI_API_KEY
#    (--backend gemini). Won't work inside network-isolated environments.
# 9. The two backends are not guaranteed to give equivalent-quality advice -
#    they're different model families with different biochemistry knowledge
#    depth. Gemini's free tier is also rate-limited (requests/minute and/or
#    /day caps that change over time) - fine for occasional interactive use,
#    not for batch-processing many structures back to back.
