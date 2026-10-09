# AI Structure Advisor — what it does, and its limitations

## What was added

**`structure_advisor.py`** (new file) — identifies a receptor structure and produces
protein-aware prep guidance, grounded in real data:

1. **RCSB PDB metadata** (title, keywords, experimental method/resolution, primary
   citation + PubMed ID) — via RCSB's public REST API, no key needed.
2. **PubMed abstract** of that primary citation, if one exists.
3. **The actual HETATM groups in your local structure file** — parsed directly off
   the coordinates (waters / ions / crystallization buffers / other groups), not guessed.
4. All of the above is handed to Claude (Anthropic API) with instructions to reason
   only from the given facts + general biochemistry knowledge, and to flag anything
   it's not confident about — returned as structured JSON (protein family, function,
   which hetero groups to keep/strip and why, protonation notes, water notes,
   flexibility warnings, literature notes, explicit caveats, confidence level).

**`prepare_structures.py`** (modified) — new flags:
```
--ai-advise               turn the advisory step on
--pdb-id 8FK4              give it a real PDB ID for full RCSB/PubMed context
--backend anthropic|gemini  which API to use (default: anthropic)
--model <model-string>     override the default model for the chosen backend
--api-key <key>            override the ANTHROPIC_API_KEY/GEMINI_API_KEY env var
```
Run: `python prepare_structures.py --receptor receptor.pdb --output-dir ./prepared --ai-advise --pdb-id 8FK4`

**Two backends, pick with `--backend`:**
- `anthropic` (default) — Claude via the Anthropic API. New Console accounts get a
  small one-time free trial credit (no card, phone verification); ongoing use is
  pay-as-you-go. Get a key at https://console.anthropic.com
- `gemini` — Google's Gemini API, which has an *ongoing* free tier (rate-limited,
  no credit card, doesn't run out like a trial credit) for the Flash/Flash-Lite
  models. Get a key free at https://aistudio.google.com/apikey. Set it as
  `GEMINI_API_KEY`. Note: on the free tier, Google's terms allow using your
  inputs/outputs to improve their models — a non-issue for public PDB structures,
  worth knowing if you ever point this at anything proprietary.

This writes `<name>_ai_advisory.md` and `.json` next to the prepared `.pdbqt`.
**It is report-only.** Nothing about the advisory step touches `receptor.pdb` or
`receptor.pdbqt` — those still go through the exact same mechanical
`prepare_receptor`/`obabel` pipeline as before. If the advisory step fails for any
reason (no key, no internet, API error), it prints a note and the mechanical prep
still completes normally — it can't break your existing workflow.

I tested the full plumbing end-to-end (arg parsing → RCSB lookup → PubMed →
prompt-building → JSON parsing → report writing) for both backends, using a
mocked LLM response, and confirmed the graceful-failure paths for no-API-key,
RCSB-unreachable, and API errors. I could not do a live end-to-end run against
either real API from this sandbox (its network egress is locked to a small
allowlist that includes api.anthropic.com but not generativelanguage.googleapis.com,
so the Gemini call was blocked at my sandbox's proxy before it even reached
Google — not a code issue, just an environment restriction on my end) — **please
do one live smoke test on your machine with a real key before relying on
either backend.**

## Limitations (so you know what to build on)

1. **Report-only by design.** It never auto-edits `receptor.pdb`/`.pdbqt`. Turning
   specific recommendation types (e.g. "strip this buffer HETATM", "use pH X for
   protonation") into actual automated flags for `prepare_receptor`/`obabel` is
   deliberately left undone — that's a real chemistry decision and I didn't want a
   bad AI call silently corrupting a prep run. This is the main thing to build next.
2. **No 3D geometry awareness.** The model sees resnames/counts from your HETATM
   groups, not actual 3D positions. It can't distinguish a genuinely catalytic Zn
   from a crystallization-artifact Zn by geometry — it's reasoning from what's
   typical for that protein family, not from your structure's specific geometry.
3. **Literature grounding is thin.** Only the single RCSB primary-citation abstract,
   if one exists. This is not a real literature search across multiple papers.
4. **Shallow RCSB metadata.** Only entry-level title/keywords are pulled, not full
   chain-level polymer entity descriptions — multi-chain hetero-complexes will get
   a shallower identification than a well-annotated single-chain entry.
5. **No sequence-based ID.** No BLAST/UniProt lookup — if you give it a local file
   with no matching PDB ID, identification relies solely on HETATM content, which
   is a much weaker signal.
6. **No batching/rate-limiting.** Fine for interactive one-off use; if you point
   this at a large structure library you'll want to add throttling/caching yourself.
7. **Still a second opinion, not ground truth.** LLMs can be confidently wrong,
   especially on obscure or poorly-annotated structures. The report always shows
   its sources so you can sanity-check it — treat it as a knowledgeable colleague's
   quick read, not a validated computational result.
8. **Needs internet + `ANTHROPIC_API_KEY`** on whatever machine actually runs it.
9. **Not wired into `redock_and_validate.py` or the GUI (`redock_gui.py`).** Only
   `prepare_structures.py` calls it right now. `structure_advisor.run_advisor()` is
   a plain importable function, so wiring it into the other scripts/GUI later is
   straightforward if you want it there too.
10. **Model strings may go stale.** Defaults to `claude-sonnet-5` (anthropic) /
    `gemini-3.5-flash` (gemini, as of July 2026 - the 2.x Gemini line this
    used to default to has since been retired for new users); if calls start
    failing with a "model not found" style error, check the provider's
    current model list (linked in `structure_advisor.py`) or pass `--model`.
11. **The two backends aren't guaranteed equivalent.** Different model
    families, different depth of biochemistry knowledge — treat a switch
    between them as a real change in who you're asking, not an interchangeable
    toggle. Gemini's free tier is also rate-limited (requests/minute and/or
    /day caps that change over time), fine for occasional use, not batching.
