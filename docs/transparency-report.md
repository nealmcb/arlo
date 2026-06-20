# Arlo Comparison Audit Transparency: Report and Guide

## Executive Summary

Arlo supports sufficient transparency for independently verifiable comparison audits: all data needed for verification can be exported by election officials in a timely manner at each stage. Some exports are already convenient and well-structured. The main gap is not missing data, but missing guidance — Arlo does not yet steer jurisdictions through the best-practice workflow of exporting, hashing, signing, and publishing each phase's artifacts before moving to the next phase.

Key design points that are correct as-is:
- **Arlo does not produce auditor transcripts** (pre-filled records of what each sampled ballot is expected to show). This is correct: producing such transcripts inside the audit tool would invite violating the blind-audit principle and could be used to script human interpretations. The right approach is for a *separate* observer tool to join the published sample list with the published CVR file *after* both are public, generating verification materials only for observers — never for audit boards.
- **Arlo runs on a jurisdiction's internal, firewalled network**, not as a public-facing service. Election officials are responsible for exporting audit artifacts and publishing them to a public site or repository. Adding public API endpoints to Arlo itself would expand its attack surface and is inappropriate for an air-gap-friendly audit tool.

Concrete improvements — ordered by priority — are described at the end of this document.

---

## Stage 1: Before the Random Seed — Uploading Inputs

### What Arlo Can Export Today

**For ballot comparison audits**, the audit admin can download (via authenticated API endpoints — Arlo runs on a firewalled internal network, so these are only accessible to authorized officials):
- **Jurisdictions file** (`GET /election/<id>/jurisdictions/file`) — the list of participating counties/jurisdictions
- **Standardized contests file** (`GET /election/<id>/standardized-contests/file`) — the set of contests to be audited
- **Ballot manifest per jurisdiction** (`GET /election/<id>/jurisdiction/<id>/ballot-manifest/csv`) — batch/tabulator structure and ballot counts
- **CVR file per jurisdiction** (`GET /election/<id>/jurisdiction/<id>/cvrs/csv`) — the cast vote record for every ballot, with imprinted IDs. Note that many jurisdictions are required to redact CVRs before publication to prevent vote revelation from rare ballot styles (see [Bernhard et al., 2025](https://www.science.org/doi/10.1126/sciadv.adt1512)); Arlo's CVR export does not currently apply any redaction.

**For batch comparison audits**, additionally:
- **Batch tallies per jurisdiction** (`GET /election/<id>/jurisdiction/<id>/batch-tallies/csv`) — reported vote counts per batch
- **Batch inventory worksheet, CVR, and ballot manifest** (generated from batch inventory sign-off flow)
- **Bundle downloads**: Arlo recently added `POST /election/<id>/batch-files/manifests-bundle` and `POST /election/<id>/batch-files/candidate-totals-bundle` — these assemble all jurisdiction manifests or batch tallies into a single ZIP, compute a **SHA-256 hash** of the inner archive, and include a `*-sha256-hash.txt` file in the outer ZIP. This is the most significant transparency artifact for the pre-seed stage.

**Audit settings** (risk limit, math type, contest definitions with vote totals) are entered before the seed and appear in the final audit report.

### Gaps

- **No equivalent bundle for CVR files** in ballot comparison audits — the most critical input for ballot comparison audits lacks the batch-bundle treatment given to batch comparison manifests.
- **No per-jurisdiction ZIP for each phase** — most observation happens at the county level. Providing a per-jurisdiction ZIP (manifest + CVR + settings) at each stage (pre-seed, pre-round, final) would be more practical than a single ZIP of all CVRs, which is rarely needed and can be unwieldy.
- **No CVR redaction on export** — Arlo exports the full CVR without redaction. Jurisdictions that are required to redact rare-style ballots must do so manually before publication.
- **No single "pre-seed inputs" bundle** that captures all inputs (contests, settings, manifests, CVRs) in one timestamped artifact.
- **No formal "checkpoint" or UI gate** that enforces exporting and acknowledging inputs before the seed is entered.
- **No Arlo-managed signing workflow** — sign-off tracking in `BatchInventoryData` stores a `signed_off_at` timestamp and `sign_off_user_id`, but this is an internal status flag, not a cryptographic signature or exported signed artifact.

### Guide: What to Do Today (Pre-Seed)

1. Once all jurisdictions have uploaded their manifests, CVRs (for ballot comparison), and batch tallies (for batch comparison), download them via the API — jurisdiction by jurisdiction, since observation happens primarily at the local (county) level. For batch comparison, use the manifests-bundle and candidate-totals-bundle endpoints to get the SHA-256-hashed ZIP bundles.
2. For ballot comparison audits, download each jurisdiction's CVR file. If your jurisdiction is required to redact CVRs before publication, apply the required redaction at this step. Compute a SHA-256 hash of each jurisdiction's CVR, manifest, and settings.
3. Export the contest settings and standardized contests file.
4. For each jurisdiction, assemble a per-jurisdiction ZIP of (manifest, CVR, relevant settings), compute a combined SHA-256 hash, and have the responsible official sign the hash (e.g., with PGP or a digital certificate). Publish these per-jurisdiction ZIPs to your public website or repository before entering the random seed.
5. Record the exact timestamp of the public publication.

---

## Stage 2: After the Seed, Before Comparisons — Publishing the Sample

### What Arlo Can Export Today

Once a round is started (seed entered, sample drawn), Arlo provides per-jurisdiction **ballot retrieval lists**:
- `GET /election/<id>/jurisdiction/<id>/round/<id>/ballots/retrieval-list` — CSV with container, tabulator, batch name, ballot position, imprinted ID, ticket numbers, already-audited flag, and audit board assignment.

For ballot comparison audits, the imprinted ID is included in the retrieval list. Since the CVR was published in Stage 1, an observer can look up the CVR interpretation for each sampled ballot using the imprinted ID as a key. Importantly, observers should do this join *after* the sample is public and *outside* of Arlo, so that audit boards never see pre-joined CVR interpretations that could script their judgments.

The **sample preview** endpoints allow the audit admin to estimate the sample size distribution across jurisdictions *before* officially launching a round. `POST /election/<id>/sample-preview` triggers a background computation given a set of proposed sample sizes; `GET /election/<id>/sample-preview` retrieves the result once complete. The response lists, for each jurisdiction, the projected number of sample draws (`numSamples`) and the number of unique ballots (`numUnique`) that would be selected. This is useful for workload planning but does not produce a binding sample — the actual sample is drawn when the round is officially started.

**The random seed** (`randomSeed`) is stored in `Election.random_seed` and is available via `GET /election/<id>/settings` once entered. It also appears in the final audit report header. The seed is typically a 20-digit number publicly drawn (e.g., from a live lottery). What is currently missing is a *structured bundle* that packages the seed together with the hashes of the pre-published inputs so observers can invoke `consistent_sampler` directly to verify the sample draw.

### Gaps

- **No combined cross-jurisdiction "all sampled ballots" export at round start time** — the retrieval list requires a separate download per jurisdiction.
- **No structured sampler-inputs bundle** — while the random seed is available from `GET /election/<id>/settings` and from the audit report, there is no single export that packages the seed together with the manifest hashes and CVR hashes used for sampling, so independent reproduction of the sample requires re-running Arlo or careful manual reconstruction.
- **Ticket numbers explain sampling but require the sampler algorithm to verify** — while `consistent_sampler` is open-source, no guided export walks observers through calling it with the exact inputs.

### Guide: What to Do Today (Post-Seed, Pre-Comparisons)

1. Immediately after starting the round, export the retrieval list for every jurisdiction (one download per jurisdiction, matching the local observation model).
2. Publish the random seed — it is available immediately via `GET /election/<id>/settings` (field `randomSeed`) and will also appear in the final audit report. The seed can be published as a plain number; it is not sensitive.
3. Publish each jurisdiction's retrieval list along with the seed and a timestamp. Observers can then independently join the retrieval list against the CVR file published in Stage 1 (using the imprinted ID as the key) to verify CVR interpretations — this join should happen externally, not inside Arlo.
4. Add a signed timestamp to the published retrieval lists so observers can verify that the sample was not changed after the seed was entered.

---

## Stage 3: After the Audit — Publishing Comparisons and Conclusions

### What Arlo Can Export Today

The primary post-audit export is the **audit admin audit report** (`GET /election/<id>/report`):
- **ELECTION INFO**: Organization, election name, state
- **CONTESTS**: Contest names, targeted/opportunistic, vote totals, number of winners
- **AUDIT SETTINGS**: Audit name, type, math type, risk limit, **random seed**, online/offline
- **AUDIT BOARDS**: Member names and affiliations (for online audits)
- **ROUNDS**: For each round and contest — sample size, risk limit met (Yes/No), **p-value**, start/end time, audited vote totals (CVR vs non-CVR for hybrid)
- **SAMPLED BALLOTS** (for ballot comparison):
  - Jurisdiction, container, tabulator, batch name, ballot position, imprinted ID
  - Ticket numbers per targeted contest
  - CVR result per contest
  - Audit result per contest
  - **Change in results** (vote delta)
  - **Change in margin** (discrepancy score, i.e., `counted_as` in supersimple terms: −2, −1, 0, +1, +2)
- **SAMPLED BATCHES** (for batch comparison):
  - Jurisdiction, batch, ballots in batch, ticket numbers, audited flag, reported results, audit results, change in results/margin

Additionally, a separate **discrepancy report** endpoint (`GET /election/<id>/discrepancy-report`) provides a focused view of the ballot or batch data without the summary header sections.

The **jurisdiction admin report** provides a per-jurisdiction view of the same sampled ballot/batch data.

### What's Included for Independent Reproducibility

For ballot comparison (supersimple), an observer can in principle independently reproduce the p-value if they have:
- The contest vote totals (in the report)
- The random seed (in the report)
- The ballot manifest and CVR files (published in Stage 1)
- The auditor's interpretations for each sampled ballot (in the report)

All of these are in the audit report or were published in earlier stages. The supersimple algorithm is open-source in Arlo's [`server/audit_math/supersimple.py`](../server/audit_math/supersimple.py).

For batch comparison (MACRO), the same is true with batch tallies instead of CVRs.

#### Sampling universe for each contest

An important detail for independent risk computation is knowing the *sampling universe* for each contest — the set of ballots (or batches) that were eligible to be sampled for it:

- For **targeted contests**, the universe is all ballots (or batches) in the contest's participating jurisdictions — those listed in the standardized contests file for that contest. The manifest for each jurisdiction defines how many ballots exist; the combined total is the denominator used in the risk calculation.
- For **opportunistic contests**, the universe is the set of ballots in the contest's participating jurisdictions that *happened to be sampled* because they were drawn for a targeted contest. Arlo internally computes this as the set of `SampledBallot`s in any jurisdiction participating in the opportunistic contest. The size of this universe determines whether the opportunistic contest achieves a meaningful risk level. This universe size is not currently exported as a standalone number.

For each contest, the participating jurisdictions are defined by the contest's entries in the standardized contests file, which is downloadable in Stage 1.

#### How `sampled_ballot_interpretations_to_cvrs` works

The internal function `sampled_ballot_interpretations_to_cvrs` (in `server/api/shared.py`) is the bridge between auditor input data and the risk computation. It converts each auditor's ballot interpretation into the `SAMPLECVRS` format expected by `supersimple.compute_risk`. For each sampled ballot it produces an entry with:
- `times_sampled`: for targeted contests, the number of ticket numbers the ballot received; for opportunistic contests, always 1
- `cvr`: a dict mapping each candidate/choice to 1 (voted) or 0 (not voted), or `None` if the ballot was `NOT_FOUND`, or `{}` if the contest was not on the ballot

The audit report's "CVR Result" and "Audit Result" columns encode this same information in human-readable form; the "Change in Margin" column encodes the `counted_as` discrepancy category (−2 to +2) that supersimple uses directly in its risk formula. An observer reproducing the p-value should map the "Change in Margin" column entries to the `counted_as` values in the supersimple algorithm.

### Gaps

- **Opportunistic contests: risk levels not computed or reported**: The audit report lists the CVR and audit results for opportunistic contests in sampled ballot rows, but Arlo does not compute or export a risk level (p-value) for them. To assess the risk achieved for an opportunistic contest, an observer must manually extract the relevant ballot rows, determine the sampling universe (see above), and run the risk computation independently.
- **Sampling universe not explicitly exported**: The number of ballots in each contest's sampling universe — needed to compute risk for both targeted and opportunistic contests — must currently be derived by the observer from the manifest files and contest-jurisdiction assignments, rather than being included directly in the audit report.
- **No p-value history for intermediate rounds**: The report shows the final p-value per round, but if the audit ran multiple rounds, the progression of risk reduction is not shown in a way that is easy to extract.
- **No machine-readable structured export**: The audit report is CSV, which is good, but it mixes section headers and data rows in a way that makes programmatic parsing non-trivial.
- **No self-contained reproducibility bundle**: There is no single export that bundles (a) the audit inputs, (b) the sample, (c) the interpretations, and (d) the risk computation result with enough metadata for a third party to run the risk computation independently without any Arlo access.
- **The "change in margin" column for ballot comparison shows the discrepancy category** (`counted_as`) not the raw vote delta: The values are −2, −1, 0, +1, +2 (supersimple categories), which is the right input for the risk formula, but an observer unfamiliar with the supersimple algorithm may not understand how to use these to reproduce the p-value.

### Guide: What to Do Today (Post-Audit)

1. Download the audit report CSV. This is the primary artifact.
2. Verify the p-values using the `supersimple.py` or `macro.py` algorithms. The report contains all the inputs: contest totals, sample size, discrepancy counts (change in margin values: −2 to +2).
3. For ballot comparison: cross-reference the "CVR Result" and "Audit Result" columns for each ballot with the original CVR files (published in Stage 1) to independently verify that the reported CVR interpretations match the CVR files.
4. For opportunistic contests: Arlo does not currently compute or export their risk levels. To compute them manually: identify all sampled ballots in jurisdictions participating in the opportunistic contest (from the audit report), determine the total ballot universe for that contest (from the manifests published in Stage 1 filtered to the contest's jurisdictions), look up CVR and audit interpretations, and run `supersimple.compute_risk` independently.
5. Publish per-jurisdiction ZIPs: for each jurisdiction, bundle the (manifest, CVR, retrieval list, audit interpretations from the report) and publish them. Together with the pre-seed Stage 1 bundle and the random seed, these form a complete, reproducible audit record.

---

## Recommendations for Improving Transparency Support

The following recommendations are ordered roughly by impact and implementation complexity.

### High Priority

1. **Add per-jurisdiction ZIP exports at each phase** — a convenience endpoint that assembles a per-jurisdiction ZIP of (manifest + CVR + relevant settings) at pre-seed, pre-round, and post-audit stages, with SHA-256 hashes for each file. This matches how observation actually works (at the county level) and gives officials a ready-to-sign artifact for each phase without requiring scripting.

2. **Add a ballot comparison CVR bundle endpoint with redaction support** (analogous to the existing `batch-files/manifests-bundle` and `batch-files/candidate-totals-bundle`). It should bundle all jurisdiction CVR files into a single ZIP with SHA-256 hashes, and optionally apply the redaction required by jurisdictions that must protect rare-style ballot privacy ([Bernhard et al., 2025](https://www.science.org/doi/10.1126/sciadv.adt1512)).

3. **Compute and export the sampling universe size for every contest** in the audit report — for both targeted and opportunistic contests. This is the denominator for the risk calculation and is the key piece of information needed to independently verify risk levels for opportunistic contests.

4. **Compute and report risk levels for opportunistic contests** in the audit report. Once the sampling universe is exported (recommendation 3), computing the p-value for opportunistic contests requires the same `supersimple.compute_risk` call already used for targeted contests, but restricted to ballots in the opportunistic contest's jurisdictions.

### Medium Priority

5. **Add a "pre-seed audit inputs" export** (`GET /election/<id>/pre-seed-bundle`) that produces a structured JSON or ZIP containing: contest settings, risk limit, jurisdictions list, and hashes of all uploaded manifest and CVR files. This would be the canonical artifact for officials to sign and publish before the seed is entered.

6. **Surface the sampler inputs as a published artifact**: After the sample is drawn, export a structured file containing the seed, the manifest hash, the CVR hash (for ballot comparison), and the `consistent_sampler` version used, so observers can call `consistent_sampler` directly to verify the sample draw.

7. **Add a structured (JSON) version of the audit report** alongside the current CSV, with clearly separated sections for easy programmatic parsing and independent risk computation, including explicit sampling universe sizes and opportunistic contest risk levels.

8. **Add a "reproducibility bundle"** post-audit that packages the audit report, the algorithm name and version, and a README explaining exactly how to reproduce the risk calculation (what inputs go into `supersimple.compute_risk` or `macro.compute_risk`, how to parse the relevant columns from the audit report, and which open-source code to use).

### Lower Priority

9. **Document the `CARD_STYLE_DATA` audit math type for public observers** — Arlo supports a card-style-data variant of ballot comparison that uses ballot style information from the CVR to compute contest-specific expected error rates (reducing sample sizes for contests that appear on only a subset of ballot styles). Currently this math type and its transparency implications are not documented for observers who need to replicate the risk computation.

10. **Integrate timestamp-authority support**: Allow officials to submit a hash of the pre-seed inputs bundle to a public timestamp service (e.g., RFC 3161) and record the receipt in Arlo's database. This provides verifiable proof of when the data was committed.

11. **Make the "change in margin" column in the audit report more transparent**: Add a companion column that explains the discrepancy in plain language (e.g., "2-vote overstatement") so that observers unfamiliar with the supersimple algorithm can understand the data.

12. **Formalize the sign-off workflow for audit admins**: Currently, the batch inventory sign-off tracks timestamps (`signed_off_at`, `sign_off_user_id`), but there is no audit-admin-level sign-off gate between the "inputs uploaded" and "seed entered" stages, nor between "sample drawn" and "comparisons entered." Adding explicit sign-off checkpoints with exported, signed receipts at each stage would support a complete chain of custody.

---

## Summary Table

| Transparency Need | Arlo Today | Gap / Note |
|---|---|---|
| Export all input data (manifests) | ✅ Per-jurisdiction CSV + bundle w/ SHA-256 (batch comparison) | ❌ No CVR bundle for ballot comparison |
| Export all input data (CVRs) | ✅ Per-jurisdiction CSV, all phases | ❌ No bundle/hash; no built-in redaction for rare-style privacy |
| Pre-seed signature artifact | ⚠️ All raw data is exportable | ❌ No structured per-jurisdiction ZIP or hash manifest; no sign-off gate |
| Post-seed: sample selection list | ✅ Per-jurisdiction retrieval list w/ imprinted ID | ❌ No combined list; observer CVR join must happen outside Arlo |
| Random seed availability | ✅ Available via settings endpoint and in audit report | ❌ No structured sampler-inputs bundle (seed + hashes) for independent verification |
| Post-audit: discrepancy/comparison export | ✅ Audit report + discrepancy report CSV | ❌ No structured JSON; sampling universe not explicit |
| Sampling universe per contest | ⚠️ Derivable from manifest + standardized contests file | ❌ Not exported as an explicit number in the audit report |
| Risk level for opportunistic contests | ❌ Not computed or exported | ❌ Significant gap; universe and p-value must be computed manually |
| CARD_STYLE_DATA (style-based selection) | ✅ Implemented in codebase | ❌ Not documented for observers; transparency implications unstated |
| Reproducibility documentation | ❌ Not provided as a bundle | ❌ Algorithm is open-source; no guided README or bundle |
| Public transparency endpoint | N/A — Arlo runs on a firewalled internal network | Officials export and publish data externally; no public endpoint needed or appropriate |
| Cryptographic signing support | ❌ No Arlo-native support | Officials sign exported artifacts externally; RFC 3161 timestamping would help |
| Activity/chain-of-custody log | ⚠️ Internal `ActivityLogRecord` table only | ❌ Not exportable |
