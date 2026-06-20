# Arlo Comparison Audit Transparency: Report and Guide

## Executive Summary

Arlo has meaningful but incomplete transparency infrastructure for comparison audits. It supports well-structured data uploads, CSV exports of inputs and outputs, and a comprehensive audit report. However, it lacks: (1) a formally-timed, bundled, pre-seed export of all audit inputs designed for public signing; (2) a mid-audit "ballot selection + CVR interpretation" export for independent observers before comparisons are recorded; and (3) cryptographic signing support built into the tool itself. Below is a stage-by-stage analysis, a guide for getting maximum transparency from what exists today, and specific recommendations for improvement.

---

## Stage 1: Before the Random Seed — Uploading Inputs

### What Arlo Can Export Today

**For ballot comparison audits**, the audit admin can download (via authenticated API endpoints):
- **Jurisdictions file** (`GET /election/<id>/jurisdictions/file`) — the list of participating counties/jurisdictions
- **Standardized contests file** (`GET /election/<id>/standardized-contests/file`) — the set of contests to be audited
- **Ballot manifest per jurisdiction** (`GET /election/<id>/jurisdiction/<id>/ballot-manifest/csv`) — batch/tabulator structure and ballot counts
- **CVR file per jurisdiction** (`GET /election/<id>/jurisdiction/<id>/cvrs/csv`) — the complete cast vote record for every ballot, with imprinted IDs

**For batch comparison audits**, additionally:
- **Batch tallies per jurisdiction** (`GET /election/<id>/jurisdiction/<id>/batch-tallies/csv`) — reported vote counts per batch
- **Batch inventory worksheet, CVR, and ballot manifest** (generated from batch inventory sign-off flow)
- **Bundle downloads**: Arlo recently added `POST /election/<id>/batch-files/manifests-bundle` and `POST /election/<id>/batch-files/candidate-totals-bundle` — these assemble all jurisdiction manifests or batch tallies into a single ZIP, compute a **SHA-256 hash** of the inner archive, and include a `*-sha256-hash.txt` file in the outer ZIP. This is the most significant transparency artifact for the pre-seed stage.

**Audit settings** (risk limit, math type, contest definitions with vote totals) are entered before the seed and appear in the final audit report.

### Gaps

- **No equivalent bundle for CVR files** in ballot comparison audits — the most critical input for ballot comparison audits lacks the batch-bundle treatment given to batch comparison manifests.
- **No single "pre-seed inputs" bundle** that captures all inputs (contests, settings, manifests, CVRs) in one timestamped artifact.
- **No formal "checkpoint" or UI gate** that enforces exporting and acknowledging inputs before the seed is entered.
- **All downloads require authenticated access** — there is no public transparency endpoint for pre-published input data.
- **No Arlo-managed signing workflow** — sign-off tracking in `BatchInventoryData` stores a `signed_off_at` timestamp and `sign_off_user_id`, but this is an internal status flag, not a cryptographic signature or exported signed artifact.

### Guide: What to Do Today (Pre-Seed)

1. Once all jurisdictions have uploaded their manifests, CVRs (for ballot comparison), and batch tallies (for batch comparison), download all of them via the API. For batch comparison, use the manifests-bundle and candidate-totals-bundle endpoints to get the SHA-256-hashed ZIP bundles.
2. For ballot comparison audits, use a script to iterate over each jurisdiction and download their CVR file, then manually compute a SHA-256 hash over the collection.
3. Export the contest settings and standardized contests file.
4. Combine all of the above into a single archive, compute a combined hash, and have election officials sign the hash (e.g., with PGP or a digital certificate). Publish this before entering the random seed.
5. Record the exact timestamp of the public publication and the Arlo audit settings (risk limit, math type, jurisdiction list, contest definitions) — these are all available before the seed is entered.

---

## Stage 2: After the Seed, Before Comparisons — Publishing the Sample

### What Arlo Can Export Today

Once a round is started (seed entered, sample drawn), Arlo provides per-jurisdiction **ballot retrieval lists**:
- `GET /election/<id>/jurisdiction/<id>/round/<id>/ballots/retrieval-list` — CSV with container, tabulator, batch name, ballot position, imprinted ID, ticket numbers, already-audited flag, and audit board assignment.

For ballot comparison audits, the imprinted ID is included in the retrieval list. Since the CVR was published in Stage 1, an observer can look up the CVR interpretation for each sampled ballot using the imprinted ID as a key.

The **sample preview** endpoint (`GET /election/<id>/sample-preview`) gives a preview of what would be sampled before officially starting the round.

### Gaps

- **No combined cross-jurisdiction "all sampled ballots" export at round start time** — the retrieval list requires a separate download per jurisdiction.
- **No export that pairs each sampled ballot with its CVR interpretation in a single file** — an observer must cross-reference the retrieval list against the full CVR files (which may each be hundreds of thousands of rows) to know what Arlo expects each ballot to show.
- **No public-facing endpoint** for the sample — publishing the sample requires the audit admin to download and publish the retrieval lists themselves.
- **Ticket numbers explain sampling but require the sampler algorithm to reproduce** — while `consistent_sampler` is open-source, no documentation or tool export bundles the exact parameters (seed + manifest hash + CVR hash) needed for an independent observer to reproduce the sample draw without re-running Arlo.

### Guide: What to Do Today (Post-Seed, Pre-Comparisons)

1. Immediately after starting the round, download the retrieval list for every jurisdiction.
2. For ballot comparison audits, programmatically join the retrieval list against the CVR CSV (using imprinted ID as the key) to produce a per-ballot file showing: jurisdiction, batch, ballot position, imprinted ID, ticket number(s), and the CVR interpretation for each audited contest. Publish this before any comparison data is entered.
3. Publish the random seed (it is visible in the audit settings section of the final report, and in the UI as soon as it is entered).
4. Publish the joined file with a timestamp and hash so observers can verify that the sample was determined entirely by the seed and pre-published inputs.

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

### Gaps

- **Opportunistic contests are not fully reported**: The audit report does include opportunistic contests in the sampled ballot rows, but the risk level for opportunistic contests is not computed or reported. A key transparency need — "what risk level did this audit achieve for each opportunistically-covered contest?" — is absent.
- **No p-value history for intermediate rounds**: The report shows the final p-value per round, but if the audit ran multiple rounds, the progression of risk reduction is not shown in a way that is easy to extract.
- **No machine-readable structured export**: The audit report is CSV, which is good, but it mixes section headers and data rows in a way that makes programmatic parsing non-trivial.
- **No self-contained reproducibility bundle**: There is no single export that bundles (a) the audit inputs, (b) the sample, (c) the interpretations, and (d) the risk computation result with enough metadata for a third party to run the risk computation independently without any Arlo access.
- **The "change in margin" column for ballot comparison shows the discrepancy category** (`counted_as`) not the raw vote delta: The values are −2, −1, 0, +1, +2 (supersimple categories), which is the right input for the risk formula, but an observer unfamiliar with the supersimple algorithm may not understand how to use these to reproduce the p-value.
- **No raw-data export of `sampled_ballot_interpretations_to_cvrs` output**: The internal function that converts all ballot interpretations to the CVR format used by `supersimple.compute_risk` is not exported directly; the audit report reconstructs enough of this but requires careful parsing to reconstruct the exact inputs to the risk calculation.

### Guide: What to Do Today (Post-Audit)

1. Download the audit report CSV. This is the primary artifact.
2. Verify the p-values using the `supersimple.py` or `macro.py` algorithms. The report contains all the inputs: contest totals, sample size, discrepancy counts (change in margin values).
3. For ballot comparison: cross-reference the "CVR Result" and "Audit Result" columns for each ballot with the original CVR files (published in Stage 1) to independently verify that the reported CVR interpretations match the CVR files.
4. For opportunistic contests: Arlo does not currently compute or export their risk levels. You would need to extract all ballots that were sampled for targeted contests but also contained opportunistic contests, look up their CVR and audit interpretations, and run the supersimple or MACRO risk calculation independently.
5. Publish the final audit report CSV along with the pre-seed bundle (Stage 1) and the ballot selection file (Stage 2). Together these form a complete, reproducible audit record.

---

## Recommendations for Improving Transparency Support

The following recommendations are ordered roughly by impact and implementation complexity.

### High Priority

1. **Add a ballot comparison CVR bundle endpoint** (analogous to the existing `batch-files/manifests-bundle` and `batch-files/candidate-totals-bundle`). It should bundle all jurisdiction CVR files into a single ZIP with a SHA-256 hash of the inner archive. This gives election officials a single artifact to sign before entering the seed.

2. **Add a "pre-seed audit inputs" export endpoint** (`GET /election/<id>/pre-seed-bundle`) that produces a single structured JSON or ZIP containing: contest settings, risk limit, jurisdictions list, and references (names + hashes) to all uploaded manifest and CVR files. This would be the canonical artifact for officials to sign and publish before the seed is entered.

3. **Add a "ballot sample + CVR expectations" export** for ballot comparison: a per-round CSV that joins each sampled ballot with its CVR interpretation for each audited contest. This should be available immediately after the round starts (before any comparisons are entered) so observers know exactly what the expected CVR interpretation is for each selected ballot.

4. **Compute and report risk levels for opportunistic contests** in the audit report. This is a significant transparency gap: observers cannot independently assess whether the audit provides meaningful risk reduction for contests that happened to appear on sampled ballots but were not targeted.

### Medium Priority

5. **Add a structured (JSON) version of the audit report** alongside the current CSV, with clearly separated sections for easy programmatic parsing and independent risk computation.

6. **Add a "reproducibility bundle"** post-audit that packages the audit report, the algorithm name and version, and a README explaining exactly how to reproduce the risk calculation (what inputs go into `supersimple.compute_risk` or `macro.compute_risk`, how to parse the relevant columns from the audit report, and which open-source code to use).

7. **Surface the sampler inputs as a published artifact**: After the sample is drawn, export the exact parameters used by `consistent_sampler` (seed, manifest interpretation, CVR interpretation) so an observer can call `consistent_sampler` directly to verify that the sample was drawn correctly from the published inputs.

8. **Add a public transparency API or public read-only audit page** that allows observers to access the audit report, sampled ballot list, and risk levels without authenticated access, gated only on the election official explicitly "publishing" the audit.

### Lower Priority

9. **Integrate timestamp-authority or cryptographic commitment support**: Allow officials to submit a hash of the pre-seed inputs bundle to a public timestamp service (e.g., RFC 3161) or a blockchain-based commitment scheme, and record the receipt URL/token in Arlo's database. This provides verifiable proof of when the data was committed.

10. **Make the "change in margin" column in the audit report more transparent**: Add a companion column that explains the discrepancy in plain language (e.g., "2-vote overstatement") so that observers unfamiliar with the supersimple algorithm can understand the data.

11. **Batch comparison: expose the full `sampled_batches_by_ticket_number` output** in the report** — the current ROUNDS section shows aggregate sampled batch counts and reported results, but the full ticket-number-to-batch mapping that is the direct input to `macro.compute_risk` is only available in the sampled batch section of the report, not in a form that directly mirrors the function's input format.

12. **Formalize the sign-off workflow for audit admins**: Currently, the batch inventory sign-off tracks timestamps, but there is no audit-admin-level sign-off gate between the "inputs uploaded" and "seed entered" stages, nor between "sample drawn" and "comparisons entered." Adding explicit sign-off checkpoints with exported, signed receipts at each stage would support a complete chain of custody.

---

## Summary Table

| Transparency Need | Arlo Today | Gap |
|---|---|---|
| Export all input data (manifests) | ✅ Per-jurisdiction + bundle w/ SHA-256 (batch comparison) | ❌ No CVR bundle for ballot comparison |
| Export all input data (CVRs) | ⚠️ Per-jurisdiction only, authenticated | ❌ No bundle/hash; no public endpoint |
| Pre-seed signature artifact | ❌ Not supported | ❌ No bundled, signable artifact |
| Post-seed: sample selection list | ✅ Per-jurisdiction retrieval list (w/ imprinted ID) | ❌ No combined list; no CVR interpretation joined |
| Post-audit: discrepancy/comparison export | ✅ Audit report + discrepancy report CSV | ❌ No structured JSON; not self-contained for reproduction |
| Risk level for opportunistic contests | ❌ Not computed or exported | ❌ Significant gap |
| Reproducibility documentation | ❌ Not provided | ❌ No bundle or README |
| Public transparency endpoint | ❌ All endpoints require auth | ❌ No public audit page |
| Cryptographic signing support | ❌ No Arlo-native support | ❌ Must be done externally |
| Activity/chain-of-custody log | ⚠️ Internal `ActivityLogRecord` table only | ❌ Not exportable |
