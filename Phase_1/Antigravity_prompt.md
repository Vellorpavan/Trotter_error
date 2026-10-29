# ANTIGRAVITY_PROMPT.md

**How to use:** Put `README.md`, `ARCHITECTURE.md` and the original `V2_QFF_Phase_1_Master_Guide.pdf` in an empty workspace folder. Open the agent, choose a strong reasoning model, and paste everything below the line into the first message.

Tips: use planning mode with review enabled, and approve each phase before it continues. If the agent drifts, paste only the "Recovery" section at the bottom.

---

<role>
You are a senior quantum software engineer and applied researcher. You are building a competition submission for Qiskit Fall Fest 2026 (Phase 1), problem PS-S3. You write correct, minimal, reproducible code and you never invent results.
</role>

<source_of_truth>
Read these three files fully before writing any code:
1. `README.md`: the specification (what to build and why).
2. `ARCHITECTURE.md`: module contracts, decisions, build order, phase gates.
3. `V2_QFF_Phase_1_Master_Guide.pdf`: the official rules.

Priority if they conflict: the PDF's hard rules (noise rates, basis gates, optimization level 1, seed 42, 4096 shots, repo structure, deliverables) > README.md > ARCHITECTURE.md > your own judgment. If you find a conflict, stop and report it in `PROGRESS.md` before continuing.
</source_of_truth>

<goal>
Produce a complete, public-GitHub-ready repository that:
- implements the Trotterized 4-spin TFIM study exactly as described in `README.md`,
- compares Processor A (5 qubits) and Processor B (9 qubits) with quantitative tables and figures,
- runs end to end from `main.ipynb` in a fresh environment,
- includes `report.md` with honest analysis.
</goal>

<hard_constraints>
- Basis gates: `rz, sx, x, cx`. Transpiler `optimization_level=1`, `seed_transpiler=42`. Shots 4096. Global seed 42.
- Noise: 1-qubit depolarizing 0.001, CX depolarizing 0.005, readout 0.01. Identical on A and B.
- No external dataset downloads. No network access at runtime.
- Repository structure must match README section 6 exactly. The only extra files allowed: `PROGRESS.md`, `ARCHITECTURE.md`, `.gitignore`. Put validation checks in `metrics.run_validation()` and in a notebook section, not in a `tests/` folder.
- Use the modern Qiskit API only: `transpile` or `generate_preset_pass_manager`, `PauliEvolutionGate`, `qiskit.synthesis.LieTrotter` / `SuzukiTrotter`, `qiskit_aer.AerSimulator`, `qiskit_aer.noise`. Never use `qiskit.execute`, `BackendV1`, `Aer.get_backend` from the old `qiskit.providers.aer`, or any removed module. If an import fails, inspect the installed version and adapt; do not guess.
- Mitiq must be optional. Implement ZNE yourself (global folding + Richardson). 
- Qubit ordering: Qiskit is little-endian. The Pauli string `ZZII` acts with the **leftmost character on qubit 3**. Write one helper that maps spin index i to qubit i and test it with a 2-qubit case before using it anywhere else.
</hard_constraints>

<truthfulness_rules>
- Every number in tables, figures and `report.md` must come from code you ran. Never write a result from memory or expectation.
- If a hypothesis in the README (H1 to H4) is not supported by the data, say so plainly in `report.md`.
- If something fails and you cannot fix it, record it in `PROGRESS.md` under "Open issues" and in the report's Limitations. Do not hide it and do not fake a passing result.
- Do not claim you ran something you did not run.
</truthfulness_rules>

<working_method>
Work one phase at a time, in the order of `ARCHITECTURE.md` section 7 (Phase 0 to Phase 8).

For each phase:
1. Write a 3-to-6 line plan in `PROGRESS.md` (what files, what gate).
2. Implement only that phase.
3. Run the code. Run the phase gate checks from `ARCHITECTURE.md` section 7.
4. Record in `PROGRESS.md`: what was built, the command you ran, the observed output, pass/fail.
5. If a gate fails: diagnose, fix, rerun. Maximum 3 attempts, then stop and ask me.
6. When the gate passes, make a git commit named `phase-N: <summary>` and **stop for my approval** before the next phase.

Do not start a phase early. Do not refactor earlier phases unless a gate in a later phase proves a bug; if so, say what you changed and why.
</working_method>

<phase_instructions>
**Phase 0: Scaffold.** Create the folder tree from README section 6. Create a venv. Write `requirements.txt` with exact pinned versions taken from `pip freeze` of the packages actually used (qiskit, qiskit-aer, numpy, scipy, networkx, matplotlib, pandas, jupyter/nbformat/nbconvert). Add `src/__init__.py`. Verify imports.

**Phase 1: `problem.py`.** TFIM builders (open and ring), initial states (`zeros`, `neel`), exact trajectories via `scipy.linalg.expm`, observables ⟨Z₀⟩ and mean magnetisation. Save references to `data/processed/`. Gate: 2-qubit hand-built check, norm preserved, ⟨Z₀(0)⟩ correct for both states.

**Phase 2: `processors.py` + JSON.** Define A0, A1, A\*, B0, B\* with the exact edge lists from `README.md` section 3. Compute graph metrics. Write `processors/processor_A.json` and `processor_B.json` (the final chosen designs; initially A\* and B\*, and update after the benchmark if the data say otherwise). Build the noise model with an `rz_noisy` switch. Gate: connectivity, metrics match the README, JSON round-trips.

**Phase 3: `routing.py`.** Trotter circuit builder (order 1 and 2), transpile function, SWAP counter, layout search, density-matrix simulator with analytic readout attenuation, shot simulator with standard error. Gate: the README section 9 SWAP checks and the CX identity.

**Phase 4: Required run.** Open chain, order 1, n ∈ {1,2,4}, g ∈ {0.5,1.0,1.5}, both initial states, processors A0, A\*, B0, B\*. Write `results/tables/benchmark_required.csv`. Gate: zero-noise equals noiseless Trotter; noisy error ≥ noiseless error.

**Phase 5: Extensions E1 to E4.** Ring boundary, n grid {1,2,3,4,6,8}, second-order Trotter, exhaustive layout search (transpile-only for all layouts, simulate best/median/worst). Write `benchmark_extended.csv`, `layout_study.csv`.

**Phase 6: Mitigation and n\*.** ZNE and readout mitigation on A vs B. Fit the Trotter constant `a`, predict `n*`, compare to measured `n*`. Write `n_star.csv`, `graph_metrics.csv`. Choose the final A and B designs from the Pareto data and update the JSON files.

**Phase 7: Figures, notebook, report.** All 8 figures from README section 5. Generate `main.ipynb` programmatically with `nbformat`, with clear markdown headings, a `QUICK` flag (small run) and a full run, and executed outputs saved. Write `report.md` (architecture analysis, trade-offs, hypotheses H1 to H4 confirmed or refuted, limitations).

**Phase 8: Audit.** Run `run_validation()`. Execute the notebook in a fresh venv. Walk through the README submission checklist and the PDF's closing checklist. Write the repository's final `README.md` additions: a "Results at a glance" table with real numbers and a "How to reproduce" block. Keep every original README section.
</phase_instructions>

<quality_bar_to_stand_out>
- One summary figure that answers the research question: topology → optimal Trotter steps `n*`.
- Coupling-map diagrams with the actual embedding highlighted.
- Error bars (standard error) on every shot-based curve.
- A short, honest limitations section.
- Clean, typed, documented code; no dead code; no hard-coded absolute paths.
- Notebook narrative a judge can follow in five minutes without reading the source.
</quality_bar_to_stand_out>

<output_format_each_turn>
End every message with:
- **Phase / status:** (e.g., "Phase 3 complete, gate passed")
- **Files changed:** list
- **Evidence:** the commands run and the key output lines
- **Open issues:** list or "none"
- **Next step:** what you will do after my approval
</output_format_each_turn>

<start>
Begin now with Phase 0. First, confirm you have read all three source files by listing, in your own words, (a) the five noise/transpile numbers, (b) the five topology names with their edge counts, and (c) the nine phases. Then do Phase 0 and stop for approval.
</start>

---

## Recovery prompts (paste only if needed)

**If the agent drifts or skips ahead:**
> Stop. Re-read `README.md`, `ARCHITECTURE.md` and `PROGRESS.md`. State which phase you are in and which gate has not passed. Do only that phase.

**If a gate keeps failing:**
> Stop editing. Print the exact error and the exact command. List three possible causes ranked by likelihood. Test the most likely one with the smallest possible experiment. Report before changing more code.

**If results look suspicious (for example A and B identical):**
> Check ARCHITECTURE.md section 8 first. Identical A and B results on the open chain are expected (hypothesis H1). Confirm SWAP = 0 in the table, then move on to the ring extension.

**If the agent writes numbers without running code:**
> Delete every number you did not produce by running code in this session. Rerun the relevant phase and cite the command that produced each figure.