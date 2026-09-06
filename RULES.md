# RULES.md

These rules are binding for the Ligora codebase and supersede prior informal guidance.

1. **Everything visible must work end-to-end.** If a feature exists in code, it has a real UI path; if a UI path exists, it runs against real data or real user input, not sample stubs.

2. **Do not reimplement shared open-source tools.** Use them, wrap them, or link out to them. Ligora is a UI and integration layer, not a new engine library.

3. **Data sources only, not hardcoded chemistry/biology intuition.** Classification, thresholds, weights, formulas, mapping rules, significance defaults, and similar numeric or categorical knowledge must come from datasets, API responses, or explicit user settings — not from buried constants in app code.

4. **Atomic and molecular constants must come from data sources, not app code.** Examples: atomic numbers, standard atomic weights, residue/ligand metadata, classification hints. Hardcoded chemical intuition is not allowed.

5. **Prefer authoritative live data.** Use RCSB, PubChem, ChEMBL, PDBBind, wwPDB CCD, UniProt, and other open endpoints where they exist. Cache responsibly, but do not invent values when a source can answer.

6. **Make rare and missing analyses usable.** If a dataset or analysis type already exists in the code but has no UI, finish the UI. If a useful analysis does not exist, add the smallest sensible layer that connects real data to the user.

7. **X86_64 desktop first, snap-ready.** Target desktop deployment and make packaging real, including data, desktop integration, and confinement semantics. No throwaway build scripts.

8. **No placeholder data in the product.** Sample data may exist only in tests, demos, or documentation. The app itself operates on user data or live sources.

9. **No fallbacks that invent values.** If a source or dataset is unavailable, the app must say so. It must not silently synthesize missing numbers or labels.

10. **Engines are integrated, not rebuilt.** Vina, gnina, PLIP, OpenMM, PySCF, and similar tools are integrated where appropriate. Their license and redistribution terms are respected and documented.

11. **Release packaging is first-class.** Snap and any other release artifacts are real, buildable, and verifiable. Packaging is not an afterthought.

12. **Industrial shipping posture.** The project targets an industrial-grade, market-shippable open-source product. Readability, correctness, maintainability, and release hygiene matter. Shortcuts that produce stubs or disposable implementation are not acceptable.
