# Zenodo DOI (GitHub release → archived record)

This project can be archived in Zenodo via the GitHub–Zenodo integration so that an exact, versioned release receives a citable DOI.

## Current DOI
- Version DOI (v1.0.4): https://doi.org/10.5281/zenodo.18733963
- Concept DOI (all versions): https://doi.org/10.5281/zenodo.18733930

## Recommended path (GitHub release → Zenodo DOI)
1. Sign in to Zenodo and connect your GitHub account.
2. In Zenodo, enable archiving for the repository:
   `Zhenghongwei11/Benchmarking-spatial-domain-identification-in-colorectal-cancer-Visium-data`.
3. Create a GitHub release (e.g., `v1.0.0`) from the repository.
4. Zenodo will archive the release and mint a DOI for that version.
5. Update the public repository documentation with the DOI link (e.g., `https://doi.org/10.5281/zenodo.XXXXXXX`), and update the journal submission materials maintained outside this repository.

## Notes
- Keep large binary artifacts (e.g., the canonical reproducibility ZIP) as release assets rather than committing them to git history.
- The manuscript DOCX can be included in the release assets if desired, but the Markdown source is the primary, diffable format.
