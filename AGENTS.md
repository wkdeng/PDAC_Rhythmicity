# PDAC_Rhythmicity publication repository

- The six `notebook_publication/` notebooks are the authoritative figure entry points. Follow the scoped publication guidance there.
- Preserve this publication-only scope. Do not add legacy notebooks, data, generated figures, caches or notebook execution outputs to Git.
- Retain self-contained notebook runtime code, documented inputs, explicit scientific thresholds, provenance and QC. Do not substitute synthetic coverage or silently change analysis methods.
- Preserve the existing `pipeline/` layout. Default workflow target is TACC/Slurm with explicit per-rule Apptainer and strict shell behavior; retain the existing BSIC alternative.
- Before workflow submission or changes, run `./run.sh -n` or `./run_rmats.sh -n` from `pipeline/` in the configured environment.
- Put figure source tables in `data/publication/<analysis_slug>/`, figures in `data/figures/publication/<analysis_slug>/`, and logs in `data/logs/publication/<analysis_slug>/`.
- Keep `data/README.md` current when input requirements change. Preserve biological parameters when adapting paths or environments.
