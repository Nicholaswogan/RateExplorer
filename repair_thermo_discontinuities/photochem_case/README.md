# TOI-1231 b Photochem test case

This folder holds the diagnostic code, its Python dependencies copied from the
original Photochem case, and the climate and spectrum input files. From the
parent repair directory, run:

```sh
conda run -n repair python photochem_case/run_case.py
```

The runner regenerates the repaired data and case mechanisms, then runs the
diagnostic here. The generated YAMLs and JSONL run log stay in this folder and
are ignored by its `.gitignore`. The two final repaired YAMLs are in the parent
directory; the pinned source checkout and validation metrics stay under
`../results/`.
