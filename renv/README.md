# renv Setup

This project does not include a pre-generated `renv.lock` because `renv` was not available in the local workspace at assembly time.

To initialize the environment:

1. Run `Rscript scripts/00_bootstrap_renv.R`
2. Allow the script to install `renv` and the declared analysis packages
3. Commit the generated `renv.lock` once package installation and snapshotting complete

The `.Rprofile` file in the project root will auto-load `renv/activate.R` after bootstrap has been completed.
