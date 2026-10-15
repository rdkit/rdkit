# GitHub Actions CI

The Azure pipeline (`azure-pipelines.yml` and `.azure-pipelines/*.yml`) is converted to [.github/workflows/build_and_test.yml](.github/workflows/build_and_test.yml). The workflow has not been run or validated yet.

## Triggers

Pushes to `master`, `Release_*` and `dev/**`, and all pull requests. Azure's default PR trigger covers every PR, so this is kept.

## Jobs

All 13 active jobs are converted, each with its original runner, timeout and package versions. The commented-out `Ubuntu_18_04_x64_fuzzer` job is not converted.

## Differences from Azure

- **ccache:** the cache template became an `actions/cache` step in every build job. It falls back from the branch's cache to the `master` cache. The Azure "Print cache variables" step was not added.
- **Test results:** `PublishTestResults` became an upload of the CTest `Test.xml` files as an artifact. GitHub has no built-in CTest results viewer.
- **Docker job:** `Ubuntu_x64_minimallib_build` runs `docker compose build` with `GET_SRC=copy` and `COMPOSE_BAKE=false`.
- **Windows jobs:** these use `cmd` and add the runner's conda `Scripts` folder to the path.
- **Base job compiler:** `CC=gcc-11` and `CXX=g++-11` are set only on `Ubuntu_x64_base`, as Azure exposed those variables as environment variables there.
- **macOS:** the jobs still install the x86_64 Miniforge on `macos-14`, as Azure did.

## Bugs fixed from the Azure originals

- Unquoted `ipykernel>=6.0` was treated as a shell redirect; it is now quoted.
- Windows `libboost-python=(boost_version)` was malformed; it is now `%BOOST_VERSION%`.
- Windows `-DCMAKE_LIBRARY_PATH` had an unbalanced quote, which was removed.
