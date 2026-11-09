# Sanitizer CI

The standalone **AddressSanitizer** workflow builds and tests one Linux x86-64
configuration with Clang 18 on Ubuntu 24.04. It runs the C++ and Boost.Python
CTest suites with address checking, including stack use-after-scope checks.
Reports fail the affected test and the job; CTest continues through the remaining
tests. Test logs and JUnit results are retained for 14 days, including on failure.

This initial configuration includes Cairo, FreeType, Eigen and XYZ2MOL. Optional
ChemDraw, CoordGen, maeparser, PubChem Shape, InChI, Avalon, FreeSASA and YAeHMOP
support is not covered. Leak checking, extra pointer-pair and initialization-order
checks, UBSan and TSan are deferred. There are no sanitizer suppressions or test
exclusions.

## Manual runs

Once `.github/workflows/asan.yml` exists on the repository's default branch,
users with write access can select **Actions → AddressSanitizer → Run workflow**
and choose a branch. The GitHub CLI equivalent is:

```sh
gh workflow run asan.yml --ref master
```

A draft PR alone does not make a new workflow manually runnable upstream. To
demonstrate it before merge, place the workflow on a fork's default branch and
run it there. This workflow does not run on pushes or pull requests.

## Optional weekly runs

Weekly execution is disabled by default. To enable it, create the repository
Actions variable `RUN_WEEKLY_SANITIZERS` with the value `true` in **Settings →
Secrets and variables → Actions → Variables**. Remove the variable or change its
value to disable weekly execution while retaining manual runs.

The schedule is Sunday at 04:23 UTC and tests the latest default-branch commit.
Until enabled, scheduled runs are skipped. Scheduled runs detect regressions
after merge and are not required PR checks. GitHub can delay or drop scheduled
runs, so the schedule is not a strict timing guarantee.

The initial resource budget is one runner, four build workers, two test workers,
a ten-minute timeout per test and a three-hour job timeout. Measure actual build
and test duration before expanding the compiler/platform matrix or sanitizer
coverage. Maintainers should agree who will triage failures before enabling the
schedule.

## Local validation

The workflow's configure, build and test commands can also be run locally from
a clean checkout with the listed dependencies installed. Set `CC`, `CXX`,
`PYTHON_EXECUTABLE` and `ASAN_SYMBOLIZER_PATH` to the matching local toolchain;
for the test command, also set `GITHUB_WORKSPACE`, `RDBASE`, `PYTHONPATH` and
`LD_LIBRARY_PATH` as shown in the workflow. Use a separate build tree when
switching compilers or sanitizer settings.

Local execution validates the instrumented build, runtime loading and tests.
GitHub trigger handling, hosted-runner dependency installation and artifact
uploads still need a GitHub Actions run.
