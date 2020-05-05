## Important Notice
The files in the corpora folders (i.e. the folders ending in `_fuzzer`) can not be directly used for purposes other than fuzzing.
This is because the fuzzer uses parts of the content for generating different information.
Consider the example `[OH3+]0`.
The first part `[OH3+]` will be used as a smiles formula, but the last part `0` will for example be used to determine
whether the fuzzer should set a certain flag to `true` or it will be used to derive an integral value.