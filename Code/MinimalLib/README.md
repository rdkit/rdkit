# RDKit MinimalLib (RDKit.js)

[![NPM Latest Version](https://img.shields.io/npm/v/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Weekly Downloads](https://img.shields.io/npm/dw/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Monthly Downloads](https://img.shields.io/npm/dm/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Yearly Downloads](https://img.shields.io/npm/dy/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Total Downloads](https://img.shields.io/npm/dt/@rdkit/rdkit?label=total%20downloads)](https://www.npmjs.com/package/@rdkit/rdkit)

RDKit.js is the official JavaScript distribution of cheminformatics functionality from [RDKit](https://github.com/rdkit/rdkit) -
a collection of cheminformatics and machine-learning software written in C++.
A new NPM package is published for every RDKit release.

The core WASM module comes from RDKit's [MinimalLib](https://github.com/rdkit/rdkit/tree/master/Code/MinimalLib).
MinimalLib is a C++ layer that wraps a subset of RDKit's API so it can be compiled to WebAssembly and used from JavaScript.
The package is built and published directly from RDKit, while keeping JavaScript documentation at [rdkitjs.com](https://rdkitjs.com).

The package itself consists of three files:

* `RDKit_minimal.js` - Standard JavaScript wrapper for loading WASM modules
* `RDKit_minimal.wasm` - The compiled RDKit MinimalLib WASM binary
* `RDKit_minimal.d.ts` - TypeScript interface generated during compilation.

That means the package has zero dependencies, and if high-level component JavaScript is needed, it needs to be implemented yourself and won't be included in the general package.
This is to ensure easy maintenance of the package.

## Install RDKit JS

You can install it using one of the many (and growing) JavaScript package managers.

```bash
npm i @rdkit/rdkit
pnpm i @rdkit/rdkit
yarn add @rdkit/rdkit
bun add @rdkit/rdkit
...
```

Or use a CDN by adding this script tag to your HTML.

```html
<script src="https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js"></script>
```

## Getting Started and Loading the WASM module

See [rdkitjs.com](https://rdkitjs.com) for example code, documentation, and demos.

Usually, a trick is needed to make the `.wasm` file available as a standalone file, but it varies from framework to framework.
We have created examples for Vanilla JS, React, Vue, Angular, Svelte, Next.js, and Node.js.
However, this is not RDKit-specific; it is generally how you support WASM for those frameworks.

## Building the MinimalLib

Go to the [MinimalLib](https://github.com/rdkit/rdkit/tree/master/Code/MinimalLib) directory, and run the following script:

```bash
scripts/build_rdkitjs.sh
```

The script results in a `build` directory with a valid JavaScript module, and includes `RDKit_minimal.{js,wasm,d.ts}`.
To set the version of the JavaScript module use `NPM_VERSION`, e.g., `export NPM_VERSION=1.2.3`.
Defaults to `0.0.0`.

Set `GET_SRC=clone` to swap the local copy for a fresh clone from GitHub.
This is relevant when you want to pass a release tag (defaults to `master`) and a git URL (defaults to `https://github.com/rdkit/rdkit.git`):

```bash
GET_SRC=clone scripts/build_rdkitjs.sh <RDKit git release tag name>
# Example: GET_SRC=clone scripts/build_rdkitjs.sh Release_2025_03_2
```

The full build will take 20-30 minutes to complete.

## License

The binary is compiled directly from RDKit, so the license is unchanged.
BSD 3-Clause.

## Citation

See [rdkit.com](https://rdkit.com) for citation.
Note the version installed.
