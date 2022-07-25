**IMPORTANT NOTE**:

The NPM release process has now moved to the official repository of [rdkit-js](https://github.com/rdkit/rdkit-js).  

The official [rdkit-js](https://github.com/rdkit/rdkit-js) repository will now be the centralized place for everything built "on top of" the core RDKit MinimalLib source code ([this repository](https://github.com/rdkit/rdkit/tree/master/Code/MinimalLib)). Please read [this](https://github.com/rdkit/rdkit-js#introduction) for more context.

# RDKit MinimalLib <!-- omit in toc -->

[![Build Status](https://dev.azure.com/rdkit-js/rdkit-js/_apis/build/status/rdkit.rdkit-js?branchName=master)](https://dev.azure.com/rdkit-js/rdkit-js/_build/latest?definitionId=1&branchName=master)
[![License](https://img.shields.io/github/license/rdkit/rdkit)](https://github.com/rdkit/rdkit-js/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/10009991.svg)](https://zenodo.org/badge/latestdoi/10009991)  
[![NPM Latest Version](https://img.shields.io/npm/v/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Weekly Downloads](https://img.shields.io/npm/dw/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Monthly Downloads](https://img.shields.io/npm/dm/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Yearly Downloads](https://img.shields.io/npm/dy/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Total Downloads](https://img.shields.io/npm/dt/@rdkit/rdkit?label=total%20downloads)](https://www.npmjs.com/package/@rdkit/rdkit)

## Table of contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Install](#install)
- [Usage](#usage)
- [Live demos](#live-demos)
- [Building the MinimalLib](#building-the-minimallib)

## Introduction

The idea of the MinimalLib is to allow the [RDKit](https://github.com/rdkit/rdkit) to be used from JavaScript so that we can add chemical capabilities to web applications.  

This initial set of functionality does not cover all of RDKit's functionality, but it is intended to be directly useful.

## Install

The most popular way of installing the MinimalLib is with NPM.

```bash
npm i @rdkit/rdkit
# yarn add @rdkit/rdkit
```  

To build the MinimalLib manually, refer to [this section](#building-the-minimallib).

## Usage  

### Using the RDKit package assets

#### Option 1: Use the npm package distribution files

Once you have the RDKit package installed in your node modules, copy the following distribution files anywhere in your deployed assets.

- `node_modules/@rdkit/rdkit/dist/RDKit_minimal.js`
- `node_modules/@rdkit/rdkit/dist/RDKit_minimal.wasm`

**NOTE: Both files must be copied at the same location in your deployed assets for the library to work properly.**

#### Option 2: Use the remote distribution files from [unpkg.com](https://unpkg.com/)

- `https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js`
- `https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.wasm`

#### Option 3: Build your own distribution files

For this method, refer to [Building the MinimalLib](#building-the-minimallib).

### Running RDKit in your JavaScript code

To use RDKit, load the javascript file and instantiate the wasm module inside the `head` tag of your `index.html`, before you run your application code:

```html
<head>
    <!-- ...other files and HTML tags... -->
    <!-- Load the RDKit JS file -->
    <script src="https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js"></script>

    <!-- Instantiate the WASM module. The inline script below could live elsewhere inside your application code. -->
    <script>
        window.initRDKitModule()
            .then(function(RDKit) {
                console.log("RDKit version: " + RDKit.version());
                window.RDKit = RDKit;
                /**
                 * The RDKit module is now loaded.
                 * You can use it anywhere.
                 */
            })
            .catch(() => {
                // handle loading errors here...
            });
    </script>
    <!-- ...your application code goes here... -->
</head>

```

## Live demos

If you are using the MinimalLib for the first time, see the getting started examples at https://www.rdkitjs.com/ .

### All live demos

- RDKit.js website: https://www.rdkitjs.com/
- RDKit.js usage with React.js: https://react.rdkitjs.com/
- Legacy examples #1: https://unpkg.com/@rdkit/rdkit/dist/GettingStartedInJS.html
- Legacy examples #2: https://unpkg.com/@rdkit/rdkit/dist/demo.html

## Building the MinimalLib

Make sure you are at the root of the [MinimalLib](https://github.com/rdkit/rdkit/tree/master/Code/MinimalLib), and run the following script.

```bash
bash scripts/build_rdkitjs.sh <RDKit git release tag name>
# Example: bash scripts/build_rdkitjs.sh Release_2021_03_1
```

This command will take several minutes to complete, and will default to using the `master` branch if no version is provided. Also, checkout the `build_rdkitjs.sh` file and the minimallib `Dockerfile` to see how things are tied together.

Once you have verified that the distribution files have been properly added in `Code/MinimalLib/dist`, refer to the [Using the RDKit package assets](#using-the-rdkit-package-assets) section for the next steps.
