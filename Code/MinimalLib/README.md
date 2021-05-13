# RDKit for JavaScript (Official)

[![Azure build Status](https://dev.azure.com/rdkit-builds/RDKit/_apis/build/status/rdkit.rdkit?branchName=master)](https://dev.azure.com/rdkit-builds/RDKit/_build/latest?definitionId=1&branchName=master)
[![Documentation Status](https://readthedocs.org/projects/rdkit/badge/?version=latest)](https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/GettingStartedInJS.html)
[![License](https://img.shields.io/github/license/rdkit/rdkit)](https://github.com/rdkit/rdkit/blob/master/license.txt)
[![NPM Latest Version](https://img.shields.io/npm/v/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Weekly Downloads](https://img.shields.io/npm/dw/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Monthly Downloads](https://img.shields.io/npm/dm/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Yearly Downloads](https://img.shields.io/npm/dy/@rdkit/rdkit)](https://www.npmjs.com/package/@rdkit/rdkit)
[![NPM Total Downloads](https://img.shields.io/npm/dt/@rdkit/rdkit?label=total%20downloads)](https://www.npmjs.com/package/@rdkit/rdkit)

## Table of contents

  - [Introduction](#introduction)
  - [Install](#install)
  - [Using the RDKit package assets](#using-the-rdkit-package-assets)
    - [Option 1: Use the npm package distribution files](#option-1-use-the-npm-package-distribution-files)
    - [Option 2: Use the remote distribution files from unpkg.com](#option-2-use-the-remote-distribution-files-from-unpkgcom)
  - [Running RDKit in your JavaScript code](#running-rdkit-in-your-javascript-code)
  - [Usage](#usage)
  - [Live demos](#live-demos)
  - [Contributing](#contributing)
    - [Preparing a new release of the package](#preparing-a-new-release-of-the-package)
    - [Releasing a new beta version of the package](#releasing-a-new-beta-version-of-the-package)

## Introduction

**Note: This package should be considered experimental. The API is not yet stable and may change from release to release.**

The idea of this package is to allow the [RDKit](https://github.com/rdkit/rdkit) to be used from JavaScript so that we can add chemical capabilities to web applications.  

Rather than attempting a comprehensive wrapper (like the old [RDKitJS](https://github.com/rdkit/RDKitjs)), this exposes a small set of key functionality. I think the general approach, including this actual library, can be useful for other wrapper projects in the future.

This initial set of functionality is not complete, but it is intended to already be directly useful.

## Install

```bash
npm i @rdkit/rdkit
```  

## Using the RDKit package assets

### Option 1: Use the npm package distribution files

Once you have the RDKit package installed in your node modules, copy the following distribution files anywhere in your deployed assets.

- `node_modules/@rdkit/rdkit/CodeMinimalLib/dist/RDKit_minimal.js`
- `node_modules/@rdkit/rdkit/CodeMinimalLib/dist/RDKit_minimal.wasm`

**NOTE: Both files must be copied at the same location in your deployed assets for the library to work properly.**

### Option 2: Use the remote distribution files from [unpkg.com](https://unpkg.com/)

- `https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js`
- `https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.wasm`

## Running RDKit in your JavaScript code

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

## Usage

See the getting started demo at https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/GettingStartedInJS.html .

Follow the examples of this page to see the various ways to use the JavaScript release of RDKit.

## Live demos

- From this npm package: https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/GettingStartedInJS.html
- From this npm package: https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/demo.html

## Contributing

### Preparing a new release of the package

Make sure you are at the root of the [RDKit](https://github.com/rdkit/rdkit) GitHub project, and on the branch and version of the project you want to release. **Note that no commits should occur during the release process.**

#### Step 1: Set the release version in package.json

```bash
npm --no-git-tag-version version <semver version matching an RDKit release>
# Example npm --no-git-tag-version version 2021.3.1
```

#### Step 2: Build the distribution files

```bash
npm run build -- <RDKit git release tag name>
# Example: npm run build -- Release_2021_03_1
```

This command will default to using the `master` branch if no version is provided. Also, checkout the `build_rdkitjs.sh` file and the minimallib `Dockerfile` to see how things are tied together.

#### Step 3: Publish the package to npm

Once you have verified that the distribution files have been properly added in `Code/MinimalLib/dist`, publish the package:

```bash
npm publish --access public
```  

#### Step 4: Set back the placeholder version in package.json

```bash
npm run resetVersion
```

And you're done!

### Releasing a new beta version of the package

The process is the same as publishing a regular version, but the version specified and the npm publish command change slightly:

```bash
npm --no-git-tag-version version 2021.3.1-beta.0 # specify beta number in version here
npm publish --beta --access public # specify npm that it's a beta version
```
