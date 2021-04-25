# RDKit

## **WARNING: THIS IS NOT AN OFFICIAL RELEASE OF RDKIT - THIS IS ONLY A PROPOSAL FROM THE PUBLISHER (GitHub user @michelml). USE THIS AT YOUR OWN RISK.**

## Table of contents

  - [Introduction](#introduction)
  - [Install](#install)
  - [Using the rdkit package assets](#using-the-rdkit-package-assets)
    - [Option 1: Inside your JavaScript build](#option-1-inside-your-javascript-build)
    - [Option 2: Using remote distribution files from unpkg.com](#option-2-using-remote-distribution-files-from-unpkgcom)
  - [Running rdkit in your JavaScript code](#running-rdkit-in-your-javascript-code)
  - [Usage](#usage)
  - [Live demos](#live-demos)
  - [Contributing](#contributing)
    - [Building the distribution files of a new release](#building-the-distribution-files-of-a-new-release)
    - [Releasing a new npm version of the package](#releasing-a-new-npm-version-of-the-package)
    - [Releasing a new beta version of the package](#releasing-a-new-beta-version-of-the-package)

## Introduction

**Note:  This package should be considered experimental. The API is not yet stable and may change from release to release.**

The idea of this package is to allow the [RDKit](https://github.com/rdkit/rdkit) to be used from JavaScript so that we can add chemical capabilities to web applications.
Rather than attempting a comprehensive wrapper (like the old [RDKitJS](https://github.com/rdkit/RDKitjs)), this exposes a small set of key functionality. I think the general approach, including this actual library, can be useful for other wrapper projects in the future.

This initial set of functionality is not complete, but it is intended to already be directly useful.

The `Dockerfile` in the `docker/` shows how to setup an appropriate environment and build the wrappers.

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

### Option 2: Use the remote distribution files from unpkg.com

- `https://unpkg.com/@rdkit/rdkit@2021.3.1-beta.0/Code/MinimalLib/dist/RDKit_minimal.js`
- `https://unpkg.com/@rdkit/rdkit@2021.3.1-beta.0/Code/MinimalLib/dist/RDKit_minimal.wasm`


## Running RDKit in your JavaScript code

To use RDKit, load the javascript file and instantiate the wasm module inside the `head` tag of your `index.html`, before you run your application code:

```html
<head>
    <!-- ...other files and html tags... -->
    <!-- Load the JS file -->
    <script src="https://unpkg.com/@rdkit/rdkit@2021.3.1-beta.0/Code/MinimalLib/dist/RDKit_minimal.js"></script>
    <!-- Instantiate the WASM module -->
    <script>
        window.initRDKitModule()
            .then(function(RDKit) {
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

The best examples on how you can use RDKit are currently at https://rdkit.org/rdkitjs/beta/GettingStartedInJS.html .

Inspect the source code of this page to see the various ways to use the JavaSciprt release of RDKit.

## Live demos

- Getting Started: https://rdkit.org/rdkitjs/beta/GettingStartedInJS.html
- From this npm package: https://unpkg.com/@rdkit/rdkit@2021.3.1-beta.0/Code/MinimalLib/dist/demo.html

## Contributing

### Building the distribution files of a new release

Make sure you are at the root of the [RDKit](https://github.com/rdkit/rdkit) GitHub project, and on the branch and version of the project you want to release. Then, run the following command:

```bash
npm run build -- <RDKit release version name>
```  

This command will default to using the `master` branch if no version is provided.

### Releasing a new npm version of the package

**Note: To release new versions, you need to be a member of the @rdkit npm organization.**

Make sure you are login to npm with a user member of the @rdkit npm organization:

```bash
npm login
```

Then, simply version and publish the package on npm:

```bash
npm version 2021.3.1
git push origin master # npm version creates a git commit
npm publish --access public
```

### Releasing a new beta version of the package

The process is the same as publishing a regular version, but the version specified and the npm publish command changes slightly:

```bash
npm version 2021.3.1-beta.0 # specify beta number in version here
git push origin master
npm publish --beta --access public # specify npm that it's a beta version
```
