# JavaScript wrappers for the RDKit

**Note** These wrappers should be considered experimental. The API is not yet stable and may change from release to release.

The idea here is to allow the RDKit to be used from JavaScript so that we can add chemical capabilities to web applications.
Rather than attempting a comprehensive wrapper (like RDKitJS), this exposes a small set of key functionality. I think the general approach, including this actual library, can be useful for other wrapper projects in the future.

This initial set of functionality is not complete, but it is intended to already be directly useful.

The `Dockerfile` in the `docker/` shows how to setup an appropriate environment and build the wrappers.

