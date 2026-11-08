# Test NPM package of MinimalLib packs and publishes well

FROM rdkit-minimallib-rdkit-src:latest
COPY --from=rdkit-minimallib:latest / /src/rdkit/Code/MinimalLib/build/
WORKDIR /src/rdkit/Code/MinimalLib/build
RUN /opt/emsdk/node/*/bin/npm pack
RUN /opt/emsdk/node/*/bin/npm publish --dry-run
