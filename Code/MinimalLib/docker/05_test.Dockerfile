# Test NPM package of MinimalLib packs and publishes well

FROM rdkit-minimallib-rdkit-src:latest
COPY --from=rdkit-minimallib:latest / /src/rdkit/Code/MinimalLib/build/
RUN apt-get install -y npm
WORKDIR /src/rdkit/Code/MinimalLib/build
RUN npm pack
RUN npm publish --dry-run
