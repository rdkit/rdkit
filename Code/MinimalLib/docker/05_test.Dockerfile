# Test NPM package of MinimalLib packs and publishes well

FROM rdkit-minimallib-rdkit-src:latest
COPY --from=rdkit-minimallib:latest / /src/rdkit/Code/MinimalLib/build/
WORKDIR /src/rdkit/Code/MinimalLib/build
RUN npm pack
RUN \
  VERSION=$(jq -r .version package.json) && \
  if [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then IS_NEXT=false; else IS_NEXT=true; fi && \
  if [ "$IS_NEXT" = true ]; then npm publish --dry-run --tag next; else npm publish --dry-run; fi
