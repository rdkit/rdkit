# include minimallib readme during npm publish process
# not main RDKit readme
# NOTE: this script goes along with the npm_postpublish.sh script
mv README.md MAIN_RDKIT_README.md
mv README MAIN_RDKIT_README
cp Code/MinimalLib/README.md README.md
