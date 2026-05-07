#include <iostream>
#include "../sketcherMinimizer.h"

int main()
{
    sketcherMinimizer minimizer;

    /* create a molecule */
    auto* min_mol = new sketcherMinimizerMolecule();

    /* add an atom and set its parameters */
    auto a1 = min_mol->addNewAtom();
    a1->setAtomicNumber(7);

    auto a2 = min_mol->addNewAtom();
    a2->setAtomicNumber(6);

    /* add a bond and set its parameters */
    auto b1 = min_mol->addNewBond(a1, a2);
    b1->setBondOrder(1);

    /* load minimizer */
    minimizer.initialize(min_mol);

    /* generate coordinates */
    minimizer.runGenerateCoordinates();

    /* print coordinates */
    auto c1 = a1->getCoordinates();
    auto c2 = a2->getCoordinates();
    std::cerr << c1 << "  " << c2 << std::endl;
    return 0;
}
