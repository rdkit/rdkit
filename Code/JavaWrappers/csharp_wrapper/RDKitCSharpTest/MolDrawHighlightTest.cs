using System.IO;
using GraphMolWrap;
using Xunit;

namespace RdkitTests;

/*
 *  Verifies that code to highlight multiple atoms/bonds is working in C#.
 *  Does not include any assertions
 */
public class MolDrawHighlightTest
{
    [Fact]
    public void TestDrawMolecule()
    {
        var testFile = "testMolHighLight.png";
        if (File.Exists(testFile))
        {
            File.Delete(testFile);
        }
        var mol = RWMol.MolFromSmiles("c1ccc(C)c(C)c1C");
        var drawer = new MolDraw2DCairo(300, 300, -1, -1, true);
        drawer.drawOptions().addAtomIndices = true;
        var red = new DrawColour(1.0, 0.0, 0.0);
        var green = new DrawColour(0.0, 1.0, 0.0);
        var blue = new DrawColour(0.0, 0.0, 1.0);
        var highlightAtoms = new Int_Vect {1, 2};
        var atomMap = new ColourPalette {new (1, red), new (2, green)};
        var highlightBonds = new Int_Vect {3, 4, 5};
        var bondMap = new ColourPalette {new (3, red), new (4, green), new (5, blue)};
        drawer.drawMolecule(mol, highlightAtoms, highlightBonds, atomMap, bondMap);
        drawer.drawMolecule(mol);
        drawer.finishDrawing();
        drawer.writeDrawingText(testFile);
    }
    
    [Fact]
    public void TestDrawMolecules()
    { 
        var testFile = "testMolsHighLight.png";
        if (File.Exists(testFile))
        {
            File.Delete(testFile);
        }
        var mol = RWMol.MolFromSmiles("c1ccc(C)c(C)c1C");
        var mol2 = RWMol.MolFromSmiles("c1ccc(C)c(C)c1");
        var mols = new ROMol_Ptr_Vect { mol, mol2 };
        var drawer = new MolDraw2DCairo(600, 300, 300, 300, true);
        drawer.drawOptions().addAtomIndices = true;
        var red = new DrawColour(1.0, 0.0, 0.0);
        var green = new DrawColour(0.0, 1.0, 0.0);
        var blue = new DrawColour(0.0, 0.0, 1.0);
        var legends = new Str_Vect { "Mol1", "Mol2" };
        var highlightAtoms = new Int_Vect_Vect {new Int_Vect {1, 2}, new Int_Vect {3, 4}};
        var atomMap = new ColourPalette_Vect
        {
            new ColourPalette{ new(1, red), new (2, green) },
            new ColourPalette{ new(2, red), new (3, green) }
        };
        var highlightBonds = new Int_Vect_Vect { new Int_Vect{3, 4, 5}, new Int_Vect {4, 5, 6}};
        var bondMap = new ColourPalette_Vect
        {
            new ColourPalette{ new(3, red), new (4, green), new (5, blue) },
            new ColourPalette{ new(4, red), new (5, green), new (6, blue) }
        };
        drawer.drawMolecules(mols, legends, highlightAtoms, highlightBonds, atomMap, bondMap);
        drawer.finishDrawing();
        drawer.writeDrawingText(testFile);
    }
}