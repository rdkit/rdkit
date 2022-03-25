using System.IO;
using System.Linq;
using GraphMolWrap;
using Xunit;

namespace RdkitTests;

public class MolDraw2DCairoTest
{
    [Fact]
    public void TestDrawMolecule()
    { 
        var mol = RWMol.MolFromSmiles("c1ccc(C)c(C)c1C");
        var drawer = new MolDraw2DCairo(300, 300, -1, -1, true);
        drawer.drawOptions().addAtomIndices = true;
        drawer.drawMolecule(mol);
        drawer.finishDrawing();
        var png1 = drawer.getImage().ToArray();
        drawer.writeDrawingText("test.png");
        byte[] png2 = File.ReadAllBytes("test.png");
        Assert.Equal(png1.Length, png2.Length);
        for (int i=0; i<png1.Length; i++) {
            Assert.Equal(png1[i], png2[i]);
        }
        File.Delete("test.png");
    }
}