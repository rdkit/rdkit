package org.RDKit;

import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import static org.junit.Assert.assertEquals;

public class MolDraw2DCairoTests extends GraphMolTest{

    @Test
    public void createImage() throws IOException {
        RWMol mol = RWMol.MolFromSmiles("c1ccc(C)c(C)c1C");
        MolDraw2DCairo drawer = new MolDraw2DCairo(300, 300, -1, -1, true);
        drawer.drawOptions().setAddAtomIndices(true);
        drawer.drawMolecule(mol);
        drawer.finishDrawing();
        byte[] png1 = drawer.toByteArray();
        // Files.write(Paths.get("test2.png"), png);
        drawer.writeDrawingText("test.png");
        byte[] png2 = Files.readAllBytes(Paths.get("test.png"));
        assertEquals(png1.length, png2.length);
        for (int i=0; i<png1.length; i++) {
            assertEquals(png1[i], png2[i]);
        }
        Files.delete(Paths.get("test.png"));
    }


    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("org.RDKit.MolDraw2DCairoTests");
    }
}
