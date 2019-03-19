package org.RDKit;

import org.junit.Test;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.stream.IntStream;

public class EditableSubstructLibraryTests extends GraphMolTest {

    File testData = new File(getRdBase(),
            "/Data/NCI/first_200.names.smi");

    private EditableSubstructLibraryTrustedSmilesWithPattern buildLibrary() {

        EditableSubstructLibraryTrustedSmilesWithPattern sssLib = new EditableSubstructLibraryTrustedSmilesWithPattern();
        int nAdded = 0;
        try (BufferedReader reader = new BufferedReader(new FileReader(testData))) {
            String smiles = null;
            while ((smiles = reader.readLine()) != null) {
                String[] terms = smiles.split("\\s+");
                assertEquals(terms.length, 2);
                String smi = terms[0];
                String id = terms[1];
                byte[] fp = sssLib.makeStringFingerprint(smi);
                sssLib.addSmiles(smiles, fp);
                nAdded++;
            }
        } catch (IOException ex) {
            fail("Unable to load test data");
        }
        //System.out.println("Added " + nAdded + " structures, library size " + sssLib.size());
        return sssLib;

    }

    private List<String> collectFromHitlist(StringChunkedHitlist hitlist) {
        List<String> allHits = new ArrayList<>();
        while (true) {
            Str_Vect hits = hitlist.next();
            if (hits.isEmpty()) break;
            //System.out.println("Got chunk size " + hits.size());
            for (int i = 0; i < hits.size(); i++) {
                allHits.add(hits.get(i));
            }
        }

        //System.out.println("Total number of hits " + allHits.size());
        return allHits;
    }

    @Test
    public void testBuildLibrary() throws IOException {
        EditableSubstructLibraryTrustedSmilesWithPattern sssLib = buildLibrary();

        ROMol query = RWMol.MolFromSmarts("[#6;$([#6]([#6])[!#6])]");

        StringChunkedHitlist hitlist = sssLib.getHitlistMatches(query, true, true, false, 10, -1);
        List<String> allHits = collectFromHitlist(hitlist);

        assertEquals(allHits.size(), 185);

        String[] javaIds = new String[]{"10", "20", "30", "40"};
        Str_Vect ids = new Str_Vect(javaIds.length);
        for (int i = 0; i < javaIds.length; i++) {
            ids.set(i, javaIds[i]);
        }

        long startSize = sssLib.size();
        String smi200 = sssLib.getMol("202").MolToSmiles(true);
        sssLib.removeMols(ids);
        long currentSize = sssLib.size();
        assertEquals(currentSize, startSize - 4);
        String currentSmi200 = sssLib.getMol("202").MolToSmiles(true);
        assertEquals(currentSmi200, smi200);
    }

    @Test
    public void testConcurrentSearch() {
        EditableSubstructLibraryTrustedSmilesWithPattern sssLib = buildLibrary();
        final CountDownLatch startGate = new CountDownLatch(1);
        final int nSearches = 5;
        final int nThreadsPerSearch = 2;
        final CountDownLatch endGate = new CountDownLatch(nSearches);
        final ROMol query = RWMol.MolFromSmarts("[#6;$([#6]([#6])[!#6])]");
        final Map<Integer, List<String>> allHits = new ConcurrentHashMap<>();
        IntStream.range(0, nSearches).forEach(searchNo ->  {
            Thread t = new Thread() {
                @Override
                public void run() {
                    try {
                        startGate.await();
                        //System.out.println("Running search "+searchNo);
                        StringChunkedHitlist hitlist = sssLib.getHitlistMatches(query, true, true, false, 10, nThreadsPerSearch);
                        List<String> hits = collectFromHitlist(hitlist);
                        allHits.put(searchNo, hits);
                    } catch (InterruptedException e) {
                        fail("Interrrupted exception");
                    } finally {
                        endGate.countDown();
                    }
                }
            };
            t.start();
        });
        startGate.countDown();
        try {
            endGate.await();
        } catch (InterruptedException e) {
            fail("Interrrupted exception");
        }

        Set<String> firstHits = new HashSet<>(allHits.get(0));
        assertEquals(allHits.size(), nSearches);
        allHits.values().forEach(hits -> {
            assertEquals(hits.size(), 185);
            assertEquals(firstHits, new HashSet<String>(hits));
        });

    }

    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.EditableSubstructLibraryTests");
    }

}
