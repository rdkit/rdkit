using System;
using GraphMolWrap;

public class rdktest
{
    // static void rxnTest() {
    //     Console.WriteLine( "Reaction tests" );
    //     var rxn = ChemicalReaction.ReactionFromSmarts("[N:1][C:2].[OH][C:3]=[O:4]>>[C:2][N:1][C:3]=[O:4]");
    //     var amine = RWMol.MolFromSmiles("CCCN"); 
    //     var acid = RWMol.MolFromSmiles("C1CC1CC(=O)O");
    //     ROMol[] rs = {amine,acid};                            
    //     ROMol_Vect rv = new ROMol_Vect(rs);
    //     for(var i=0;i<100000;i++){
    //         var ps=rxn.runReactants(rv);
    //         if(i%100 == 0) {
    //             Console.WriteLine( "\t{0}", i );
    //         }
    //     }
    //     Console.WriteLine( "Goodbye" );
    // }
    // static void smiTest() {
    //     Console.WriteLine( "repeatedly from smiles" );
    //     for(var i=0;i<1000000;i++){
    //         ROMol m1=RDKFuncs.MolFromSmiles("c1ccccc1");
    //         if(i%1000 == 0) {
    //             Console.WriteLine( "\t{0}", i );
    //         }
    //     }

    //     Console.WriteLine( "Goodbye" );
    // }
    
    static void morganTest() 
    {
        // ----- Object creation -----

        Console.WriteLine( "Creating some objects:" );

        ROMol m1=RWMol.MolFromSmiles("c1ccccc1");
        Console.WriteLine(" mol: "+m1+" "+m1.getNumAtoms());
        ROMol m2=RWMol.MolFromSmiles("c1ccccn1");

        var fp1=RDKFuncs.MorganFingerprintMol(m1,2);
        var fp2=RDKFuncs.MorganFingerprintMol(m2,2);

        Console.WriteLine(" sim: "+RDKFuncs.DiceSimilarity(fp1,fp2));
    }

    static void Main() 
    {
        // ----- Object creation -----

        Console.WriteLine( "Creating some objects:" );

        ROMol m1=RWMol.MolFromSmiles("c1ccccc1");
        Console.WriteLine(" mol: "+m1+" "+m1.getNumAtoms());
        ROMol m2=RWMol.MolFromSmiles("c1ccccn1");
        Console.WriteLine(" smi: "+m1+" "+m1.MolToSmiles());
        Console.WriteLine(" smi2: "+m2+" "+m2.MolToSmiles());            


        ExplicitBitVect fp1=RDKFuncs.LayeredFingerprintMol(m1);
        ExplicitBitVect fp2=RDKFuncs.LayeredFingerprintMol(m2);

        Console.WriteLine(" sim: "+RDKFuncs.TanimotoSimilarityEBV(fp1,fp2));

        //rxnTest();
        //smiTest();
        morganTest();

        Console.WriteLine( "Goodbye" );
    }
}
