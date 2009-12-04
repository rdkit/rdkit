using System;

public class rdktest
{
    static void Main() 
    {
        // ----- Object creation -----

        Console.WriteLine( "Creating some objects:" );

        using (ROMol m = RDKFuncs.MolFromSmiles("c1ccccc1"))
        {
           Console.WriteLine(" mol: "+m+" "+m.getNumAtoms());
        }
        Console.WriteLine( "Goodbye" );
    }
}
