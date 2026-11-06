-- load the library, so that all rdkit.* parameters are defined
SELECT rdkit_version() IS NOT NULL AS loaded;

-- rdkit.torsion_fp_size and rdkit.atompair_fp_size were documented by
-- mistake; they are accepted with a warning but have no effect
SET rdkit.torsion_fp_size = 512;
SET rdkit.atompair_fp_size = 512;
SELECT length(bfp_to_binary_text(torsionbv_fp('c1ccccc1CCO'::mol))) * 8 AS torsion_bits,
       length(bfp_to_binary_text(atompairbv_fp('c1ccccc1CCO'::mol))) * 8 AS atompair_bits;
SET rdkit.hashed_torsion_fp_size = 512;
SET rdkit.hashed_atompair_fp_size = 512;
SELECT length(bfp_to_binary_text(torsionbv_fp('c1ccccc1CCO'::mol))) * 8 AS torsion_bits,
       length(bfp_to_binary_text(atompairbv_fp('c1ccccc1CCO'::mol))) * 8 AS atompair_bits;
RESET ALL;

-- deprecated names are hidden from SHOW ALL / pg_settings
SELECT name FROM pg_settings
 WHERE name IN ('rdkit.hashed_torsion_fp_size', 'rdkit.hashed_atompair_fp_size',
                'rdkit.torsion_fp_size', 'rdkit.atompair_fp_size')
 ORDER BY name;
