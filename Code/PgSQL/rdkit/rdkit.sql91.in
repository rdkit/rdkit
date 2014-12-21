CREATE OR REPLACE FUNCTION rdkit_version()
RETURNS text 
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_in(cstring)
RETURNS mol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_out(mol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_recv(internal)
RETURNS mol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_send(mol)
RETURNS bytea
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION qmol_recv(internal)
RETURNS qmol
AS 'MODULE_PATHNAME', 'mol_recv'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION qmol_send(qmol)
RETURNS bytea
AS 'MODULE_PATHNAME', 'mol_send'
LANGUAGE C STRICT IMMUTABLE;

CREATE TYPE mol (
        INTERNALLENGTH = -1,
        INPUT = mol_in,
        OUTPUT = mol_out,
        RECEIVE = mol_recv,
        SEND = mol_send,
        STORAGE = extended
);


CREATE OR REPLACE FUNCTION qmol_in(cstring)
RETURNS qmol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION qmol_out(qmol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE TYPE qmol (
        INTERNALLENGTH = -1,
        INPUT = qmol_in,
        OUTPUT = qmol_out,
        RECEIVE = qmol_recv,
        SEND = qmol_send,
        STORAGE = extended
);

CREATE OR REPLACE FUNCTION is_valid_smiles(cstring)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_from_smiles(cstring)
RETURNS mol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_from_smarts(cstring)
RETURNS mol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_to_smiles(mol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_to_smarts(mol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_to_smiles(qmol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_to_smarts(qmol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION is_valid_smarts(cstring)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION is_valid_ctab(cstring)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_from_ctab(cstring,bool default false)
RETURNS mol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_to_ctab(mol,bool default true)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION is_valid_mol_pkl(bytea)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_from_pkl(bytea)
RETURNS mol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_to_pkl(mol)
RETURNS bytea
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_in(cstring)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_out(bfp)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE TYPE bfp (
        INTERNALLENGTH = -1,
		INPUT = bfp_in,
		OUTPUT = bfp_out,
		STORAGE = extended
);

CREATE OR REPLACE FUNCTION sfp_in(cstring)
RETURNS sfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION sfp_out(sfp)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE TYPE sfp (
        INTERNALLENGTH = -1,
		INPUT = sfp_in,
		OUTPUT = sfp_out,
		STORAGE = extended
);

CREATE OR REPLACE FUNCTION bfp_from_binary_text(bytea)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_to_binary_text(bfp)
RETURNS bytea
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION layered_fp(mol)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION rdkit_fp(mol)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION morganbv_fp(mol,int default 2)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION morgan_fp(mol,int default 2)
RETURNS sfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION featmorganbv_fp(mol,int default 2)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION featmorgan_fp(mol,int default 2)
RETURNS sfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION atompair_fp(mol)
RETURNS sfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION torsion_fp(mol)
RETURNS sfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION atompairbv_fp(mol)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION torsionbv_fp(mol)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION maccs_fp(mol)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION avalon_fp(mol,bool default false,int default 15761407)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;




CREATE OR REPLACE FUNCTION tanimoto_sml(bfp, bfp)
RETURNS float8 
AS 'MODULE_PATHNAME', 'bfp_tanimoto_sml'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION dice_sml(bfp, bfp)
RETURNS float8 
AS 'MODULE_PATHNAME', 'bfp_dice_sml'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION tversky_sml(bfp, bfp, float4, float4)
RETURNS float8 
AS 'MODULE_PATHNAME', 'bfp_tversky_sml'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION tanimoto_dist(bfp, bfp)
RETURNS float8
AS 'MODULE_PATHNAME', 'bfp_tanimoto_dist'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION dice_dist(bfp, bfp)
RETURNS float8
AS 'MODULE_PATHNAME', 'bfp_dice_dist'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION tanimoto_sml_op(bfp, bfp)
RETURNS bool 
AS 'MODULE_PATHNAME', 'bfp_tanimoto_sml_op'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OPERATOR % (
	LEFTARG = bfp,
	RIGHTARG = bfp,
	PROCEDURE = tanimoto_sml_op(bfp, bfp),
	COMMUTATOR = '%',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION dice_sml_op(bfp, bfp)
RETURNS bool 
AS 'MODULE_PATHNAME', 'bfp_dice_sml_op'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OPERATOR # (
	LEFTARG = bfp,
	RIGHTARG = bfp,
	PROCEDURE = dice_sml_op(bfp, bfp),
	COMMUTATOR = '#',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION size(bfp)
RETURNS int4
AS 'MODULE_PATHNAME', 'bfp_size'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION tanimoto_sml(sfp, sfp)
RETURNS float8 
AS 'MODULE_PATHNAME', 'sfp_tanimoto_sml'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION dice_sml(sfp, sfp)
RETURNS float8 
AS 'MODULE_PATHNAME', 'sfp_dice_sml'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION tanimoto_sml_op(sfp, sfp)
RETURNS bool 
AS 'MODULE_PATHNAME', 'sfp_tanimoto_sml_op'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OPERATOR % (
	LEFTARG = sfp,
	RIGHTARG = sfp,
	PROCEDURE = tanimoto_sml_op(sfp, sfp),
	COMMUTATOR = '%',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION dice_sml_op(sfp, sfp)
RETURNS bool 
AS 'MODULE_PATHNAME', 'sfp_dice_sml_op'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OPERATOR # (
	LEFTARG = sfp,
	RIGHTARG = sfp,
	PROCEDURE = dice_sml_op(sfp, sfp),
	COMMUTATOR = '#',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION add(sfp, sfp)
RETURNS sfp
AS 'MODULE_PATHNAME', 'sfp_add'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION subtract(sfp, sfp)
RETURNS sfp
AS 'MODULE_PATHNAME', 'sfp_subtract'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION all_values_gt(sfp, int) 
RETURNS bool
AS 'MODULE_PATHNAME', 'sfp_allvals_gt'
LANGUAGE C STRICT IMMUTABLE COST 10;

CREATE OR REPLACE FUNCTION all_values_lt(sfp, int)
RETURNS bool
AS 'MODULE_PATHNAME', 'sfp_allvals_lt'
LANGUAGE C STRICT IMMUTABLE COST 10;




CREATE OR REPLACE FUNCTION mol_amw(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_logp(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_hba(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_hbd(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numrotatablebonds(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numatoms(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numheavyatoms(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numheteroatoms(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numrings(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numaromaticrings(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numaliphaticrings(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numsaturatedrings(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numaromaticheterocycles(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;
CREATE OR REPLACE FUNCTION mol_numaliphaticheterocycles(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;
CREATE OR REPLACE FUNCTION mol_numsaturatedheterocycles(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_numaromaticcarbocycles(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;
CREATE OR REPLACE FUNCTION mol_numaliphaticcarbocycles(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;
CREATE OR REPLACE FUNCTION mol_numsaturatedcarbocycles(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;
CREATE OR REPLACE FUNCTION mol_numheterocycles(mol)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;



CREATE OR REPLACE FUNCTION mol_fractioncsp3(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_tpsa(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi0n(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi1n(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi2n(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi3n(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi4n(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi0v(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi1v(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi2v(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi3v(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_chi4v(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_kappa1(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_kappa2(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_kappa3(mol)
RETURNS real
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;


CREATE OR REPLACE FUNCTION mol_formula(mol,bool default false, bool default true)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;



CREATE OR REPLACE FUNCTION mol_inchi(mol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_inchikey(mol)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_murckoscaffold(mol)
RETURNS mol
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;



CREATE OR REPLACE FUNCTION substruct(mol, mol)
RETURNS bool
AS 'MODULE_PATHNAME', 'mol_substruct'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION substruct_count(mol, mol,bool default true)
RETURNS int
AS 'MODULE_PATHNAME', 'mol_substruct_count'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR @> (
	LEFTARG = mol,
	RIGHTARG = mol,
	PROCEDURE = substruct(mol, mol),
	COMMUTATOR = '<@',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION substruct(mol, qmol)
RETURNS bool
AS 'MODULE_PATHNAME', 'mol_substruct'
LANGUAGE C STRICT IMMUTABLE;
CREATE OPERATOR @> (
	LEFTARG = mol,
	RIGHTARG = qmol,
	PROCEDURE = substruct(mol, qmol),
	COMMUTATOR = '<@',
	RESTRICT = contsel,
	JOIN = contjoinsel
);


CREATE OR REPLACE FUNCTION rsubstruct(mol, mol)
RETURNS bool
AS 'MODULE_PATHNAME', 'mol_rsubstruct'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION rsubstruct(qmol, mol)
RETURNS bool
AS 'MODULE_PATHNAME', 'mol_rsubstruct'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR <@ (
	LEFTARG = mol,
	RIGHTARG = mol,
	PROCEDURE = rsubstruct(mol, mol),
	COMMUTATOR = '@>',
	RESTRICT = contsel,
	JOIN = contjoinsel
);
CREATE OPERATOR <@ (
	LEFTARG = qmol,
	RIGHTARG = mol,
	PROCEDURE = rsubstruct(qmol, mol),
	COMMUTATOR = '@>',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION mol_cmp(mol,mol)
RETURNS int4
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_lt(mol,mol)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_le(mol,mol)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_eq(mol,mol)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_ge(mol,mol)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_gt(mol,mol)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION mol_ne(mol,mol)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR < (
        LEFTARG = mol,
    RIGHTARG = mol,
    PROCEDURE = mol_lt,
        COMMUTATOR = '>',
    NEGATOR = '>=',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR <= (
        LEFTARG = mol,
    RIGHTARG = mol,
    PROCEDURE = mol_le,
        COMMUTATOR = '>=',
    NEGATOR = '>',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR >= (
        LEFTARG = mol,
    RIGHTARG = mol,
    PROCEDURE = mol_ge,
        COMMUTATOR = '<=',
    NEGATOR = '<',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR > (
        LEFTARG = mol,
    RIGHTARG = mol,
    PROCEDURE = mol_gt,
        COMMUTATOR = '<',
    NEGATOR = '<=',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR = (
        LEFTARG = mol,
    RIGHTARG = mol,
    PROCEDURE = mol_eq,
        COMMUTATOR = '=',
    NEGATOR = '<>',
        RESTRICT = eqsel,
    JOIN = eqjoinsel,
        SORT1 = '<',
    SORT2 = '<'
);

CREATE OPERATOR <> (
        LEFTARG = mol,
    RIGHTARG = mol,
    PROCEDURE = mol_ne,
        COMMUTATOR = '<>',
    NEGATOR = '=',
        RESTRICT = neqsel,
    JOIN = neqjoinsel
);

CREATE OPERATOR @= (
    LEFTARG = mol,
    RIGHTARG = mol,
    PROCEDURE = mol_eq,
    COMMUTATOR = '=',
    NEGATOR = '<>',
    RESTRICT = eqsel,
    JOIN = eqjoinsel
);

CREATE OR REPLACE FUNCTION bfp_cmp(bfp,bfp)
RETURNS int4
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_lt(bfp,bfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_le(bfp,bfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_eq(bfp,bfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_ge(bfp,bfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_gt(bfp,bfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION bfp_ne(bfp,bfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR < (
        LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = bfp_lt,
        COMMUTATOR = '>',
    NEGATOR = '>=',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR <= (
        LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = bfp_le,
        COMMUTATOR = '>=',
    NEGATOR = '>',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR >= (
        LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = bfp_ge,
        COMMUTATOR = '<=',
    NEGATOR = '<',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR > (
        LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = bfp_gt,
        COMMUTATOR = '<',
    NEGATOR = '<=',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR = (
        LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = bfp_eq,
        COMMUTATOR = '=',
    NEGATOR = '<>',
        RESTRICT = eqsel,
    JOIN = eqjoinsel,
        SORT1 = '<',
    SORT2 = '<'
);

CREATE OPERATOR <> (
        LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = bfp_ne,
        COMMUTATOR = '<>',
    NEGATOR = '=',
        RESTRICT = neqsel,
    JOIN = neqjoinsel
);

CREATE OR REPLACE FUNCTION sfp_cmp(sfp,sfp)
RETURNS int4
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION sfp_lt(sfp,sfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION sfp_le(sfp,sfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION sfp_eq(sfp,sfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION sfp_ge(sfp,sfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION sfp_gt(sfp,sfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION sfp_ne(sfp,sfp)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR < (
        LEFTARG = sfp,
    RIGHTARG = sfp,
    PROCEDURE = sfp_lt,
        COMMUTATOR = '>',
    NEGATOR = '>=',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR <= (
        LEFTARG = sfp,
    RIGHTARG = sfp,
    PROCEDURE = sfp_le,
        COMMUTATOR = '>=',
    NEGATOR = '>',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR >= (
        LEFTARG = sfp,
    RIGHTARG = sfp,
    PROCEDURE = sfp_ge,
        COMMUTATOR = '<=',
    NEGATOR = '<',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR > (
        LEFTARG = sfp,
    RIGHTARG = sfp,
    PROCEDURE = sfp_gt,
        COMMUTATOR = '<',
    NEGATOR = '<=',
        RESTRICT = contsel,
    JOIN = contjoinsel
);

CREATE OPERATOR = (
        LEFTARG = sfp,
    RIGHTARG = sfp,
    PROCEDURE = sfp_eq,
        COMMUTATOR = '=',
    NEGATOR = '<>',
        RESTRICT = eqsel,
    JOIN = eqjoinsel,
        SORT1 = '<',
    SORT2 = '<'
);

CREATE OPERATOR <> (
        LEFTARG = sfp,
    RIGHTARG = sfp,
    PROCEDURE = sfp_ne,
        COMMUTATOR = '<>',
    NEGATOR = '=',
        RESTRICT = neqsel,
    JOIN = neqjoinsel
);

CREATE OPERATOR CLASS mol_ops
    DEFAULT FOR TYPE mol USING btree AS
        OPERATOR        1       < ,
        OPERATOR        2       <= ,
        OPERATOR        3       = ,
        OPERATOR        4       >= ,
        OPERATOR        5       > ,
        FUNCTION        1       mol_cmp(mol, mol);

CREATE OPERATOR CLASS bfp_ops
    DEFAULT FOR TYPE bfp USING btree AS
        OPERATOR        1       < ,
        OPERATOR        2       <= ,
        OPERATOR        3       = ,
        OPERATOR        4       >= ,
        OPERATOR        5       > ,
        FUNCTION        1       bfp_cmp(bfp, bfp);

CREATE OPERATOR CLASS sfp_ops
    DEFAULT FOR TYPE sfp USING btree AS
        OPERATOR        1       < ,
        OPERATOR        2       <= ,
        OPERATOR        3       = ,
        OPERATOR        4       >= ,
        OPERATOR        5       > ,
        FUNCTION        1       sfp_cmp(sfp, sfp);

CREATE OPERATOR CLASS mol_ops
    DEFAULT FOR TYPE mol USING hash AS
	OPERATOR	1	=,
	FUNCTION	1	hashvarlena(internal);

CREATE OPERATOR CLASS bfp_ops
    DEFAULT FOR TYPE bfp USING hash AS
	OPERATOR	1	=,
	FUNCTION	1	hashvarlena(internal);

CREATE OPERATOR CLASS sfp_ops
    DEFAULT FOR TYPE sfp USING hash AS
	OPERATOR	1	=,
	FUNCTION	1	hashvarlena(internal);

CREATE OR REPLACE FUNCTION gmol_consistent(bytea,internal,int4)
    RETURNS bool
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gmol_compress(internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gmol_decompress(internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gmol_penalty(internal,internal,internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION gmol_picksplit(internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gmol_union(bytea, internal)
    RETURNS _int4
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gmol_same(bytea, bytea, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OPERATOR CLASS mol_ops
DEFAULT FOR TYPE mol USING gist
AS
	OPERATOR	3	@> (mol, mol),
	OPERATOR	4	<@ (mol, mol),
	OPERATOR	3	@> (mol, qmol),
	OPERATOR	4	<@ (qmol, mol),
	OPERATOR	6	@= (mol, mol),
    FUNCTION    1   gmol_consistent (bytea, internal, int4),
    FUNCTION    2   gmol_union (bytea, internal),
    FUNCTION    3   gmol_compress (internal),
    FUNCTION    4   gmol_decompress (internal),
    FUNCTION    5   gmol_penalty (internal, internal, internal),
    FUNCTION    6   gmol_picksplit (internal, internal),
    FUNCTION    7   gmol_same (bytea, bytea, internal),
STORAGE         bytea;

CREATE OR REPLACE FUNCTION gbfp_consistent(bytea,internal,int4)
    RETURNS bool
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gbfp_distance(internal, bytea, smallint, oid)
RETURNS float8
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION gbfp_compress(internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gbfp_penalty(internal,internal,internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR <%> (
    LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = tanimoto_dist(bfp, bfp),
    COMMUTATOR = '<%>'
);

CREATE OPERATOR <#> (
    LEFTARG = bfp,
    RIGHTARG = bfp,
    PROCEDURE = dice_dist(bfp, bfp),
    COMMUTATOR = '<#>'
);

CREATE OPERATOR CLASS bfp_ops
DEFAULT FOR TYPE bfp USING gist
AS
    OPERATOR    1   % (bfp, bfp),
    OPERATOR    2   # (bfp, bfp),
    OPERATOR    3   <%> FOR ORDER BY pg_catalog.float_ops,
    OPERATOR    4   <#> FOR ORDER BY pg_catalog.float_ops,
    FUNCTION    1   gbfp_consistent (bytea, internal, int4),
    FUNCTION    2   gmol_union (bytea, internal),
    FUNCTION    3   gbfp_compress (internal),
    FUNCTION    4   gmol_decompress (internal),
    FUNCTION    5   gbfp_penalty (internal, internal, internal),
    FUNCTION    6   gmol_picksplit (internal, internal),
    FUNCTION    7   gmol_same (bytea, bytea, internal),
    FUNCTION    8   (bfp, bfp) gbfp_distance(internal, bytea, smallint, oid),
STORAGE         bytea;

CREATE OR REPLACE FUNCTION gsfp_consistent(bytea,internal,int4)
    RETURNS bool
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gsfp_compress(internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OPERATOR CLASS sfp_ops
DEFAULT FOR TYPE sfp USING gist
AS
    OPERATOR    1   % (sfp, sfp),
    OPERATOR    2   # (sfp, sfp),
    FUNCTION    1   gsfp_consistent (bytea, internal, int4),
    FUNCTION    2   gmol_union (bytea, internal),
    FUNCTION    3   gsfp_compress (internal),
    FUNCTION    4   gmol_decompress (internal),
    FUNCTION    5   gmol_penalty (internal, internal, internal),
    FUNCTION    6   gmol_picksplit (internal, internal),
    FUNCTION    7   gmol_same (bytea, bytea, internal),
STORAGE         bytea;

CREATE OR REPLACE FUNCTION gslfp_consistent(bytea,internal,int4)
    RETURNS bool
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gslfp_compress(internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gslfp_decompress(internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gslfp_penalty(internal,internal,internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION gslfp_picksplit(internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gslfp_union(bytea, internal)
    RETURNS _int4
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION gslfp_same(bytea, bytea, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OPERATOR CLASS sfp_low_ops
FOR TYPE sfp USING gist
AS
    OPERATOR    1   % (sfp, sfp),
    OPERATOR    2   # (sfp, sfp),
    FUNCTION    1   gslfp_consistent (bytea, internal, int4),
    FUNCTION    2   gslfp_union (bytea, internal),
    FUNCTION    3   gslfp_compress (internal),
    FUNCTION    4   gslfp_decompress (internal),
    FUNCTION    5   gslfp_penalty (internal, internal, internal),
    FUNCTION    6   gslfp_picksplit (internal, internal),
    FUNCTION    7   gslfp_same (bytea, bytea, internal),
STORAGE         bytea;



CREATE OR REPLACE FUNCTION reaction_in(cstring)
RETURNS reaction
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_out(reaction)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_recv(internal)
RETURNS reaction
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_send(reaction)
RETURNS bytea
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE TYPE reaction (
        INTERNALLENGTH = -1,
        INPUT = reaction_in,
        OUTPUT = reaction_out,
        RECEIVE = reaction_recv,
        SEND = reaction_send,
        STORAGE = extended
);

CREATE OR REPLACE FUNCTION reaction_from_smiles(cstring)
RETURNS reaction
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_from_smarts(cstring)
RETURNS reaction
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_to_smiles(reaction)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_to_smarts(reaction)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_from_ctab(cstring)
RETURNS reaction
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_to_ctab(reaction)
RETURNS cstring
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_numreactants(reaction)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_numproducts(reaction)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_numagents(reaction)
RETURNS integer
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION substruct(reaction, reaction)
RETURNS bool
AS 'MODULE_PATHNAME', 'reaction_substruct'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR @> (
	LEFTARG = reaction,
	RIGHTARG = reaction,
	PROCEDURE = substruct(reaction, reaction),
	COMMUTATOR = '<@',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION substructFP(reaction, reaction)
RETURNS bool
AS 'MODULE_PATHNAME', 'reaction_substructFP'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR ?> (
	LEFTARG = reaction,
	RIGHTARG = reaction,
	PROCEDURE = substructFP(reaction, reaction),
	COMMUTATOR = '?<',
	RESTRICT = contsel,
	JOIN = contjoinsel
);


CREATE OR REPLACE FUNCTION rsubstruct(reaction, reaction)
RETURNS bool
AS 'MODULE_PATHNAME', 'reaction_rsubstruct'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR <@ (
	LEFTARG = reaction,
	RIGHTARG = reaction,
	PROCEDURE = rsubstruct(reaction, reaction),
	COMMUTATOR = '@>',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION rsubstructFP(reaction, reaction)
RETURNS bool
AS 'MODULE_PATHNAME', 'reaction_rsubstructFP'
LANGUAGE C STRICT IMMUTABLE;


CREATE OPERATOR ?< (
	LEFTARG = reaction,
	RIGHTARG = reaction,
	PROCEDURE = rsubstructFP(reaction, reaction),
	COMMUTATOR = '?>',
	RESTRICT = contsel,
	JOIN = contjoinsel
);

CREATE OR REPLACE FUNCTION reaction_ne(reaction,reaction)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR <> (
        LEFTARG = reaction,
    RIGHTARG = reaction,
    PROCEDURE = reaction_ne,
        COMMUTATOR = '<>',
    NEGATOR = '=',
        RESTRICT = neqsel,
    JOIN = neqjoinsel
);

CREATE OR REPLACE FUNCTION reaction_eq(reaction,reaction)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OPERATOR @= (
    LEFTARG = reaction,
    RIGHTARG = reaction,
    PROCEDURE = reaction_eq,
    COMMUTATOR = '=',
    NEGATOR = '<>',
    RESTRICT = eqsel,
    JOIN = eqjoinsel
);

CREATE OR REPLACE FUNCTION greaction_consistent(bytea,internal,int4)
    RETURNS bool
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION greaction_compress(internal)
    RETURNS internal
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE;

CREATE OPERATOR CLASS reaction_ops
DEFAULT FOR TYPE reaction USING gist
AS
	OPERATOR	3	@> (reaction, reaction),
	OPERATOR	4	<@ (reaction, reaction),
	OPERATOR	6	@= (reaction, reaction),
	OPERATOR	7	?> (reaction, reaction),
	OPERATOR	8	?< (reaction, reaction),	
    FUNCTION    1   greaction_consistent (bytea, internal, int4),
    FUNCTION    2   gmol_union (bytea, internal),
    FUNCTION    3   greaction_compress (internal),
    FUNCTION    4   gmol_decompress (internal),
    FUNCTION    5   gmol_penalty (internal, internal, internal),
    FUNCTION    6   gmol_picksplit (internal, internal),
    FUNCTION    7   gmol_same (bytea, bytea, internal),
STORAGE         bytea;

CREATE OR REPLACE FUNCTION has_reaction_substructmatch(queryreaction char, tablename regclass, columnname text) 
  RETURNS SETOF reaction AS
$BODY$
DECLARE
	nof_all_entries real;
	nof_index_matches real;
	match_ratio real;
BEGIN
  SET enable_seqscan=off;
  SET enable_bitmapscan=on;
  SET enable_indexscan=on;
  RAISE NOTICE 'Your query: %', queryreaction;
  EXECUTE 'SELECT COUNT(*) FROM ' || tablename INTO nof_all_entries;
  RAISE NOTICE 'Number of reactions in table: %', nof_all_entries;
  EXECUTE 'SELECT COUNT(*) FROM ' || tablename || ' WHERE ' || quote_ident(columnname) || '?>' || quote_literal(queryreaction) INTO nof_index_matches;
  RAISE NOTICE 'Number of matched reactions in the index: %', nof_index_matches;
  match_ratio := nof_index_matches/nof_all_entries;
  RAISE NOTICE 'Match ratio: %', match_ratio;
  IF match_ratio > 0.7 THEN
    SET enable_seqscan=on;
    SET enable_bitmapscan=off;
    SET enable_indexscan=off;
    IF match_ratio >= 1.0 THEN
      RAISE NOTICE 'Your query matches % percent of the index. You are sure you already have build an index?', match_ratio*100.0;
    END IF;
    RAISE NOTICE 'Your query matches % percent of the index. Executing strategy: SequentialScan. Starting substructure matching..', match_ratio*100.0;
  ELSE
    SET enable_seqscan=off;
    SET enable_bitmapscan=on;
    SET enable_indexscan=on;
    RAISE NOTICE 'Executing strategy: IndexScan and BitMapHeapScan. % matches have to be rechecked. Starting substructure matching...', nof_index_matches;
  END IF;
  RETURN QUERY EXECUTE 'SELECT * FROM ' || tablename || ' WHERE ' || quote_ident(columnname) || '@>' || quote_literal(queryreaction);
END
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION reaction_difference_fp(reaction, int default 1)
RETURNS sfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION reaction_structural_bfp(reaction,int default 5)
RETURNS bfp
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT IMMUTABLE;
  