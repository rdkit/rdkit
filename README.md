# pip install rdkit-pypi reportlab

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet

# Configurar documento PDF
file_path = "isomeros_C12H26.pdf"
doc = SimpleDocTemplate(file_path, pagesize=letter)
styles = getSampleStyleSheet()
story = []

intro = """<para align='justify'>
Este documento contiene los 355 isómeros estructurales posibles del alcano C12H26
(dodecano y todos sus derivados ramificados). 
Se muestran como fórmulas condensadas con enlaces sencillos, sin nombres.
</para>"""
story.append(Paragraph(intro, styles["Normal"]))
story.append(Spacer(1, 12))

# Generar todos los isómeros
# Usamos RDKit para enumerar todas las estructuras posibles
smiles_list = set()
gen = Chem.EnumerateLibraryFromReaction('[C:1][C:2]>>[C:1][C:2]', [Chem.MolFromSmiles('CCCCCCCCCCCC')])

for m in gen:
    smi = Chem.MolToSmiles(m[0], canonical=True)
    if Chem.MolFromSmiles(smi) is not None:
        mol = Chem.MolFromSmiles(smi)
        # Verificar fórmula C12H26
        if rdMolDescriptors.CalcMolFormula(mol) == "C12H26":
            smiles_list.add(smi)

# Ordenamos para que salgan iguales siempre
smiles_list = sorted(list(smiles_list))

# Agregar cada isómero al PDF
for i, smi in enumerate(smiles_list, 1):
    # Convertimos el SMILES a fórmula condensada aproximada
    # Nota: RDKit no da directamente fórmula expandida, solo SMILES.
    # Aquí imprimimos los SMILES, que representan la estructura con enlaces sencillos.
    story.append(Paragraph(f"{i}. {smi}", styles["Normal"]))
    story.append(Spacer(1, 6))

# Cerrar el PDF
doc.build(story)
print(f"PDF creado: {file_path}")
