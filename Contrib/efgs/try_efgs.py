# -*- coding: utf-8 -*-
"""
Implementation of Ertl's Functional Group Algorithm in RDKit

Ertl, Peter. 

An algorithm to identify functional groups in organic molecules

J Cheminform (2017) 9:36; DOI: 10.1186/s13321-017-0225-z

@author: gonzalo.colmenarejo
"""


from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from chembl_structure_pipeline import standardizer as sdz
from datetime import datetime
from IPython.display import Image
from PIL import Image as pilImage
from PIL import ImageDraw, ImageFont
from rdkit.Chem.Draw import rdMolDraw2D
import efgs
import pandas as pd
import io



ertl_smis = [
    'Cc1nc(NS(=O)(=O)c2ccc(N)cc2)nc(C)c1', 
    'NC(=N)c1ccc(C=Cc2ccc(cc2O)C(=N)N)cc1', 
    'CC(=O)Nc1nnc(s1)S(=O)(=O)N',
    'NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl',
    'CNC1=Nc2ccc(Cl)cc2C(=N(=O)C1)c3ccccc3',
    'Cc1onc(c1C(=O)NC2C3SC(C)(C)C(N3C2=O)C(=O)O)c4ccccc4',
    'Clc1ccccc1C2=NCC(=O)Nc3ccc(cc23)N(=O)=O',
    'COc1cc(cc(C(=O)NCC2CCCN2CC=C)c1OC)S(=O)(=O)N',
    'Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl',
    'Clc1ccc2Oc3ccccc3N=C(N4CCNCC4)c2c1',
    'FC(F)(F)CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13',
    'OCC1OC(CC1O)n2cnc3C(O)CNC=Nc32',
    'CCNC1CC(C)S(=O)(=O)c2sc(cc12)S(=O)(=O)N',
    'CC(O)C1C2C(C)C(=C(N2C1=O)C(=O)O)SC3CNC(C3)C(=O)N(C)C',
    'CC1CN(CC(C)N1)c2c(F)c(N)c3c(=O)c(cn(C4CC4)c3c2F)C(=O)O',
    'CC(=CCC1C(=O)N(N(C1=O)c2ccccc2)c3ccccc3)C',
    'Clc1ccc2N=C3NC(=O)CN3Cc2c1Cl',
    'CC(=O)NC1C(NC(=N)N)C=C(OC1C(O)C(O)CO)C(=O)O',
    'CC(O)C(O)C1CNc2nc(N)nc(O)c2N1',
    'NC1CCCCN(C1)c2c(Cl)cc3c(=O)c(cn(C4CC4)c3c2Cl)C(=O)O',
]



def mindex(mol):
    # Utility function to draw molecule with idx
    for atom in mol.GetAtoms():
        try:
            atom.SetAtomMapNum(atom.GetIdx())
        except:
            None
    return mol


def demindex(mol):
    # Utility function to return mindex-ed molecule to normal
    for atom in mol.GetAtoms():
        atom.ClearProp('molAtomMapNumber')
    return mol


def molimg(x):
    # Utility function to draw image of labeled mol from PNG binary string
    try:
        return pilImage.open(io.BytesIO(x))
    except:
        return None
    

def img_grid(images, titles = None, num_columns = 5,font = "arial.ttf", font_size = 18, title_height = 30, cell_width = 300, cell_height = 300, bg_color = (255, 255, 255)):
    """
    Creates a grid of images with optional titles under each image.

    Parameters:
    - images (list): List of PIL Image objects.
    - titles (list, optional): List of titles as strings for each image. If None, no titles will be added.
    - num_columns (int): Number of columns in the grid.
    - font_sizes: Font size of titles (if provided).
    - title_height: Height of titles (if provieded).
    - cell_width (int): Width of each image cell.
    - cell_height (int): Height of each image cell (not including title space).
    - bg_color (tuple): Background color as an RGB tuple.

    Returns:
    - PIL Image object containing the image grid.
    """
    # Set title height if titles are provided; otherwise, set to zero
    title_height = title_height if titles else 0
    cell_height_with_title = cell_height + title_height + 2

    # Calculate grid dimensions
    num_images = len(images)
    num_rows = (num_images + num_columns - 1) // num_columns  # Round up for last row

    # Create a blank background image for the grid
    grid_width = num_columns * cell_width
    grid_height = num_rows * cell_height_with_title
    grid_image = pilImage.new("RGB", (grid_width, grid_height), bg_color)

    # Set up font for titles if needed
    if titles:
        try:
            font = ImageFont.truetype(font, font_size)  # Change font and size as desired
        except IOError:
            font = ImageFont.load_default()

    # Iterate through each cell in the grid
    for i in range(num_rows * num_columns):
        row = i // num_columns
        col = i % num_columns
        x = col * cell_width
        y = row * cell_height_with_title

        if i < num_images:
            # Resize the image to fit in the cell
            img_resized = images[i].resize((cell_width, cell_height))
            grid_image.paste(img_resized, (x, y))

            # Draw title if provided
            if titles and i < len(titles):
                draw = ImageDraw.Draw(grid_image)
                title = titles[i]

                # Calculate title position and center it
                text_bbox = draw.textbbox((0, 0), title, font=font)
                text_width = text_bbox[2] - text_bbox[0]
                title_x = x + (cell_width - text_width) // 2  # Center title horizontally
                title_y = y + cell_height  # Position directly beneath the image

                # Draw title text
                draw.text((title_x, title_y), title, fill="black", font=font)
        else:
            # Leave blank space if no image for this cell
            blank_cell = pilImage.new("RGB", (cell_width, cell_height), bg_color)
            grid_image.paste(blank_cell, (x, y))

    return grid_image



# ## Mols from file
start_time = datetime.now()
d = pd.read_csv("../data/ch33_inchis.csv", sep = ";", low_memory = False)
d = d[pd.notna(d.inchi)].reset_index(drop = True) # make sure there are no missing inchis
d = d[d.inchi.apply(lambda x: "." not in x.split("/")[1])].reset_index(drop = True) # make sure we have no complex molecules
d.shape # (839063, 1)
d = d.sample(frac=1, random_state=4210).reset_index(drop = True) # scramble the df
n = d.shape[0] # Here one can put an alternative much smaller n 
d = d[d.index.isin(range(n))].reset_index(drop = True) # keep only n
d["mol"] = d.inchi.apply(lambda x: sdz.standardize_mol(Chem.MolFromInchi(x)))
d = d[pd.notna(d.mol)].reset_index(drop = True) # remove entries without molecule object
d.shape[0]  # n
end_time = datetime.now()
print('Duration to analyze {} molecules: {}'.format(n, end_time - start_time))

start_time = datetime.now()
d[["img_text","fgs","fgsmi","fgmol"]] = d.apply(lambda x: pd.Series(efgs.get_dec_fgs(x["mol"])), axis = 1)
end_time = datetime.now()
print('Duration to analyze {} molecules: {}'.format(n, end_time - start_time))


# Generate distribution of FGs
fg_dist = d.explode(["fgsmi","fgmol"]).groupby("fgsmi", as_index = False).apply(lambda x: 
                                                                                pd.DataFrame({"fgsmi": [x["fgsmi"].iloc[0]], 
                                                                                              "n": [x.shape[0]], 
                                                                                              "mol": [x["fgmol"].iloc[0]]})).sort_values("n", ascending = False)

fg_dist.shape  # 5829, 3 (for all ChEMBL comps)
fg_dist.to_csv("../results/fg_dist.csv", sep = ";", index = False)  
# fg_dist = pd.read_csv("../results/fg_dist.csv", sep = ";", index = False) # To recover data
    
mols = fg_dist.mol.iloc[:64].tolist()
legs = fg_dist.fgsmi.iloc[:64].tolist()


## High-resolution grid with top-48 FGs 
images = []
img_width = 250
img_height = 250
for m in mols[:48]:
    Chem.rdDepictor.Compute2DCoords(m)
    Chem.rdDepictor.StraightenDepiction(m)
    d2d = rdMolDraw2D.MolDraw2DCairo(img_width,img_height)
    dopts = d2d.drawOptions()
    dopts.bondLineWidth = 5
    dopts.padding = 0.2
    d2d.DrawMolecule(m)
    d2d.FinishDrawing()
    images.append(molimg(d2d.GetDrawingText()))
   
fgs_grid = img_grid(images[:48], titles = legs[:48], num_columns = 8, font_size = 28, cell_width = 300, cell_height = 300, bg_color=(255, 255, 255))
fgs_grid = fgs_grid.resize((6000, int(6000*fgs_grid.height/fgs_grid.width)))  
fgs_grid.save("../results/fg_grid.png", dpi = (300,300))


# High-resolution grid with 20 first molecules with colored functional groups
images = [molimg(x) for x in d.img_text.iloc[0:20].tolist()]
cell_width = images[0].width
cell_height = images[0].height
molfg_grid = img_grid(images, num_columns = 5, cell_width = cell_width, cell_height = cell_height, bg_color=(255, 255, 255))
molfg_grid
molfg_grid = molfg_grid.resize((6000, int(6000*molfg_grid.height/molfg_grid.width)))
molfg_grid.save("../RESULTS/molfg_grid.png", dpi = (300,300))


# High-resolution grid with 20 Ertl molecules in Figure 1
ertmols = [Chem.MolFromSmiles(s) for s in ertl_smis[:20]]
images = [molimg(efgs.get_dec_fgs(ertmol)[0]) for ertmol in ertmols] 
cell_width = images[0].width
cell_height = images[0].height
molfg_grid = img_grid(images, num_columns = 5, cell_width = cell_width, cell_height = cell_height, bg_color=(255, 255, 255))
molfg_grid
molfg_grid = molfg_grid.resize((6000, int(6000*molfg_grid.height/molfg_grid.width)))
molfg_grid.save("../RESULTS/ertl_molfg_grid.png", dpi = (300,300))
