
"""RDKit Conformer Browser
Derived by Malitha Humayun Kabir from Greg Landrum's code as a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck's suggestions helped to improve many lines of codes
Date: 28th July 2017
Email# malitha12345@gmail.com
"""

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from ipywidgets import Layout, Label, Button, Box, HBox, VBox
from ipywidgets import Dropdown, SelectMultiple, IntSlider, HTML, Checkbox, Button, Text
from IPython.display import display
import time


##### Visualization

BGCOLORS_3D = ('0x000000', '0xeeeeee', '0xffffff')

PROP_RDKIT = tuple(sorted(prop for prop, _ in Descriptors._descList))

DRAWING_LIGAND_3D=('line', 'cross', 'stick', 'cartoon', 'sphere', 'surface', 'ballstick')

DRAWING_PROTEIN_3D=('line', 'cartoon', 'surface')

COLOR_SCHEME_3D=('default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon', 
                 'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                 'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                 'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop')


# I think this function is no longer required
# User should be able to supply appropriate input and notebooks are there to help
def ProcessMolContainingObj(mol):
    """This function checks whether the object type fits the requirements for rendering.
    If the oject doesn't have necessary attributes, it takes action to include that.
    """
    if isinstance(mol, dict):
        moldict = dict(mol)
    elif type(mol) is tuple:
        moldict = {mol[0]: mol[1]}
    elif hasattr(mol, '__iter__') is False:
        moldict = {'0': mol}
    elif type(mol) is list:
        mol_keys=[str(x) for x in range(len(mol))]
        mol_2=[(mol_keys[i],mol[i]) for i in range(len(mol))]
        moldict = dict(mol_2)
    else:
        moldict = None
    return moldict
    
    
    
def update3D(model_id):
    """ This function invoked whenever user interacts with widgets.
    It runs first time through handle_button() when the start button clicked.
    """
    uid=globals()['rdkit_wg_dict'][model_id]
    
    if globals()['rdkit_3dviewer_start_button_clicked_'+uid]:
        
        view = globals()['view_'+uid]
        view.removeAllModels()
        view.removeAllSurfaces()
        view.removeAllLabels()


        molDictKey=globals()['molId_'+uid].value

        if globals()['selectAllMols_'+uid].value:
            globals()['rdkit_mol_selected_'+uid] = set(globals()['moldict_'+uid])
        elif globals()['selectMultiMols_'+uid].value:
            globals()['rdkit_mol_selected_'+uid].add(molDictKey)
        else:
            globals()['rdkit_mol_selected_'+uid] = {molDictKey}


        molNames = set(globals()['rdkit_mol_selected_'+uid])
        if type(molDictKey) != 'str':
            globals()['selected_mols_view_'+uid].value = ', '.join([str(x) for x in molNames])
        else:
            globals()['selected_mols_view_'+uid].value = ', '.join(molNames)
        
        for i in molNames:
            mol = globals()['moldict_'+uid][i]
            
            if mol.GetNumConformers()>1:
                
                allConfIds = list(range(mol.GetNumConformers()))
                globals()['confId_'+uid].options = allConfIds
                
                if globals()['selectAllConfs_'+uid].value:
                    globals()['rdkit_conf_selected_'+uid] = set(allConfIds)
                elif globals()['selectMultiConfs_'+uid].value:
                    globals()['rdkit_conf_selected_'+uid].add(globals()['confId_'+uid].value)
                else:
                    globals()['rdkit_conf_selected_'+uid] = {globals()['confId_'+uid].value}
                    
            elif mol.GetNumConformers()==1:
                globals()['confId_'+uid].options=[0]
                globals()['rdkit_conf_selected_'+uid].add(globals()['confId_'+uid].value)
                globals()['selected_confs_view_'+uid].value = str(globals()['confId_'+uid].value)
                
            
        confIds = set(globals()['rdkit_conf_selected_'+uid])
        globals()['selected_confs_view_'+uid].value = ', '.join([str(x) for x in confIds])
        
        
        # Add molecule to viewer
        LigModelNumber = -1
        for m in molNames:
            for confId in confIds:
                LigModelNumber = LigModelNumber + 1
                mol = globals()['moldict_'+uid][m]
                mb = Chem.MolToMolBlock(mol,confId=confId)
                view.addModel(mb,'sdf')
                
                
                
        if len(molNames)==1:
            mol = globals()['moldict_'+uid][molDictKey]
            # For precalculated property
            try:
                all_prop_from_mol=list(mol.GetPropNames())
                if len(all_prop_from_mol)>0:
                    # update widget
                    globals()['prop_precalc_wg_'+uid].options=all_prop_from_mol
                    # Get selected property
                    prop_name=eval('prop_precalc_wg_'+uid+'.value')
                    prop_value=mol.GetProp(prop_name)
                    prop_str = prop_name + ' : ' + str(prop_value)
                    # Update viewer
                    globals()['prop_precalc_view_'+uid].value = prop_str
                else:
                    globals()['prop_precalc_view_'+uid].value = 'No precalculated property found!'
            except:
                pass
                
            # For calculating rdkit supported property
            try:
                # Descriptor calculation schema eval("Descriptors.TPSA(mol)")
                prop_name=globals()['prop_calc_wg_'+uid].value
                prop_calc_cmd="Descriptors."+ prop_name + "(mol)"
                prop_str = prop_name + ' : ' + str(eval(prop_calc_cmd))
                # Update viewer
                globals()['prop_calc_view_'+uid].value = prop_str
            except:
                pass
        elif len(molNames)>1:
            try:
                globals()['prop_precalc_view_'+uid].value = 'single molecule selection required!'
                globals()['prop_calc_view_'+uid].value = 'single molecule selection required!'
            except:
                pass
                
                
        try:
            drawAs=globals()['drawAs_wg_'+uid].value
        except:
            drawAs = globals()['drawAs_no_wg_'+uid]

        try:
            ligand_color=globals()['colorScheme_'+uid].value
        except:
            ligand_color = 'default'
            
            
        if drawAs == 'surface':
            view.addSurface('SES', {});
        elif drawAs == 'ballstick':
            view.setStyle({},{'stick':{'radius':'0.2','colorscheme': ligand_color},
                              'sphere':{'radius':'0.4', 'colorscheme': ligand_color}});
        else:
            view.setStyle({},{drawAs:{'colorscheme': ligand_color}})
            
            
        if len(molNames)==1 and len(confIds)==1:
            mol = globals()['moldict_'+uid][molDictKey]
            confId=globals()['confId_'+uid].value
            sconf = mol.GetConformer(confId)
            xyz = sconf.GetPositions()
            OwningMol = sconf.GetOwningMol()
            
            try:
                labelConf=globals()['confLabel_'+uid].value
                labelAtom=globals()['atomLabel_'+uid].value
                if labelConf:
                    if type(molDictKey) != 'str':
                        molDictKey = str(molDictKey)
                    label = molDictKey + ':' + str(confId)
                    view.addLabel(label, {'backgroundColor':'gray', 'fontColor':'white',
                                          'showBackground':'true', 'alignment':'bottomCenter'})
                    
                if labelAtom:
                    label_create=[OwningMol.GetAtomWithIdx(i).GetSymbol()+
                                  str(OwningMol.GetAtomWithIdx(i).GetIdx()+1) 
                                  for i in range(sconf.GetNumAtoms())
                                 ]
                    i = None
                    for i in range(sconf.GetNumAtoms()):
                        view.addLabel(label_create[i], {'inFront' : 'false', 
                                                        'fontSize' : '12',
                                                        'fontColor':'gray',
                                                        'showBackground':'false',
                                                        'position' : {'x' : xyz[i][0],
                                                                      'y' : xyz[i][1],
                                                                      'z' : xyz[i][2]
                                                                   }
                                                       })
            except:
                pass
                
        
        # Add protein in viewer
        try:
            pVisibility = globals()['proteinVisible_'+uid].value
            
            if pVisibility:
                try:
                    pStyle=globals()['pStyle_wg_'+uid].value
                except:
                    pStyle = globals()['pStyle_no_wg_'+uid]
                    
                if 'protein_'+uid in globals():
                    pdb = Chem.MolToPDBBlock(globals()['protein_'+uid])
                    globals()['view_'+uid].addModel(pdb,'pdb')
                    
                    if pStyle == 'cartoon':
                        try:
                            # helicesAsTubes variable is possible only if useDrawAs is True
                            helicesAsTubes = globals()['pStyle_tube_wg_'+uid].value
                            if helicesAsTubes:
                                view.setStyle({'model':LigModelNumber+1},{'cartoon':{'color': 'spectrum',
                                                                               'arrows': 'true',
                                                                               'tubes' : 'true'}})
                            else:
                                view.setStyle({'model':LigModelNumber+1},{'cartoon':{'color': 'spectrum',
                                                                               'arrows': 'true'}})
                        except:
                            view.setStyle({'model':LigModelNumber+1},{'cartoon':{'color': 'spectrum',
                                                                                 'arrows': 'true'}})
                    elif pStyle == 'surface':
                        view.addSurface('SES', {'model':LigModelNumber+1});
                    elif pStyle == 'line':
                        view.setStyle({'model':LigModelNumber+1},{'line':{}});
        except:
            pass
        
        
        # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
        #view.zoomTo()
        view.setBackgroundColor(globals()['background_'+uid].value)
        display(view.update())
        
    
def handle_change(change):
    """This function handles all the interactive widgets except buttons and 3Dmol.js viewer"""
    update3D(change.owner._model_id)
    
def handle_start_button(b):
    """This function handles start button."""
    b.icon='check'
    b.description="Done!"
    uid=globals()['rdkit_wg_dict'][b._model_id]
    globals()['rdkit_3dviewer_start_button_clicked_'+uid] = True
    update3D(b._model_id)
    
def handle_zoomTo_button(b):
    """This function handles zoomTo button"""
    uid=globals()['rdkit_wg_dict'][b._model_id]
    globals()['view_'+uid].zoomTo()
    display(globals()['view_'+uid].update())
    
    
def ChangeActiveLigand(uid, molId, confId, keepExistingModels = False):
    """This function handles ligand select through python code
    uid : string 
    molId : string 
    confId : intiger
    keepExistingModels : logical (if False - all models are removed - molId and confId associated model rendered)
    """
    if type(uid) != 'str':
        uid = str(uid)
    if type(molId) != 'str':
        molId = str(molId)
        
    if globals()['molId_'+uid].value != molId:
        globals()['molId_'+uid].value = molId
        
    if globals()['confId_'+uid].value != confId:
        globals()['confId_'+uid].value = confId
    
    if keepExistingModels is False:
        globals()['selectMultiMols_'+uid].value=False
        globals()['selectAllMols_'+uid].value=False
        globals()['selectMultiConfs_'+uid].value=False
        globals()['selectAllConfs_'+uid].value=False
    
    update3D(globals()['molId_'+uid]._model_id)
    
    
    
def ShowConformers3D(uid = None,
                     moldict = None, protein = None,
                     useDrawAs = False, 
                     drawAs=None, pStyle=None,
                     propPanel = False, 
                     colorPanel = False, 
                     labelPanel = False):
    
    """This function initiates required widgets and 3Dmol.js viewer
    uid : an unique id for 3DMol.js viewer, each of the ipywidgets, and all the global objects
    moldict : it is a molecule dictionary where is keys are string
    protein : it is rdkit standard mol object
    drawAs : representation type for moldict or small molecule
    pStyle : protein representation style
    ** the rest of the arguments are panels (user may or may not need those)
    """
    
    # Common suffix for widgets
    if uid is None:
        uid=str(time.time()).replace('.','')
    
    # Required global objects
    globals()['rdkit_mol_selected_'+uid] = set()
    globals()['rdkit_conf_selected_'+uid] = set()
    globals()['rdkit_3dviewer_start_button_clicked_'+uid] = False

    if 'rdkit_wg_dict' not in globals():
        globals()['rdkit_wg_dict'] = dict()
    
    
    # Right hand panel (widgets)
    
    itemLayout=Layout(display='flex', flex_flow='row', justify_content='space-between')
    
    wgListBox=list()
    
    wgListBox.append(Box([Label(value='uid'), HTML(description='', value=str(uid))], layout=itemLayout))
    
    # Add protein
    if protein is not None:
        #print 'adding protein'
        globals()['protein_'+uid] = protein
    
    
    # molecules and conformers
    globals()['moldict_'+uid] = moldict
    keys=list(globals()['moldict_'+uid].keys())
    
    # Users are expected to have keys int. 
    # If keys are not str then still the rendering is possible but ChangeActiveLigand() wont't work currently
    if type(keys[0])=='str': 
        mol_id_ini = keys[0]
    else:
        mol_id_ini = str(keys[0])
        
    globals()['selected_mols_view_'+uid] = HTML(description='', value=mol_id_ini)
    globals()['selected_confs_view_'+uid] = HTML(description='', value='0')
    globals()['molId_'+uid] = Dropdown(description='', options=keys,value=keys[0])
    globals()['selectMultiMols_'+uid] = Checkbox(description='selectMultiMols', value=False)
    globals()['selectAllMols_'+uid] = Checkbox(description='selectAllMols', value=False)
    globals()['confId_'+uid] = Dropdown(description='', options=list(range(9)), value=0)
    globals()['selectMultiConfs_'+uid] = Checkbox(description='selectMultiConfs', value=False)
    globals()['selectAllConfs_'+uid] = Checkbox(description='selectAllConfs', value=False)
    
    
    wgListBox.append(Box([Label(value='Mols'),globals()['selected_mols_view_'+uid]], layout=itemLayout))
    wgListBox.append(Box([Label(value='Confs'),globals()['selected_confs_view_'+uid]], layout=itemLayout))
    
    wgListBox.append(Box([Label(value='molId'),globals()['molId_'+uid]], layout=itemLayout))
    
    cbMolSelect=Box([Label(value=''),
                     HBox([globals()['selectMultiMols_'+uid], 
                           globals()['selectAllMols_'+uid]])
                    ], layout=itemLayout)
    
    wgListBox.append(cbMolSelect)
    
    wgListBox.append(Box([Label(value='confId'),globals()['confId_'+uid]], layout=itemLayout))
    
    cbConfSelect=Box([Label(value=''),
                      HBox([globals()['selectMultiConfs_'+uid], 
                            globals()['selectAllConfs_'+uid]])
                     ], layout=itemLayout)
    
    wgListBox.append(cbConfSelect)
    
    if protein is not None:
        globals()['proteinVisible_'+uid] = Checkbox(description='proteinVisible', value=True)
        wgListBox.append(Box([Label(value=''),globals()['proteinVisible_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['proteinVisible_'+uid]._model_id:uid})
        globals()['proteinVisible_'+uid].observe(handle_change, names='value')
    
    
    globals()['rdkit_wg_dict'].update({globals()['molId_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectMultiMols_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectAllMols_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['confId_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectMultiConfs_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['selectAllConfs_'+uid]._model_id:uid})
    
    
    globals()['molId_'+uid].observe(handle_change, names='value')
    globals()['selectMultiMols_'+uid].observe(handle_change, names='value')
    globals()['selectAllMols_'+uid].observe(handle_change, names='value')
    globals()['confId_'+uid].observe(handle_change, names='value')
    globals()['selectMultiConfs_'+uid].observe(handle_change, names='value')
    globals()['selectAllConfs_'+uid].observe(handle_change, names='value')
    
    
    # prop
    if propPanel is True:
        globals()['prop_precalc_view_'+uid] = HTML(description='', value='initializing...')
        globals()['prop_precalc_wg_'+uid] = Dropdown(description='',options=['select'], value='select')
        globals()['prop_calc_view_'+uid] = HTML(description='', value='initializing...')
        globals()['prop_calc_wg_'+uid] = Dropdown(description='',options=globals()['PROP_RDKIT'],value='MolLogP')
        wgListBox.append(Box([Label(value='calc'),globals()['prop_calc_view_'+uid]], layout=itemLayout))
        wgListBox.append(Box([Label(value='calc'),globals()['prop_calc_wg_'+uid]], layout=itemLayout))
        wgListBox.append(Box([Label(value='precalc'),globals()['prop_precalc_view_'+uid]], layout=itemLayout))
        wgListBox.append(Box([Label(value='precalc'),globals()['prop_precalc_wg_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['prop_calc_wg_'+uid]._model_id:uid})
        globals()['rdkit_wg_dict'].update({globals()['prop_precalc_wg_'+uid]._model_id:uid})
        globals()['prop_calc_wg_'+uid].observe(handle_change, names='value')
        globals()['prop_precalc_wg_'+uid].observe(handle_change, names='value')
        
        
    # drawAs and useDrawAs
    if drawAs is None:
        drawAs = 'stick'
    
    globals()['drawAs_no_wg_'+uid] = drawAs
    
    if pStyle is None:
        pStyle = 'cartoon'
    
    globals()['pStyle_no_wg_'+uid] = pStyle
    
    if useDrawAs:
        
        if moldict is not None:
            globals()['drawAs_wg_'+uid] = Dropdown(description='', options=DRAWING_LIGAND_3D, value=drawAs)
            wgListBox.append(Box([Label(value='drawAs'),globals()['drawAs_wg_'+uid]], layout=itemLayout))
            globals()['rdkit_wg_dict'].update({globals()['drawAs_wg_'+uid]._model_id:uid})
            globals()['drawAs_wg_'+uid].observe(handle_change, names='value')

        if protein is not None:
            globals()['pStyle_wg_'+uid] = Dropdown(description='', options=DRAWING_PROTEIN_3D, value=pStyle)
            wgListBox.append(Box([Label(value='pStyle'),globals()['pStyle_wg_'+uid]], layout=itemLayout))
            globals()['rdkit_wg_dict'].update({globals()['pStyle_wg_'+uid]._model_id:uid})
            globals()['pStyle_wg_'+uid].observe(handle_change, names='value')

            globals()['pStyle_tube_wg_'+uid] = Checkbox(description='helicesAsTubes', value=False)
            wgListBox.append(Box([Label(value=''),globals()['pStyle_tube_wg_'+uid]], layout=itemLayout))
            globals()['rdkit_wg_dict'].update({globals()['pStyle_tube_wg_'+uid]._model_id:uid})
            globals()['pStyle_tube_wg_'+uid].observe(handle_change, names='value')
    
    
    # colors 
    if colorPanel is True:
        globals()['colorScheme_'+uid] = Dropdown(description='', options=COLOR_SCHEME_3D, value='default')
        wgListBox.append(Box([Label(value='ligand color'),globals()['colorScheme_'+uid]], layout=itemLayout))
        globals()['rdkit_wg_dict'].update({globals()['colorScheme_'+uid]._model_id:uid})
        globals()['colorScheme_'+uid].observe(handle_change, names='value')
        
    
    
    # labels
    if labelPanel is True:
        globals()['confLabel_'+uid] = Checkbox(description='confLabel', value=False)
        globals()['atomLabel_'+uid] = Checkbox(description='atomLabel', value=False)
        
        cbLabel=Box([Label(value=''),
                     HBox([globals()['confLabel_'+uid], 
                           globals()['atomLabel_'+uid]])
                    ], layout=itemLayout)
        
        wgListBox.append(cbLabel)
        
        globals()['rdkit_wg_dict'].update({globals()['confLabel_'+uid]._model_id:uid})
        globals()['rdkit_wg_dict'].update({globals()['atomLabel_'+uid]._model_id:uid})
        globals()['confLabel_'+uid].observe(handle_change, names='value')
        globals()['atomLabel_'+uid].observe(handle_change, names='value')
    
    
    # background
    globals()['background_'+uid] = Dropdown(description='', options= BGCOLORS_3D, value=BGCOLORS_3D[0])
    wgListBox.append(Box([Label(value='background'),globals()['background_'+uid]], layout=itemLayout))
    globals()['rdkit_wg_dict'].update({globals()['background_'+uid]._model_id:uid})
    globals()['background_'+uid].observe(handle_change, names='value')
    
    # buttons
    globals()['start_'+uid] = Button(description="Start!", button_style='success')
    globals()['zoomTo_'+uid] = Button(description="zoomTo", button_style='success')
    buttons=Box([Label(value=''),HBox([globals()['zoomTo_'+uid], globals()['start_'+uid]])],layout=itemLayout)
    globals()['rdkit_wg_dict'].update({globals()['start_'+uid]._model_id:uid})
    globals()['rdkit_wg_dict'].update({globals()['zoomTo_'+uid]._model_id:uid})
    globals()['start_'+uid].on_click(handle_start_button)
    globals()['zoomTo_'+uid].on_click(handle_zoomTo_button)
    
    wgListBox.append(buttons)
    
    
    
    # left panel (container holding table for viewer)
    size = (435, 485)
    viewerLayout=Layout(width=str(size[0]+5)+'px', height=str(size[1]+5)+'px', border='solid')
    wgLeftWidget = HTML('''<table><tr><td id="%s"></td></tr></table>'''%uid, layout=viewerLayout)
    
    # right panel
    wgRightBox=VBox(wgListBox, layout=Layout(border='solid', width='45%'))
    
    # Combining left and right panel
    all_wg = HBox([wgLeftWidget, wgRightBox])
    
    # displaying everything
    display(all_wg)
    
    
    # inserting 3DMol.js viewer in existing container (table)
    
    globals()['view_'+uid] = py3Dmol.view(width=size[0],height=size[1])
    globals()['view_'+uid].setBackgroundColor('0x000000')
    globals()['view_'+uid].zoomTo()
    
    display(globals()['view_'+uid].insert(uid))
    