
"""RDKit Conformer Browser
Derived by Malitha Humayun Kabir from Greg Landrum's code as a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck's suggestions helped to improve many lines of codes
Date: 31th July 2017
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


# I think this function is not required
# User should be able to supply appropriate input
# A notebook showing 3Dmol visualization should be present in rdkit doc
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
    
    
class MolViewState(object):
    def __init__(self, molecules, uid):
        """ Molecules is dictionary of molecules """
        if uid is None:
            self.uid = str(time.time()).replace('.','')
        else:
            self.uid = uid
            
        self.moldict = molecules
        
        # These should have reasonable initial values
        self.rdkit_mol_select = set()
        self.rdkit_conf_select = set()
        self.allConfIDs = []
        
    def selectMolecules(self, selectAllMols, selectMultiMols, selectMol):
        """ Select either all moleculs or add selectMol or show only selectMol """
        if selectAllMols:
            self.rdkit_mol_select = set(self.moldict)
        elif selectMultiMols:
            self.rdkit_mol_select.add(selectMol)
        else:
            self.rdkit_mol_select = {selectMol}
        
    def selectConformations(self, selectAllConfs, selectMultiConfs, selectConf):
        """ For all selected molecules, select either all conformations or add selectConf or show only selectConf """
        for mol in self.selectedMolecules:
            nConformers = mol.GetNumConformers()
            if nConformers > 1:
                if selectAllConfs:
                    self.rdkit_conf_select = set(range(nConformers))
                elif selectMultiConfs:
                    self.rdkit_conf_select.add(selectConf)
                else:
                    self.rdkit_conf_select = {selectConf}
            elif mol.GetNumConformers() == 1:
                self.rdkit_conf_select = {0}
    
    @property
    def selectedMolNames(self):
        """ Return the names of all selected molecules """
        return self.rdkit_mol_select
    
    @property
    def selectedConfIds(self):
        """ Return the names of all selected molecules """
        return self.rdkit_conf_select
    
    @property
    def selectedMolecules(self):
        """ Return the selected molecules """
        return [self.moldict[name] for name in self.selectedMolNames]
    
    @property
    def selectedModels(self):
        """ Iterator over all selected models (molecules/conformations) """
        for mol in self.selectedMolecules:
            for confId in self.selectedConfIds:
                yield Chem.MolToMolBlock(mol, confId=confId)
                
    @property
    def allConfIds(self):
        """ Return the number of conformations - use the first selected molecule to determine """
        nconfIds = self.selectedMolecules[0].GetNumConformers()
        return list(range(nconfIds))
    @property
    def getPropPrecalculated(self):
        """ Return the precalculated properties """
        if len(self.selectedMolNames) == 1:
            mol = self.moldict[list(self.selectedMolNames)[0]]
            if len(mol.GetPropNames()) == 0:
                return 'Not found'
            else:
                prop = [(prop,mol.GetProp(prop)) for prop in mol.GetPropNames()]
                return dict(prop)
        else:
            return None
    @property
    def getPropRDKit(self):
        """ Return the RDKit calculated properties """
        if len(self.selectedMolNames) == 1:
            mol = self.moldict[list(self.selectedMolNames)[0]]
            prop = list()
            for prop_name in PROP_RDKIT:
                calculator = Descriptors.__dict__[prop_name]
                prop.append((prop_name,calculator(mol)))
            return dict(prop)
        else:
            return None
    
    
def handleMolAndConf(molViewState):
    """ Handles molecules and conformers in update3D function """
    uid = molViewState.uid
    molViewState.selectMolecules(globals()['selectAllMols_'+uid].value, 
                                 globals()['selectMultiMols_'+uid].value,
                                 globals()['molId_'+uid].value)
    sMolNames = molViewState.selectedMolNames
    globals()['selected_mols_view_'+uid].value = ', '.join(sMolNames)
    molViewState.selectConformations(globals()['selectAllConfs_'+uid].value, 
                                     globals()['selectMultiConfs_'+uid].value,
                                     globals()['confId_'+uid].value)
    globals()['confId_'+uid].options = molViewState.allConfIds
    sConfIds = molViewState.selectedConfIds
    globals()['selected_confs_view_'+uid].value = ', '.join([str(x) for x in sConfIds])
    return molViewState
    
def handleProperty(molViewState):
    """ Handles property in update3D function """
    # For retrieving precalculated property
    uid = molViewState.uid
    preCalcProp = molViewState.getPropPrecalculated
    if isinstance(preCalcProp, dict):
        try:
            globals()['prop_precalc_wg_'+uid].options = preCalcProp.keys()
            prop=globals()['prop_precalc_wg_'+uid].value
            globals()['prop_precalc_view_'+uid].value = prop + ' : ' + str(preCalcProp[prop])
        except:
            pass
    elif preCalcProp == 'Not found':
        try:
            globals()['prop_precalc_view_'+uid].value = 'No precalculated property found!'
        except:
            pass
    else:
        try:
            globals()['prop_precalc_view_'+uid].value = 'Single molecule selection required!'
        except:
            pass
        
    # For retrieving rdkit calculated property
    propRDKit = molViewState.getPropRDKit
    if isinstance(propRDKit, dict):
        try:
            prop=globals()['prop_calc_wg_'+uid].value
            globals()['prop_calc_view_'+uid].value = prop + ' : ' + str(propRDKit[prop])
        except:
            pass
    else:
        try:
            globals()['prop_calc_view_'+uid].value = 'Single molecule selection required!'
        except:
            pass
    return 'done'
    
def handleLigandStyle(uid, view):
    """ Handles ligand drawing style and color in update3D function """
    try:
        drawAs=globals()['drawAs_wg_'+uid].value
    except:
        drawAs = globals()['drawAs_no_wg_'+uid]
    
    try:
        color=globals()['colorScheme_'+uid].value
    except:
        color = 'default'
        
    if drawAs == 'surface':
        view.addSurface('SES', {});
    elif drawAs == 'ballstick':
        view.setStyle({},{'stick':{'radius':'0.2','colorscheme': color},
                          'sphere':{'radius':'0.4', 'colorscheme': color}
                         }
                     );
    else:
        view.setStyle({},{drawAs:{'colorscheme': color}})
            
    return view
    

def handleLigandLabeling(molViewState, view):
    """ Handles ligand labeling (conf label and atom label) in update3D function """
    uid = molViewState.uid
    sMolNames = molViewState.selectedMolNames
    sConfIds = molViewState.selectedConfIds
    if len(sMolNames)==1 and len(sConfIds)==1:
        mol = molViewState.moldict[list(sMolNames)[0]]
        confId=globals()['confId_'+uid].value
        sconf = mol.GetConformer(confId)
        xyz = sconf.GetPositions()
        try:
            if globals()['confLabel_'+uid].value:
                label = list(sMolNames)[0] + ':' + str(confId)
                view.addLabel(label, {'backgroundColor':'gray', 'fontColor':'white',
                                      'showBackground':'true', 'alignment':'bottomCenter'})
            if globals()['atomLabel_'+uid].value:
                labels=list()
                for i in list(range(sconf.GetNumAtoms())):
                    perAtomLabel = mol.GetAtomWithIdx(i).GetSymbol()+str(mol.GetAtomWithIdx(i).GetIdx())
                    labels.append(perAtomLabel)
                
                for j in range(sconf.GetNumAtoms()):
                    view.addLabel(labels[j], {'inFront' : 'false', 
                                              'fontSize' : '12',
                                              'fontColor':'gray',
                                              'showBackground':'false',
                                              'position' : {'x' : xyz[j][0], 'y' : xyz[j][1], 'z' : xyz[j][2]}
                                             })
        except:
            pass
    return view
    
def handleProteinStyle(uid, view, pModelNum):
    """ Handles protein drawing style in update3D function """
    try:
        pStyle=globals()['pStyle_wg_'+uid].value
    except:
        pStyle = globals()['pStyle_no_wg_'+uid]
        
    if pStyle == 'cartoon':
        try:
            # helicesAsTubes variable is possible only if useDrawAs is True
            helicesAsTubes = globals()['pStyle_tube_wg_'+uid].value
            if helicesAsTubes:
                view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum',
                                                              'arrows': 'true', 'tubes' : 'true'}})
            else:
                view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum', 'arrows': 'true'}})
        except:
            view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum', 'arrows': 'true'}})
    elif pStyle == 'surface':
        view.addSurface('SES', {'model':pModelNum});
    elif pStyle == 'line':
        view.setStyle({'model':pModelNum},{'line':{}});
    return view
    
    
def update3D(model_id):
    """ This function invoked whenever user interacts with widgets.
    It runs first time through handle_button() when the start button clicked.
    """
    uid=globals()['ipy_wgs'][model_id]
    
    if globals()['rdkit_3dviewer_start_button_clicked_'+uid]:
        
        view = globals()['view_'+uid]
        view.removeAllModels()
        view.removeAllSurfaces()
        view.removeAllLabels()
        
        uid = globals()['ipy_wgs'][model_id]
        molViewState  = globals()['mol_views'][uid]
        
        
        # handling mol and conf
        molViewState = handleMolAndConf(molViewState)
        
        # handling properties
        handleProperty(molViewState)
        
        # Add models (molecules/conformations) to viewer
        pModelNum = 0
        for model in molViewState.selectedModels:
            pModelNum = pModelNum+1
            globals()['view_'+uid].addModel(model, 'sdf')
            
        # Ligand Drawing style and colors
        view=handleLigandStyle(uid, view)
        # Ligand labeling
        view = handleLigandLabeling(molViewState, view)
        
        # Add model (protein) to viewer
        if 'protein_'+uid in globals():
            # Show protein if visibitity is true
            if globals()['proteinVisible_'+uid].value:
                pdb = Chem.MolToPDBBlock(globals()['protein_'+uid])
                view.addModel(pdb,'pdb')
                # Protein Drawing style and colors
                view = handleProteinStyle(uid, view, pModelNum)
                
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
    uid=globals()['ipy_wgs'][b._model_id]
    globals()['rdkit_3dviewer_start_button_clicked_'+uid] = True
    update3D(b._model_id)
    
def handle_zoomTo_button(b):
    """This function handles zoomTo button"""
    uid=globals()['ipy_wgs'][b._model_id]
    globals()['view_'+uid].zoomTo()
    display(globals()['view_'+uid].update())
    
    
def ChangeActiveLigand(uid, molId, confId, keepExistingModels = True):
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
                     drawAs='stick', pStyle='cartoon',
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
    
    if 'ipy_wgs' not in globals():
        globals()['ipy_wgs'] = dict()
        
    rdkitWG = globals()['ipy_wgs']
    
    if 'mol_views' not in globals():
        globals()['mol_views'] = dict()
        
    molViewState = MolViewState(moldict, uid)
    globals()['mol_views'][molViewState.uid] = molViewState
    
    uid = molViewState.uid
    
    # Required global objects
    globals()['rdkit_3dviewer_start_button_clicked_'+uid] = False
    globals()['drawAs_no_wg_'+uid] = drawAs
    globals()['pStyle_no_wg_'+uid] = pStyle
    
    keys=molViewState.moldict.keys()
    
    # Right hand panel (widgets)
    
    wgListBox=list()
    
    itemLayout=Layout(display='flex', flex_flow='row', justify_content='space-between')
    
    wgListBox.append(Box([Label(value='uid'), HTML(description='', value=str(uid))], layout=itemLayout))
    
    # Add protein
    if protein is not None:
        #print 'adding protein'
        globals()['protein_'+uid] = protein
    
    
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
    
    
    
    rdkitWG[globals()['molId_'+uid]._model_id] = uid
    rdkitWG[globals()['selectMultiMols_'+uid]._model_id] = uid
    rdkitWG[globals()['selectAllMols_'+uid]._model_id] = uid
    rdkitWG[globals()['confId_'+uid]._model_id] = uid
    rdkitWG[globals()['selectMultiConfs_'+uid]._model_id] = uid
    rdkitWG[globals()['selectAllConfs_'+uid]._model_id] = uid
    
    
    globals()['molId_'+uid].observe(handle_change, names='value')
    globals()['selectMultiMols_'+uid].observe(handle_change, names='value')
    globals()['selectAllMols_'+uid].observe(handle_change, names='value')
    globals()['confId_'+uid].observe(handle_change, names='value')
    globals()['selectMultiConfs_'+uid].observe(handle_change, names='value')
    globals()['selectAllConfs_'+uid].observe(handle_change, names='value')
    
    if protein is not None:
        globals()['proteinVisible_'+uid] = Checkbox(description='proteinVisible', value=True)
        wgListBox.append(Box([Label(value=''),globals()['proteinVisible_'+uid]], layout=itemLayout))
        rdkitWG[globals()['proteinVisible_'+uid]._model_id] = uid
        globals()['proteinVisible_'+uid].observe(handle_change, names='value')
    
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
        
        rdkitWG[globals()['prop_calc_wg_'+uid]._model_id] = uid
        rdkitWG[globals()['prop_precalc_wg_'+uid]._model_id] = uid
        globals()['prop_calc_wg_'+uid].observe(handle_change, names='value')
        globals()['prop_precalc_wg_'+uid].observe(handle_change, names='value')
        
        
    # useDrawAs
    if useDrawAs:
        
        if moldict is not None:
            globals()['drawAs_wg_'+uid] = Dropdown(description='', options=DRAWING_LIGAND_3D, value=drawAs)
            wgListBox.append(Box([Label(value='drawAs'),globals()['drawAs_wg_'+uid]], layout=itemLayout))
            rdkitWG[globals()['drawAs_wg_'+uid]._model_id] = uid
            globals()['drawAs_wg_'+uid].observe(handle_change, names='value')

        if protein is not None:
            globals()['pStyle_wg_'+uid] = Dropdown(description='', options=DRAWING_PROTEIN_3D, value=pStyle)
            wgListBox.append(Box([Label(value='pStyle'),globals()['pStyle_wg_'+uid]], layout=itemLayout))
            rdkitWG[globals()['pStyle_wg_'+uid]._model_id] = uid
            globals()['pStyle_wg_'+uid].observe(handle_change, names='value')
            
            globals()['pStyle_tube_wg_'+uid] = Checkbox(description='helicesAsTubes', value=False)
            wgListBox.append(Box([Label(value=''),globals()['pStyle_tube_wg_'+uid]], layout=itemLayout))
            rdkitWG[globals()['pStyle_tube_wg_'+uid]._model_id] = uid
            globals()['pStyle_tube_wg_'+uid].observe(handle_change, names='value')
    
    
    # colors 
    if colorPanel is True:
        globals()['colorScheme_'+uid] = Dropdown(description='', options=COLOR_SCHEME_3D, value='default')
        wgListBox.append(Box([Label(value='ligand color'),globals()['colorScheme_'+uid]], layout=itemLayout))
        rdkitWG[globals()['colorScheme_'+uid]._model_id] = uid
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
        
        rdkitWG[globals()['confLabel_'+uid]._model_id] = uid
        rdkitWG[globals()['atomLabel_'+uid]._model_id] = uid
        globals()['confLabel_'+uid].observe(handle_change, names='value')
        globals()['atomLabel_'+uid].observe(handle_change, names='value')
    
    
    # background
    globals()['background_'+uid] = Dropdown(description='', options= BGCOLORS_3D, value=BGCOLORS_3D[0])
    wgListBox.append(Box([Label(value='background'),globals()['background_'+uid]], layout=itemLayout))
    rdkitWG[globals()['background_'+uid]._model_id] = uid
    globals()['background_'+uid].observe(handle_change, names='value')
    
    # buttons
    globals()['start_'+uid] = Button(description="Start!", button_style='success')
    globals()['zoomTo_'+uid] = Button(description="zoomTo", button_style='success')
    buttons=Box([Label(value=''),HBox([globals()['zoomTo_'+uid], globals()['start_'+uid]])],layout=itemLayout)
    rdkitWG[globals()['start_'+uid]._model_id] = uid
    rdkitWG[globals()['zoomTo_'+uid]._model_id] = uid
    globals()['start_'+uid].on_click(handle_start_button)
    globals()['zoomTo_'+uid].on_click(handle_zoomTo_button)
    
    wgListBox.append(buttons)
    
    globals()['ipy_wgs'] = rdkitWG
    
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
    