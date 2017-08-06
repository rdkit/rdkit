
"""RDKit Conformer Browser
Derived by Malitha Humayun Kabir from Greg Landrum's code as a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck's suggestions helped to improve many lines of codes
Date: 6th August 2017
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
import sys

##### Visualization

BGCOLORS_3D = ('0x000000', '0xeeeeee', '0xffffff')

PROP_RDKIT = tuple(sorted(prop for prop, _ in Descriptors._descList))

DRAWING_LIGAND_3D=('line', 'cross', 'stick', 'sphere', 'surface', 'ballstick')

DRAWING_PROTEIN_3D=('line', 'cartoon', 'surface')

LIGAND_COLOR_SCHEME_3D=('default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon',
                        'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                        'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                        'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop')

REPRESENTATION_DEFAULT = {'ligand' : 'stick', 'protein' : 'cartoon'}

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
    def __init__(self, uid, molecules, protein):
        """ Molecules is dictionary of molecules """
        self.uid = uid
        self.moldict = molecules
        self.protein = protein
        
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
    wgData = globals()['wgData'][uid]
    molViewState.selectMolecules(wgData['selectAllMols'].value, 
                                 wgData['selectMultiMols'].value,
                                 wgData['molId'].value)
    sMolNames = molViewState.selectedMolNames
    wgData['selected_mols_view'].value = ', '.join(sMolNames)
    molViewState.selectConformations(wgData['selectAllConfs'].value, 
                                     wgData['selectMultiConfs'].value,
                                     wgData['confId'].value)
    wgData['confId'].options = molViewState.allConfIds
    sConfIds = molViewState.selectedConfIds
    wgData['selected_confs_view'].value = ', '.join([str(x) for x in sConfIds])
    return molViewState
    
def handleProperty(molViewState):
    """ Handles property in update3D function """
    # For retrieving precalculated property
    uid = molViewState.uid
    wgData = globals()['wgData'][uid]
    preCalcProp = molViewState.getPropPrecalculated
    if 'prop_precalc_wg' in wgData:
        if isinstance(preCalcProp, dict):
            wgData['prop_precalc_wg'].options = list(preCalcProp.keys())
            prop=wgData['prop_precalc_wg'].value
            wgData['prop_precalc_view'].value = prop + ' : ' + str(preCalcProp[prop])
        elif preCalcProp == 'Not found':
            wgData['prop_precalc_view'].value = 'No precalculated property found!'
        else:
            wgData['prop_precalc_view'].value = 'Single molecule selection required!'
        
    # For retrieving rdkit calculated property
    propRDKit = molViewState.getPropRDKit
    if 'prop_calc_wg' in wgData:
        if isinstance(propRDKit, dict):
            prop=wgData['prop_calc_wg'].value
            wgData['prop_calc_view'].value = prop + ' : ' + str(propRDKit[prop])
        else:
            wgData['prop_calc_view'].value = 'Single molecule selection required!'
            
    return 'done'
    
def handleLigandStyle(uid, view):
    """ Handles ligand drawing style and color in update3D function """
    
    wgData = globals()['wgData'][uid]
    
    if 'drawAs_wg' in wgData:
        drawAs=wgData['drawAs_wg'].value
    else:
        drawAs = REPRESENTATION_DEFAULT['ligand']
    
    if 'colorScheme' in wgData:
        color=wgData['colorScheme'].value
    else:
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
    wgData = globals()['wgData'][uid]
    sMolNames = molViewState.selectedMolNames
    sConfIds = molViewState.selectedConfIds
    
    if len(sMolNames)==1 and len(sConfIds)==1:
        
        mol = molViewState.moldict[list(sMolNames)[0]]
        confId=wgData['confId'].value
        sconf = mol.GetConformer(confId)
        xyz = sconf.GetPositions()
        
        if 'confLabel' in wgData:
            
            if wgData['confLabel'].value:
                label = list(sMolNames)[0] + ':' + str(confId)
                view.addLabel(label, {'backgroundColor':'gray', 'fontColor':'white',
                                      'showBackground':'true', 'alignment':'bottomCenter'})
        
        if 'atomLabel' in wgData:
            
            if wgData['atomLabel'].value:
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
    
    return view
    
def handleProteinStyle(molViewState, view, pModelNum):
    """ Handles protein drawing style in update3D function """
    uid = molViewState.uid
    wgData = globals()['wgData'][uid]
    
    # 'pStyle_wg' and 'pStyle_tube_wg' rendered as pair
    if 'pStyle_wg' in wgData:
        pStyle=wgData['pStyle_wg'].value
    else:
        pStyle = REPRESENTATION_DEFAULT['protein']
        
    
    if pStyle == 'surface':
        view.addSurface('SES', {'model':pModelNum});
    elif pStyle == 'line':
        view.setStyle({'model':pModelNum},{'line':{}});
    elif pStyle == 'cartoon' and 'pStyle_tube_wg' not in wgData:
        view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum', 'arrows': 'true'}})
    else:
        if wgData['pStyle_tube_wg'].value:
            view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum',
                                                          'arrows': 'true', 'tubes' : 'true'}})
        else:
            view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum', 'arrows': 'true'}})
    
    return view
    
    
def update3D(model_id):
    """ This function invoked whenever user interacts with widgets.
    It runs first time through handle_button() when the start button clicked.
    """
    uid = globals()['wgMIds'][model_id]
    viewInstantiated = globals()['viewInstantiated'][uid]
    wgData = globals()['wgData'][uid]
    molViewState  = globals()['molData'][uid]
    
    if viewInstantiated:
        
        view = wgData['view']
        view.removeAllModels()
        view.removeAllSurfaces()
        view.removeAllLabels()
        
        
        # handling mol and conf
        molViewState = handleMolAndConf(molViewState)
        
        # handling properties
        handleProperty(molViewState)
        
        # Add models (molecules/conformations) to viewer
        pModelNum = 0
        for model in molViewState.selectedModels:
            pModelNum = pModelNum+1
            view.addModel(model, 'sdf')
            
        # Ligand Drawing style and colors
        view=handleLigandStyle(uid, view)
        # Ligand labeling
        view = handleLigandLabeling(molViewState, view)
        
        # Add model (protein) to viewer
        if molViewState.protein is not None:
            # Show protein if visibitity is true
            if 'proteinVisible' in wgData:
                if wgData['proteinVisible'].value:
                    pdb = Chem.MolToPDBBlock(molViewState.protein)
                    view.addModel(pdb,'pdb')
                    view = handleProteinStyle(molViewState, view, pModelNum)
        
        # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
        #view.zoomTo()
        view.setBackgroundColor(wgData['background'].value)
        
        display(view.update())
        
    
def handle_change(change):
    """This function handles all the interactive widgets except buttons and 3Dmol.js viewer"""
    update3D(change.owner.model_id)
    
def handle_start_button(b):
    """This function handles start button."""
    b.icon='check'
    b.description="Done!"
    uid=globals()['wgMIds'][b.model_id]
    globals()['viewInstantiated'][uid] = True
    update3D(b.model_id)
    
def handle_zoomTo_button(b):
    """This function handles zoomTo button"""
    uid=globals()['wgMIds'][b.model_id]
    wgData = globals()['wgData'][uid]
    wgData['view'].zoomTo()
    display(wgData['view'].update())
    
    
def ChangeActiveLigand(uid, molId, confId, removeAllModels = False):
    """This function handles ligand select through python code
    uid : string 
    molId : string 
    confId : intiger
    removeAllModels : logical (default False,  if True - all models removed)
    """
    
    if isinstance(uid, str if sys.version_info[0] >= 3 else basestring) is False:
        raise TypeError("uid must be a string")
        
    try:
        wgData = globals()['wgData'][uid]
    except:
        raise KeyError("invalid uid")
        
    if isinstance(molId, str if sys.version_info[0] >= 3 else basestring) is False:
        raise TypeError("molId must be a string")
        
    if molId not in globals()['molData'][uid].moldict:
        raise KeyError("invalid molId")
        
    # To Do: handle invalid confId
    
    wgData['molId'].value = molId
    wgData['confId'].value = confId
    
    
    if removeAllModels:
        wgData['selectMultiMols'].value=False
        wgData['selectAllMols'].value=False
        wgData['selectMultiConfs'].value=False
        wgData['selectAllConfs'].value=False
    
    update3D(wgData['molId'].model_id)
    
    
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
    
    if uid is None:
        uid = str(time.time()).replace('.','')
        
    if isinstance(uid, str if sys.version_info[0] >= 3 else basestring) is False:
        raise TypeError("uid must be str")
    
    if isinstance(moldict, dict) is False:
        raise TypeError("moldict must be dict")
    
    if all(isinstance(key, str if sys.version_info[0] >= 3 else basestring) for key in moldict.keys()) is False:
        raise TypeError("keys of moldict must be str")
        
    
    if 'wgMIds' in globals():
        wgMIds = globals()['wgMIds']
    else:
        wgMIds = {}
        
    if 'wgData' in globals():
        wgData = globals()['wgData']
        wgData[uid] = {}
    else:
        wgData = {}
        wgData[uid] = {}
        
    if 'molData' in globals():
        molData = globals()['molData']
    else:
        molData = {}
        
    if 'viewInstantiated' in globals():
        viewInstantiated = globals()['viewInstantiated']
    else:
        viewInstantiated = {}
    
    
    viewInstantiated[uid] = False
    
    molViewState = MolViewState(uid, moldict, protein)
    molData[uid] = molViewState
    
    keys = list(molViewState.moldict.keys())
    
    
    # Right hand panel (widgets)
    
    wgBox=list()
    
    itemLayout=Layout(display='flex', flex_flow='row', justify_content='space-between')
    
    wgBox.append(Box([Label(value='uid'), HTML(description='', value=str(uid))], layout=itemLayout))
    
    wgData[uid]['selected_mols_view'] = HTML(description='', value=keys[0])
    wgData[uid]['selected_confs_view'] = HTML(description='', value='0')
    wgData[uid]['molId'] = Dropdown(description='', options=keys,value=keys[0])
    wgData[uid]['selectMultiMols'] = Checkbox(description='selectMultiMols', value=False)
    wgData[uid]['selectAllMols'] = Checkbox(description='selectAllMols', value=False)
    wgData[uid]['confId'] = Dropdown(description='', options=list(range(9)), value=0)
    wgData[uid]['selectMultiConfs'] = Checkbox(description='selectMultiConfs', value=False)
    wgData[uid]['selectAllConfs'] = Checkbox(description='selectAllConfs', value=False)
    
    
    wgBox.append(Box([Label(value='Mols'),wgData[uid]['selected_mols_view']], layout=itemLayout))
    wgBox.append(Box([Label(value='Confs'),wgData[uid]['selected_confs_view']], layout=itemLayout))
    
    wgBox.append(Box([Label(value='molId'),wgData[uid]['molId']], layout=itemLayout))
    
    cbMolSelect=Box([Label(value=''),
                     HBox([wgData[uid]['selectMultiMols'], 
                           wgData[uid]['selectAllMols']])
                    ], layout=itemLayout)
    
    wgBox.append(cbMolSelect)
    
    wgBox.append(Box([Label(value='confId'),wgData[uid]['confId']], layout=itemLayout))
    
    cbConfSelect=Box([Label(value=''),
                      HBox([wgData[uid]['selectMultiConfs'], 
                            wgData[uid]['selectAllConfs']])
                     ], layout=itemLayout)
    
    wgBox.append(cbConfSelect)
    
    wgData[uid]['molId'].observe(handle_change, names='value')
    wgData[uid]['selectMultiMols'].observe(handle_change, names='value')
    wgData[uid]['selectAllMols'].observe(handle_change, names='value')
    wgData[uid]['confId'].observe(handle_change, names='value')
    wgData[uid]['selectMultiConfs'].observe(handle_change, names='value')
    wgData[uid]['selectAllConfs'].observe(handle_change, names='value')
    
    wgMIds[wgData[uid]['molId'].model_id] = uid
    wgMIds[wgData[uid]['selectMultiMols'].model_id] = uid
    wgMIds[wgData[uid]['selectAllMols'].model_id] = uid
    wgMIds[wgData[uid]['confId'].model_id] = uid
    wgMIds[wgData[uid]['selectMultiConfs'].model_id] = uid
    wgMIds[wgData[uid]['selectAllConfs'].model_id] = uid
    
    
    if protein is not None:
        wgData[uid]['proteinVisible'] = Checkbox(description='proteinVisible', value=True)
        wgBox.append(Box([Label(value=''),wgData[uid]['proteinVisible']], layout=itemLayout))
        wgMIds[wgData[uid]['proteinVisible'].model_id] = uid
        wgData[uid]['proteinVisible'].observe(handle_change, names='value')
    
    # prop
    if propPanel is True:
        wgData[uid]['prop_precalc_view'] = HTML(description='', value='initializing...')
        wgData[uid]['prop_precalc_wg'] = Dropdown(description='',options=['select'], value='select')
        wgData[uid]['prop_calc_view'] = HTML(description='', value='initializing...')
        wgData[uid]['prop_calc_wg'] = Dropdown(description='',options=PROP_RDKIT,value='MolLogP')
        
        wgBox.append(Box([Label(value='calc'),wgData[uid]['prop_calc_view']], layout=itemLayout))
        wgBox.append(Box([Label(value='calc'),wgData[uid]['prop_calc_wg']], layout=itemLayout))
        wgBox.append(Box([Label(value='precalc'),wgData[uid]['prop_precalc_view']], layout=itemLayout))
        wgBox.append(Box([Label(value='precalc'),wgData[uid]['prop_precalc_wg']], layout=itemLayout))
        
        wgMIds[wgData[uid]['prop_calc_wg'].model_id] = uid
        wgMIds[wgData[uid]['prop_precalc_wg'].model_id] = uid
        wgData[uid]['prop_calc_wg'].observe(handle_change, names='value')
        wgData[uid]['prop_precalc_wg'].observe(handle_change, names='value')
        
        
    # useDrawAs
    if useDrawAs:
        
        if moldict is not None:
            wgData[uid]['drawAs_wg'] = Dropdown(description='', options=DRAWING_LIGAND_3D, value=drawAs)
            wgBox.append(Box([Label(value='drawAs'),wgData[uid]['drawAs_wg']], layout=itemLayout))
            wgMIds[wgData[uid]['drawAs_wg'].model_id] = uid
            wgData[uid]['drawAs_wg'].observe(handle_change, names='value')

        if protein is not None:
            wgData[uid]['pStyle_wg'] = Dropdown(description='', options=DRAWING_PROTEIN_3D, value=pStyle)
            wgBox.append(Box([Label(value='pStyle'),wgData[uid]['pStyle_wg']], layout=itemLayout))
            wgMIds[wgData[uid]['pStyle_wg'].model_id] = uid
            wgData[uid]['pStyle_wg'].observe(handle_change, names='value')
            
            wgData[uid]['pStyle_tube_wg'] = Checkbox(description='helicesAsTubes', value=False)
            wgBox.append(Box([Label(value=''),wgData[uid]['pStyle_tube_wg']], layout=itemLayout))
            wgMIds[wgData[uid]['pStyle_tube_wg'].model_id] = uid
            wgData[uid]['pStyle_tube_wg'].observe(handle_change, names='value')
    
    
    # colors
    if colorPanel is True:
        wgData[uid]['colorScheme'] = Dropdown(description='', options=LIGAND_COLOR_SCHEME_3D, value='default')
        wgBox.append(Box([Label(value='ligand color'),wgData[uid]['colorScheme']], layout=itemLayout))
        wgMIds[wgData[uid]['colorScheme'].model_id] = uid
        wgData[uid]['colorScheme'].observe(handle_change, names='value')
        
    # labels
    if labelPanel is True:
        wgData[uid]['confLabel'] = Checkbox(description='confLabel', value=False)
        wgData[uid]['atomLabel'] = Checkbox(description='atomLabel', value=False)
        
        cbLabel=Box([Label(value=''),
                     HBox([wgData[uid]['confLabel'], 
                           wgData[uid]['atomLabel']])
                    ], layout=itemLayout)
        
        wgBox.append(cbLabel)
        
        wgMIds[wgData[uid]['confLabel'].model_id] = uid
        wgMIds[wgData[uid]['atomLabel'].model_id] = uid
        wgData[uid]['confLabel'].observe(handle_change, names='value')
        wgData[uid]['atomLabel'].observe(handle_change, names='value')
    
    
    # background
    wgData[uid]['background'] = Dropdown(description='', options= BGCOLORS_3D, value=BGCOLORS_3D[0])
    wgBox.append(Box([Label(value='background'),wgData[uid]['background']], layout=itemLayout))
    wgData[uid]['background'].observe(handle_change, names='value')
    wgMIds[wgData[uid]['background'].model_id] = uid
    
    # buttons
    wgData[uid]['start'] = Button(description="Start!", button_style='success')
    wgData[uid]['zoomTo'] = Button(description="zoomTo", button_style='success')
    
    buttons=Box([Label(value=''),HBox([wgData[uid]['zoomTo'], wgData[uid]['start']])],layout=itemLayout)
    wgBox.append(buttons)
    
    wgData[uid]['start'].on_click(handle_start_button)
    wgData[uid]['zoomTo'].on_click(handle_zoomTo_button)
    wgMIds[wgData[uid]['start'].model_id] = uid
    wgMIds[wgData[uid]['zoomTo'].model_id] = uid
    
    
    # right panel
    wgRightBox=VBox(wgBox, layout=Layout(border='solid', width='45%'))
    
    # left panel (container holding table for viewer)
    size = (435, 485)
    viewerLayout=Layout(width=str(size[0]+5)+'px', height=str(size[1]+5)+'px', border='solid')
    wgLeftWidget = HTML('''<table><tr><td id="%s"></td></tr></table>'''%uid, layout=viewerLayout)
    
    # Combining left and right panel
    viewer = HBox([wgLeftWidget, wgRightBox])
    
    # displaying everything
    display(viewer)
    
    # inserting 3DMol.js viewer in existing container (table)
    wgData[uid]['view'] = py3Dmol.view(width=size[0],height=size[1])
    wgData[uid]['view'].setBackgroundColor('0x000000')
    wgData[uid]['view'].zoomTo()
    
    globals()['wgMIds'] = wgMIds
    globals()['wgData'] = wgData
    globals()['molData'] = molData
    globals()['viewInstantiated'] = viewInstantiated
    
    display(globals()['wgData'][uid]['view'].insert(uid))
    