
"""RDKit Conformer Browser
Derived by Malitha Humayun Kabir from Greg Landrum's code as a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck's suggestions helped to improve many lines of codes
Date: 13th August 2017
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

PY3 = sys.version_info[0] >= 3

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
    def __init__(self, uid, molecules, protein, calculateProperties = False):
        """ Molecules is dictionary of molecules """
        if calculateProperties:
            for i in molecules:
                mol = molecules[i]
                for prop_name in PROP_RDKIT:
                    calculator = Descriptors.__dict__[prop_name]
                    mol.SetProp(prop_name, str(calculator(mol)))
                    
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
        if len(self.selectedMolNames) != 1:
            return None
        else:
            mol = self.moldict[list(self.selectedMolNames)[0]]
            if len(mol.GetPropNames()) == 0:
                return 'Not found'
            else:
                return {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}
            
    
    
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
    if 'prop_wg' in wgData:
        if isinstance(preCalcProp, dict):
            wgData['prop_wg'].options = sorted(list(preCalcProp.keys()))
            prop=wgData['prop_wg'].value
            wgData['prop_view'].value = prop + ' : ' + str(preCalcProp[prop])
        elif preCalcProp == 'Not found':
            wgData['prop_view'].value = 'No precalculated property found!'
        else:
            wgData['prop_view'].value = 'Single molecule selection required!'
        
    
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
    
    if not len(sMolNames)==1 and len(sConfIds)==1:
        return None
    else:
        if 'confLabel' in wgData:
            if wgData['confLabel'].value:
                confLabel = list(sMolNames)[0] + ':' + str(wgData['confId'].value)
                view.addLabel(confLabel, {'backgroundColor':'gray', 'fontColor':'white',
                                          'showBackground':'true', 'alignment':'bottomCenter'})
        
        if 'atomLabel' in wgData:
            
            mol = molViewState.moldict[list(sMolNames)[0]]
            confId=wgData['confId'].value
            sconf = mol.GetConformer(confId)
            xyz = sconf.GetPositions()
            
            if wgData['atomLabel'].value:
                for i in list(range(sconf.GetNumAtoms())):
                    atomLabel = '{}{}'.format(mol.GetAtomWithIdx(i).GetSymbol(), i)
                    view.addLabel(atomLabel, {'inFront' : 'false', 
                                              'fontSize' : '12',
                                              'fontColor':'gray',
                                              'showBackground':'false',
                                              'position' : {'x':xyz[i][0], 'y':xyz[i][1], 'z':xyz[i][2]}
                                             })
    
    return view
    
def handleProteinStyle(molViewState, view, pModelNum):
    """ Handles protein drawing style in update3D function """
    uid = molViewState.uid
    wgData = globals()['wgData'][uid]
    
    # 'pStyle_wg' and 'pStyle_helicesAsTubes_wg' rendered as pair
    if 'pStyle_wg' in wgData:
        pStyle=wgData['pStyle_wg'].value
    else:
        pStyle = REPRESENTATION_DEFAULT['protein']
        
    
    if pStyle == 'surface':
        view.addSurface('SES', {'model':pModelNum});
    elif pStyle == 'line':
        view.setStyle({'model':pModelNum},{'line':{}});
    elif pStyle == 'cartoon' and 'pStyle_helicesAsTubes_wg' not in wgData:
        view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum', 'arrows': 'true'}})
    else:
        if wgData['pStyle_helicesAsTubes_wg'].value:
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
        for model in molViewState.selectedModels:
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
                    pModelNum = len(list(molViewState.selectedModels))
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
    
    if isinstance(uid, str if PY3 else basestring) is False:
        raise TypeError("uid must be a string")
        
    try:
        wgData = globals()['wgData'][uid]
    except:
        raise KeyError("invalid uid")
        
    if isinstance(molId, str if PY3 else basestring) is False:
        raise TypeError("molId must be a string")
        
    if molId not in globals()['molData'][uid].moldict:
        raise KeyError("invalid molId")
        
    # TODO: handle invalid confId
    
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
                     propPanel = False, calculateProperties = False,
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
        
    if isinstance(uid, str if PY3 else basestring) is False:
        raise TypeError("uid must be str")
        
    if isinstance(moldict, dict) is False:
        raise TypeError("moldict must be dict")
        
    if all(isinstance(key, str if PY3 else basestring) for key in moldict.keys()) is False:
        raise TypeError("keys of moldict must be str")
        
    wgMIds = globals().setdefault('wgMIds', {})
    
    wgData = globals().setdefault('wgData', {})
    wgData[uid] = {}
    
    molData = globals().setdefault('molData', {})
    
    viewInstantiated = globals().setdefault('viewInstantiated', {})
    viewInstantiated[uid] = False
    
    molViewState = MolViewState(uid, moldict, protein, calculateProperties)
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
    
    widgetsFor3DView = ['molId', 'selectMultiMols', 'selectAllMols', 
                        'confId', 'selectMultiConfs', 'selectAllConfs']
    
    if protein is not None:
        wgData[uid]['proteinVisible'] = Checkbox(description='proteinVisible', value=True)
        wgBox.append(Box([Label(value=''),wgData[uid]['proteinVisible']], layout=itemLayout))
        widgetsFor3DView.append('proteinVisible')
    
    # prop
    if propPanel is True:
        wgData[uid]['prop_view'] = HTML(description='', value='initializing...')
        wgData[uid]['prop_wg'] = Dropdown(description='',options=['select'], value='select')
        wgBox.append(Box([Label(value='prop'),wgData[uid]['prop_view']], layout=itemLayout))
        wgBox.append(Box([Label(value='prop'),wgData[uid]['prop_wg']], layout=itemLayout))
        widgetsFor3DView.append('prop_wg')
        
    # useDrawAs
    if useDrawAs:
        if moldict is not None:
            wgData[uid]['drawAs_wg'] = Dropdown(description='', options=DRAWING_LIGAND_3D, value=drawAs)
            wgBox.append(Box([Label(value='drawAs'),wgData[uid]['drawAs_wg']], layout=itemLayout))
            widgetsFor3DView.append('drawAs_wg')
        if protein is not None:
            wgData[uid]['pStyle_wg'] = Dropdown(description='', options=DRAWING_PROTEIN_3D, value=pStyle)
            wgBox.append(Box([Label(value='pStyle'),wgData[uid]['pStyle_wg']], layout=itemLayout))
            wgData[uid]['pStyle_helicesAsTubes_wg'] = Checkbox(description='helicesAsTubes', value=False)
            wgBox.append(Box([Label(value=''),wgData[uid]['pStyle_helicesAsTubes_wg']], layout=itemLayout))
            widgetsFor3DView.extend(['pStyle_wg', 'pStyle_helicesAsTubes_wg'])
            
    # colors
    if colorPanel is True:
        wgData[uid]['colorScheme'] = Dropdown(description='', options=LIGAND_COLOR_SCHEME_3D, value='default')
        wgBox.append(Box([Label(value='ligand color'),wgData[uid]['colorScheme']], layout=itemLayout))
        widgetsFor3DView.append('colorScheme')
        
    # labels
    if labelPanel is True:
        wgData[uid]['confLabel'] = Checkbox(description='confLabel', value=False)
        wgData[uid]['atomLabel'] = Checkbox(description='atomLabel', value=False)
        cbLabel=Box([Label(value=''), HBox([wgData[uid]['confLabel'], 
                                            wgData[uid]['atomLabel']
                                           ])
                    ], layout=itemLayout)
        wgBox.append(cbLabel)
        widgetsFor3DView.extend(['confLabel', 'atomLabel'])
        
    # background
    wgData[uid]['background'] = Dropdown(description='', options= BGCOLORS_3D, value=BGCOLORS_3D[0])
    wgBox.append(Box([Label(value='background'),wgData[uid]['background']], layout=itemLayout))
    widgetsFor3DView.append('background')
    
    
    # buttons
    wgData[uid]['start'] = Button(description="Start!", button_style='success')
    wgData[uid]['zoomTo'] = Button(description="zoomTo", button_style='success')
    
    buttons=Box([Label(value=''),HBox([wgData[uid]['zoomTo'], wgData[uid]['start']])],layout=itemLayout)
    wgBox.append(buttons)
    
    
    # Observe and update dict with model_id
    wgData[uid]['start'].on_click(handle_start_button)
    wgData[uid]['zoomTo'].on_click(handle_zoomTo_button)
    wgMIds[wgData[uid]['start'].model_id] = uid
    wgMIds[wgData[uid]['zoomTo'].model_id] = uid
    
    wgDataUid = wgData[uid]
    for i in widgetsFor3DView:
        selectedWidget = wgDataUid[i]
        selectedWidget.observe(handle_change, names='value')
        wgMIds[selectedWidget.model_id] = uid
    
    
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
    