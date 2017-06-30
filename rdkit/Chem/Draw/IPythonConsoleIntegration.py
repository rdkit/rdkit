
"""RDKit Conformer Browser
Derived from Greg Landrum's code by Malitha Humayun Kabir
As a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Date: 28th June 2017
"""

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from ipywidgets import interact, interactive, fixed
import ipywidgets as widgets

molSize_3d = (400, 400)
bgcolor_3d = '0xeeeeee'

def startViewer(size=None,bgcolor=None):
    if size is None:
        size=molSize_3d
    if bgcolor is None:
        bgcolor=bgcolor_3d
    view = py3Dmol.view(width=size[0],height=size[1])
    view.setBackgroundColor(bgcolor)
    return view


def processSuppliedMolFor3DViewer(ms):
    
    try:
        # list of tuple (name + mol obj) # dict key possible
        # ms = [('m1', <rdkit.Chem.rdchem.Mol at 0x7246d88>),
        # ('m2', <rdkit.Chem.rdchem.Mol at 0x7246ea0>),
        # ('m3', <rdkit.Chem.rdchem.Mol at 0x7246ed8>),
        # ('m4', <rdkit.Chem.rdchem.Mol at 0x7246f10>)]
        moldict = dict(ms)
        
    except TypeError:
        
        if type(ms) is tuple:
            # Not list single tuple (name + mol obj) # dict key possible
            # ms = ('m1', <rdkit.Chem.rdchem.Mol at 0x7246d88>)
            moldict=list()
            moldict.append(ms)
            moldict = dict(moldict)
        elif hasattr(ms, '__iter__') is False:
            # rdkit.Chem.rdchem.Mol
            # Not list... mol obj only ... no name... no dict key possible from such list...
            # So, generate dict key
            # ms = <rdkit.Chem.rdchem.Mol object at 0x07246D88>
            moldict=list()
            moldict.append(('m0', ms))
            moldict = dict(moldict)
        elif type(ms) is list:
            # list of mol obj only ... no name... no dict key possible from such list...
            # So, generate dict key
            # ms = [<rdkit.Chem.rdchem.Mol at 0x7246d88>,
            # <rdkit.Chem.rdchem.Mol at 0x7246ea0>,
            # <rdkit.Chem.rdchem.Mol at 0x7246ed8>,
            # <rdkit.Chem.rdchem.Mol at 0x7246f10>]
            ms_name=['m'+str(x) for x in range(len(ms))]
            ms2=[(ms_name[i],ms[i]) for i in range(len(ms))]
            moldict = dict(ms2)
    return moldict
    
def addMolToViewForScrolling(# fixed objects
                             molecule_obj_list, 
                             view,
                             widget_uniqueid,
                             # property viewer
                             property_view_existing,
                             property_view_rdkit_offered,
                             # molecule and conformer selector
                             molecule_select,
                             confId,
                             # value from property select dropdown
                             existing_descriptors,
                             rdkit_offered_descriptors,
                             # 3DMol.js viewer related stuffs Dropdowm
                             useDrawAs,
                             drawAs,
                             molColorScheme,
                             # value from labeling checkbox True/False
                             showConfLabel,
                             showAtomLabel
                            ):
    
    #
    # Two widgets got NO effect those are just to view 
    # property_view_existing, property_view_rdkit_offered
    # 
    # existing_descriptors is for getting the selected value from dropdown
    # existing_descriptors_dropdown is handled as global object
    # 
    
    
    #### Molecule and conformer selection chunk
    
    # Get mol from supplied list object
    mol = molecule_obj_list[molecule_select]
    
    if mol.GetNumConformers()>0:
        # Update conformer select slider
        globals()['confId_slider_'+widget_uniqueid].max=mol.GetNumConformers()-1
        # For conformers
        conf_selected=mol.GetConformer(confId)
        xyz=conf_selected.GetPositions()
        # For owning Mol
        OwningMol=conf_selected.GetOwningMol()
        
    #### To Do: If the Mol object does not have any conformers????
    
    
    #### removing descriptor widget description
    globals()['existing_descriptors_viewer_'+widget_uniqueid].description=''
    globals()['rdkit_offered_descriptors_viewer_'+widget_uniqueid].description=''
    
    #### Descriptor handling chunkGetting descriptors
    
    # Get available properties and update dropdown widget 
    available_properties=list(OwningMol.GetPropNames())
    
    if len(available_properties)>0:
        # update dropdown with existing properties
        globals()['existing_descriptors_dropdown_'+widget_uniqueid].options=available_properties
        
        # Get selected property name and associated value
        des_name_existing=eval('existing_descriptors_dropdown_'+widget_uniqueid+'.value')
        des_val_existing=OwningMol.GetProp(des_name_existing)
        
        # Update viewer for existing descriptor with selected property name and value
        globals()['existing_descriptors_viewer_'+widget_uniqueid].value='Existing descriptor: '\
                                            + des_name_existing \
                                            + ': ' \
                                            + str(des_val_existing)
        #
        #
    #
    #
    # Update viewer for rdkit offered descriptor with selected property name and value
    # This is real time calculation
    #
    # Descriptor calculation schema eval("Descriptors.TPSA(OwningMol)")
    # In above line descriptor_list is TPSA
    #
    des_cmd="Descriptors."+rdkit_offered_descriptors+"(OwningMol)"
    des_val=str(eval(des_cmd))
    globals()['rdkit_offered_descriptors_viewer_'+widget_uniqueid].value=\
                                                        'RDKit offered descriptor: '\
                                                        +rdkit_offered_descriptors\
                                                        +': '\
                                                        +des_val
    
    
    ##### 3Dmol.js viewer realated stuffs below
    
    #
    # Clearing previous 3Dmol objects withoiut resetting view
    #
    view.removeAllModels()
    view.removeAllSurfaces()
    view.removeAllLabels()
    
    #
    #### Adding model to viewer
    #
    if mol.GetNumAtoms()>=999 or drawAs == 'cartoon':
        # py3DMol is happier with TER and MASTER records present
        pdb = Chem.MolToPDBBlock(mol,flavor=0x20|0x10)
        view.addModel(pdb,'pdb')
    else:
        # py3Dmol does not currently support v3k mol files, so
        # we can only provide those with "smaller" molecules
        mb = Chem.MolToMolBlock(mol,confId=confId)
        view.addModel(mb,'sdf')
        
    #
    #### Making decision about model rendering style and color
    #
    
    if useDrawAs is False:
        #use from globalStyle
        view.setStyle({},{globals()['myGlobalStyle_'+widget_uniqueid]:{'colorscheme': molColorScheme}})
    else:
        #update global style and use that
        globals()['myGlobalStyle_'+widget_uniqueid] = drawAs
        view.setStyle({},{globals()['myGlobalStyle_'+widget_uniqueid]:{'colorscheme': molColorScheme}})
    
    #
    # This is exception for surface
    #
    if drawAs is 'surface':
        view.addSurface({}, '$3Dmol.SurfaceType.VDW');
        
    if drawAs is 'ballstick':
        view.setStyle({},{'stick':{'radius':'0.2','colorscheme': molColorScheme},
                          'sphere':{'radius':'0.4', 'colorscheme': molColorScheme}}
                     );
    
    #
    # Labeling conformer
    #
    if showConfLabel is True:
        label=molecule_select+':'+str(confId)
        view.addLabel(label, {'backgroundColor':'gray', 'fontColor':'white',
                              'showBackground':'true', 'alignment':'bottomCenter'})
    
    #
    # Labeling atom
    #
    if showAtomLabel is True:
        
        label_create=[OwningMol.GetAtomWithIdx(i).GetSymbol()+
                      str(OwningMol.GetAtomWithIdx(i).GetIdx()+1) 
                      for i in range(conf_selected.GetNumAtoms())
                     ]
        
        i = None
        for i in range(conf_selected.GetNumAtoms()):
            view.addLabel(label_create[i], {'inFront' : 'false', 
                                            'fontSize' : '12',
                                            'fontColor':'gray',
                                            'showBackground':'false',
                                            'position' : {'x' : xyz[i][0],
                                                          'y' : xyz[i][1],
                                                          'z' : xyz[i][2]
                                                       }
                                           })
    
    #print(drawAs)
    # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
    #view.zoomTo()
    return view.update()


def browseMolConformers(mol_obj_list,
                        view, 
                        confId=None, 
                        useDrawAs=False, 
                        drawAs=None):
    
    
    # Listing all widgets for tracking purpose
    # 
    # molecule_list_dropdown
    # confId_dropdown
    # 
    # existing_descriptors_viewer
    # existing_descriptors_dropdown
    # rdkit_offered_descriptors_viewer
    # rdkit_offered_descriptors_dropdown
    # 
    # useDrawAs_dropdown
    # drawAs_dropdown
    # molColorScheme_dropdown
    # 
    # confLabelCheckBok
    # atomLabelCheckBox
    # 
    
    
    # processing supplied object that contains molecules
    moldict = processSuppliedMolFor3DViewer(mol_obj_list)
    
    # Creating molecule listing widget # default value molecule id = 0 or m0
    molecule_list=list(moldict.keys())
    molecule_list_dropdown = widgets.Dropdown(options=molecule_list,value=molecule_list[0])
    
    # Creating conformer listing widget # default value confId=0
    globals()['confId_slider_'+view.uniqueid] = widgets.IntSlider(min=0,max=9,step=1,value=0)
    #IntSlider(min=-10,max=30,step=1,value=10)
    
    
    #### A set of widgets for descriptor handling
    
    globals()['existing_descriptors_viewer_'+view.uniqueid]=widgets.HTML(value='Trying to get existing property', disabled=False)
    
    globals()['existing_descriptors_dropdown_'+view.uniqueid]=widgets.Dropdown(options=['please select'],value='please select')
    
    globals()['rdkit_offered_descriptors_viewer_'+view.uniqueid]=widgets.HTML(value='Trying to get rdkit offered property', disabled=False)
    
    
    
    
    descriptors_rdkit=['MolLogP', 'MolMR', 'MolWt', 'ExactMolWt', 'HeavyAtomCount', 
                       'HeavyAtomMolWt', 'NHOHCount', 'NOCount', 'NumHAcceptors', 
                       'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumValenceElectrons']
    rdkit_offered_descriptors_dropdown = widgets.Dropdown(options=descriptors_rdkit,value='MolLogP')
    
    
    # default myGlobalStyle = 'stick'
    # possible other options line cross stick cartoon sphere
    globals()['myGlobalStyle_'+view.uniqueid] = 'stick'
    
    if useDrawAs is False:
        # Then go with default settings
        drawAs = globals()['myGlobalStyle_'+view.uniqueid]
        useDrawAs_dropdown = widgets.Dropdown(options=[False, True],value=False)
        
    else:
        useDrawAs_dropdown = widgets.Dropdown(options=[False, True],value=True)
        # Use user supplied drawAS
        if drawAs is None:
            # User forgot to provide drawAs argument
            drawAs = globals()['myGlobalStyle_'+view.uniqueid]
        else:
            # User supplied drawAs argument while useDrawAs is True
            globals()['myGlobalStyle_'+view.uniqueid] = drawAs
        
    
    all_drawing_types_3d=['line', 'cross', 'stick', 'cartoon', 'sphere', 'surface', 'ballstick']
    # This is widget for model style
    drawAs_dropdown = widgets.Dropdown(options=all_drawing_types_3d,value=drawAs)
    
    
    all_color_scheme_3d=['default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon', 
                     'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                     'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                     'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop']
    
    # This is widget for model style
    
    molColorScheme_dropdown = widgets.Dropdown(options=all_color_scheme_3d,value='default')
    
    
    
    # This is widget for conformer label
    
    confLabelCheckBok = widgets.Checkbox(description='confLabelCheckBok', value=False)
    
    # This is widget for atom label of each conformers
    
    atomLabelCheckBox = widgets.Checkbox(description='atomLabelCheckBox', value=False)
    
    
    # Now start interacting
    result=interact(addMolToViewForScrolling,
                    
                    ####
                    molecule_obj_list=fixed(moldict),
                    view=fixed(view),
                    widget_uniqueid=fixed(view.uniqueid),
                    ####
                    
                    
                    property_view_existing=eval('existing_descriptors_viewer_'+view.uniqueid),
                    property_view_rdkit_offered=eval('rdkit_offered_descriptors_viewer_'+view.uniqueid),
                    
                    molecule_select=molecule_list_dropdown,
                    confId=eval('confId_slider_'+view.uniqueid),
                    
                    existing_descriptors=eval('existing_descriptors_dropdown_'+view.uniqueid),
                    rdkit_offered_descriptors=rdkit_offered_descriptors_dropdown,
                    
                    useDrawAs=useDrawAs_dropdown,
                    drawAs=drawAs_dropdown,
                    molColorScheme=molColorScheme_dropdown,
                    
                    showConfLabel=confLabelCheckBok,
                    showAtomLabel=atomLabelCheckBox
                    ####
                   );
    return result

