Documentation for p_con (https://github.com/pzc/rdkit/blob/master/Contrib/pzc/p_con.html)

from p_con import p_con


pco = p_con("P43088")
pco.verbous = True

pco.step_0_get_chembl_data() # Download Compounds for P43088 from ChEMBL
(or pco.load_mols("sdf-file.sdf.gz"))

pco.step_1_keeplargestfrag() # remove small Fragments from compounds

pco.step_2_remove_dupl()     # remove duplicate-Entries

pco.step_3_merge_IC50()      # merge IC50 from Entries with same canonical smiles into one compound

pco.step_4_set_TL(4000,ic50_tag="value") # set TrafficLights, value > 4000nm: 0, else 1

pco.step_5_remove_descriptors() # remove Descriptors from compounds

pco.step_6_calc_descriptors() # calculate new Descriptors which are used to create prediction-models

pco.step_7_train_models() # train up to 10 models

pco.save_model_info("model_info.csv",mode="csv")   # create csv with data for each model
pco.save_model_info("model_info.html",mode="html") # create html -#-

for i in range(len(pco.model)):
    pco.save_model("model_%d.pkl" % i,i)

for i in range(len(pco.model)):
    act,inact = pco.predict(i)
    print "Model %d active: %d\tinactive: %d" % (i,act,inact)



# to Check compounds using Models

pco2 = p_con("P43088")
pco2.verbous = True
pco2.load_mols("P43088.sdf.gz")
models = ["model1.pkl","model2.pkl"]
pco2.load_models(models)

print "\n#Model\tActive\tInactive"

for i in range(len(self.model)):
    act,inact = pco2.predict(i)
    print "%d\t%d\t%d" % (i,act,inact)
