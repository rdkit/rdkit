Glare Algorithm.

Implementation of

    GLARE: A New Approach for Filtering Large Reagent Lists in 
           Combinatorial Library Design Using Product Properties
    Jean-Francois Truchon* and Christopher I. Bayly

    http://pubs.acs.org/doi/pdf/10.1021/ci0504871

    Usage:
       # somehow make sidechains1/2 with props [mw, alogp, tpsa]
       r1 = RGroups(sidechains1)
       r2 = RGroups(sidechains2)
       lib = Library([r1, r2])
       props = [
         Property("mw", 0, 0, 500),
         Property("alogp", 1, -2.4, 5),
         Property("tpsa", 2, 0, 90)
       ] 

      glare = Glare()
      glare.optimize(lib, props)
      
Notes:
Some nomenclature:

 A Libary is made of RGroups
 RGroups are a collection of sidechains (the paper uses Fragments)
  that can populate the rgroup position.

 We desire to optimize the Library so that we have a good chance
  of making the desired products.

 From the testing code, using Fake data:

    r1 = RGroups(makeFakeSidechains("aldehydes", num=1000))
    r2 = RGroups(makeFakeSidechains("boronic_acids", num=1500))
    
    libs = Library([r1,r2])
    props = [
        Property("mw",    propIdx=0, minValue=0, maxValue=500),
        Property("alogp", propIdx=1, minValue=-2.4, maxValue=5),
        Property("tpsa",  propIdx=2, minValue=0, maxValue=90)
    ]
    
    glare = Glare()
    # optimize the library...
    glare.optimize(libs, props)
