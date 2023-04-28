"""
Glare Algorithm

Some nomenclature:

 GLARE: A New Approach for Filtering Large Reagent Lists in 
           Combinatorial Library Design Using Product Properties
    Jean-Francois Truchon* and Christopher I. Bayly

    http://pubs.acs.org/doi/pdf/10.1021/ci0504871

 A Libary is made of RGroups
 RGroups are a collection of sidechains (the paper uses Fragments)
  that can populate the rgroup position.

 We desire to optimize the Library so that we have a good chance
  of making the products in the desired property range.

 Example From the testing code, using Fake data:

    r1 = RGroups(makeFakeSidechains("aldehydes", num=1000))
    r2 = RGroups(makeFakeSidechains("boronic_acids", num=1500))
    
    lib = Library([r1,r2])
    props = [
        Property("mw",    propIdx=0, minValue=0, maxValue=500),
        Property("alogp", propIdx=1, minValue=-2.4, maxValue=5),
        Property("tpsa",  propIdx=2, minValue=0, maxValue=90)
    ]
    
    glare = Glare()
    # optimize the library...
    glare.optimize(lib, props)
    for reactant_idx, rgroup in enumerate(lib.rgroups):
        print(f"Reactants for reactant {reactant_idx}")
        for reactant in rgroup.sidechhains:
            print(reactant.name)
"""
import itertools
import math
import operator
import random
from functools import reduce


class Property:

  def __init__(self, name, propIdx, minValue, maxValue, scaffoldoffset=0.0):
    """name, propIdx, minValue, maxValue, scaffoldoffset -> initial a Property
        name is the name of the property.
         propIdx:  the index of the property in the property vector
         minValue: the minimum acceptable value for the property
         maxValue: the maximum acceptable value for the property
         scaffoldoffset: any offset from the reaction scaffold (defaults to 0)
        """
    self.name = name
    self.propIdx = propIdx
    self.minValue = minValue
    self.maxValue = maxValue
    self.offset = scaffoldoffset

  def evaluate(self, sidechains):
    """sidechains -> Evaluate a list of sidechains to see if they
        pass the property values.

        Each sidechain must have a property vector e.g. (s.props for s in sidechains)
        which is a vector of values where s.props[propIdx] is the property
        being inspected
        """
    product = self.offset
    propIdx = self.propIdx
    for s in sidechains:
      product += s.props[propIdx]
    return self.minValue <= product <= self.maxValue


class Sidechain:
  """Holds the name (identifier) and property list for the
    given sidechain/fragment.  Properties are assumed to
    be numerical values"""

  def __init__(self, name, props, goodCount=0, **extra_data):
    """name, props, goodCount=0 -> initialize a Sidechain
        initialize a sidechain.
        name: the unique name for the sidechain
        props: the property vector (see Properties class for details)
        goodCount: the number of times this reagent belongs to
            a good product, where good is a product that is in the desired
            property space.
        """
    self.name = name
    self.props = props
    self.good_count = goodCount  # shared variable
    self.dropped = False  # shared variable
    self.extra_data = extra_data

  def data(self):
    return self.extra_data

  def __str__(self):
    return "Sidechain %s(%s, goodCount=%s, **%r)" % (self.name, self.props, self.good_count,
                                                     self.extra_data)

  def __repr__(self):
    return "Sidechain(%r, %r, %s, **%r)" % (self.name, self.props, self.good_count, self.extra_data)


class RGroups:
  """Holds a collection of sidechains for the given RGroup"""

  def __init__(self, sidechains):
    """Sidechains -> RGroups
         sidechains: the list of Sidechains that make up the potential
                     sidechains at this rgroup position"""
    self.sidechains = sidechains

    self.rejected = []  # list of rejected sidechains
    self.initial_size = len(sidechains)

  def count(self):
    """Returns the number of possible sidechains"""
    return len(self.sidechains)

  def randomize(self):
    """Randomly shuffles the sidechains and reset the goodness counts"""
    random.shuffle(self.sidechains)
    for s in self.sidechains:
      s.good_count = 0

  def effectiveness(self):
    """-> return the current effectiveness of this collection
        effectiveness is the number of items left divided by the 
        initial amount"""

    return len(self.sidechains) / float(self.initial_size)

  def chunk_size(self, num_chunks):
    """num_chunks -> return the number of sidechains in each chunk
        if the sidechains are split into num_chunks chunks"""
    return int(math.ceil(float(len(self.sidechains)) / num_chunks))

  def chunk(self, chunk_idx, num_chunks):
    """chunk_idx, num)chunks -> RGroups
        return the chunk_idxth chunk given num_chunks total chunks"""
    assert chunk_idx >= 0 and chunk_idx < num_chunks, "%s %s" % (chunk_idx, num_chunks)

    n = self.chunk_size(num_chunks)
    return RGroups(self.sidechains[chunk_idx * n:(chunk_idx + 1) * n])

  def prune(self, fractionToKeep):
    """fractionToKeep -> Sort the sidechains from the most often 
        found if good products to the least, and keep the best 
        fractionToKeep percentage"""
    assert 0 < fractionToKeep <= 1.0, "fractionToKeep: %s" % fractionToKeep

    self.sidechains.sort(key=lambda x: x.good_count, reverse=True)
    fragment_index = int(len(self.sidechains) * fractionToKeep + 0.5)

    # update rejected set
    self.rejected += self.sidechains[fragment_index:]
    self.sidechains = self.sidechains[:fragment_index]


class Library:
  """A library is a collection of RGroups that need to be combinitorially
    combined"""

  def __init__(self, rgroups):
    """rgroups -> Initialize the Library.
        rgroups: the list of possible RGroups that is combinitorially
                 combined to make the library"""
    self.rgroups = rgroups

  def isValid(self):
    """If we have an empty set for any rgroup, return False"""
    for rg in self.rgroups:
      if len(rg.sidechains) == 0:
        return False
    return True

  def randomize(self):
    """randomize the order of the sidechains"""
    for rg in self.rgroups:
      rg.randomize()

  def getSidechainsPerPartition(self, total_num_partitions_per_rgroup):
    """total_num_partitions -> [num_fragments/partition for rgroup1, 
                                    num_fragments/partition for rgroup2]
        return the number of sidechains in a partition
        for each rgroup"""

    sizes = [(libIdx, max(rg.count() // total_num_partitions_per_rgroup, 1))
             for libIdx, rg in enumerate(self.rgroups)]

    # "optimally" apportion the partitions according the
    #  the glare paper see Appendix eq (8) and (9)
    # sort by size
    sizes.sort(key=lambda sz: sz[1])
    last_size = 1
    opt_sizes = []
    for libIdx, current_size in sizes[:-1]:
      opt_sizes.append((libIdx, current_size - (current_size % last_size)))
      last_size = current_size

    # From the Glare paper:
    # the last library size is set equal to the second to last
    # From Table 3, it is easy to understand that, if the fourth dimension
    # was split in 24 instead of 12, a factor of 2 would be gained from the
    # reduced size of the sublibraries. However, twice as many sublibraries
    # would be needed, and the net speedup would be null, hence, the decision to
    # set p4=p3. (p4 here is the last library)
    libIdx, current_size = sizes[-1]
    opt_sizes.append((libIdx, last_size))
    # back to the original library order
    opt_sizes.sort()
    res = [size for libIdx, size in opt_sizes]
    return res

  def chunk(self, num_partitions):
    """num_partitions -> [Library(..), Library(...)]

        Return new libraries that are chunks of this one.
        These are the libraries that get sampled to see of
        sidechains participate in good products.
        """
    partitions = self.getSidechainsPerPartition(num_partitions)
    max_subsets = max(partitions)
    enumeration_indices = []
    for i in range(max_subsets):
      combinations = []
      for size in partitions:
        combinations.append(i % size)
      enumeration_indices.append(combinations)

    library_sets = []
    for subset_index, combinations in enumerate(enumeration_indices):
      libs = []
      partitioned_rgroups = []
      for lib_index, libpart_index in enumerate(combinations):
        lib = self.rgroups[lib_index]
        num_chunks = partitions[lib_index]
        partitioned_rgroups.append(lib.chunk(chunk_idx=libpart_index, num_chunks=num_chunks))
      lib = Library(partitioned_rgroups)
      if lib.isValid():
        library_sets.append(lib)

    return library_sets

  def effectiveness(self):
    """-> returns the average effectiveness of this library set"""
    sum = 0.0
    for rg in self.rgroups:
      sum += rg.effectiveness()
    return sum / len(self.rgroups)

  def evaluate(self, props):
    """props -> num_good_enumerations, total_enumerations

        props: a list of Property evaluators for the fragments.

        returns the number of good enumerations and the total number of
        enumerations for this Library"""
    frags = [rg.sidechains for rg in self.rgroups]
    good = 0
    bad = 0
    for i, frag in enumerate(itertools.product(*frags)):
      for p in props:
        if not p.evaluate(frag):
          bad += 1
          break
      else:
        good += 1
        for sidechain in frag:
          sidechain.good_count += 1
    return good, i + 1


class Glare:
  """Glare Algorithm.  Implementation of

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
    """

  def __init__(
      self,
      desiredFinalGoodness=0.95,
      maxIterations=100,
      rgroupScale=6.0,  # None if no scaling
      initialFraction=None,  #None=auto -100.,
      numPartitions=16):
    self.fractionGood = self.desiredFinalGoodness = desiredFinalGoodness
    self.maxIterations = maxIterations
    self.rgroupScale = rgroupScale

    if initialFraction is not None:
      self.initialFraction = initialFraction / 100.
    else:
      self.initialFraction = initialFraction
    self.numPartitions = numPartitions

  def optimize(self, library, props):
    """library, props
        Given a Library and the list of Propery evaluators,
        optimize the library.
        The library is modified in place by removing building blocks
        (sidechains) that are not likely to pass the property
        criteria.
        """
    # attempt to generate report like glare application
    print("------- PARAMETERS: --------------")
    print("GOOODNESS THRESHOLD : %s%%" % (self.desiredFinalGoodness * 100))
    print("MIN PARTITION SIZE : %s" % self.numPartitions)
    if self.initialFraction is None or self.initialFraction > 0.999:
      print("INITIAL FRACTION TO KEEP : AUTOMATIC")
    else:
      print("INITIAL FRACTION TO KEEP : %s%%" % (self.initialFraction * 100))
    print("Actual SIZE : %s = %s" %
          (" x ".join([str(len(rg.sidechains)) for rg in library.rgroups
                       ]), reduce(operator.mul, [len(rg.sidechains) for rg in library.rgroups])))

    running_total = 0.0
    Gt = self.desiredFinalGoodness

    for iteration in range(1, self.maxIterations + 1):
      # chunk of the total library into smaller more manageable sets
      #  and run combinatorial analysis on the sub libraries
      #  each of these records the number of times a sidechain is used
      #  in a successful enumeration which is then used to prune the
      #  library at the end
      #
      for rg in library.rgroups:
        rg.randomize()

      good = total = 0.0
      chunked_libs = library.chunk(self.numPartitions)
      # for each chunk, do the combinatorial check to see
      #  if reagents make good products
      for libidx, chunk in enumerate(chunked_libs):
        g, t = chunk.evaluate(props)
        good += g
        total += t
      running_total += total
      Gi = good / total  # current goodness

      if Gi < 1e-12:
        # I think we're done here :)
        fraction = 0.0
      elif iteration == 1:
        G0 = Gi  # Goodness at first iteration

        # the first time, use the initalFraction or a "good enough"
        #  value
        if self.initialFraction is not None:
          fraction = K0 = self.initialFraction
        else:
          # auto choose the fraction based on the current good percentage
          #  and the desired
          fraction = K0 = min(-1.1 * (Gt - G0) + 1.2, 0.9)
      else:
        # the second time, gradually eliminate reagents slowing
        #  down as the number of iterations increases
        #  see equation (5) in reference
        if abs(Gt - G0) < 1e-4:
          Ki = 1.0
        else:
          Ki = (1.0 - K0) * (Gi - G0) / (Gt - G0) + K0
        fraction = min(1.0, Ki)

      # prune the library to keep the highest occurring sidechains
      #  note that even if all sidechains are acceptable,
      #  some will always get pruned

      max_size = float(max([len(rg.sidechains) for rg in library.rgroups]))
      for rg in library.rgroups:
        scale = 1.0
        if self.rgroupScale is not None:
          # scale differently size rgroups via equation (6)  in paper
          numSidechains = len(rg.sidechains)
          numer = 1.0
          denom = 1.0 + math.exp(-self.rgroupScale * ((numSidechains / max_size) - 0.5))
          scale = numer / denom
        fraction_to_reject = (1.0 - fraction) * scale
        # keep the best fraction...
        rg.prune(1.0 - fraction_to_reject)

      print("-------------- ITERATION : %s ----------------------" % iteration)
      print("GOODNESS      : %s%%" % (Gi * 100))
      print("NUMBER EVAL   : %s" % (total))
      print("CUMUL EVAL    : %s" % (running_total))
      print("KEPT IN STEP  : %s%%" % (fraction * 100.))
      if not iteration:
        print("GOODNESS THRESHOLD : %s" % self.desiredFinalGoodness)
        print("MIN PARTITION SIZE : %s" % self.numPartitions)
        print("INITIAL FRACTION TO KEEP : ")
        if self.fractionToKeep > 0.999:
          print("AUTOMATIC")
        else:
          print("%s%%" % self.fractionToKeep)

      print("Actual SIZE : %s = %s" %
            (" x ".join([str(len(rg.sidechains)) for rg in library.rgroups
                         ]), reduce(operator.mul, [len(rg.sidechains) for rg in library.rgroups])))
      print("EFFECTIVENESS : %s%%" % (library.effectiveness() * 100.))

      # stopping criteria
      if iteration and Gi < 1e-12:
        return
      elif abs(Gi - self.desiredFinalGoodness) < 0.001 or \
           Gi > self.desiredFinalGoodness:
        return


######################################################################
# testing codes
def makeFakeProps():
  mw = random.randint(10, 500)
  alogp = random.randint(-10, 10)
  tpsa = random.randint(0, 180)
  return [mw, alogp, tpsa]


def makeFakeSidechains(lib, num):
  res = []
  for i in range(num):
    res.append(Sidechain(lib + "_" + str(i), makeFakeProps()))
  return res


def testGlare():
  a = RGroups(makeFakeSidechains("aldehydes", 1000))
  b = RGroups(makeFakeSidechains("boronic_acids", 1500))

  lib = Library([a, b])
  props = [
    Property("mw", 0, 0, 500, 230.1419),
    Property("alogp", 1, -2.4, 5, 2.212749),
    Property("tpsa", 2, 0, 90, 24.5)
  ]

  glare = Glare()
  glare.optimize(lib, props)
  # print out the selected reactants
  for reactant_idx, rgroup in enumerate(lib.rgroups):
    print(f"Reactants for reactant {reactant_idx}")
    for reactant in rgroup.sidechains:
      print(reactant.name)


if __name__ == "__main__":
  testGlare()
