def buildDistMat(sz):
  mask = 0x1
  ds = []
  for i in range(sz):
    for j in range(sz):
      dist=0
      a=i
      b=j
      while a or b:
        ta = a&mask
        tb = b&mask
        if ta^tb:
          dist+=1
        a >>= 1
        b >>= 1
      ds.append(dist)
  return ds
def numOnBits(sz):
  mask = 0x1
  ds = []
  for i in range(sz):
    nOn = 0
    a = i
    while a:
      if a&mask:
        nOn+=1
      a >>= 1
    ds.append(nOn)
  return ds
if __name__=='__main__':
  sz = 256
  #ds = buildDistMat(sz)
  ds = numOnBits(sz)

  step = 32
  i=0
  while i<len(ds):
    print '  %s,'%(','.join([str(x) for x in ds[i:i+step]]))
    i+=step
  

        
