import python_objs

def lentest():
  assert python_objs.seq_len([1,2])==2
  assert python_objs.seq_len((1,2))==2
  assert python_objs.seq_len("fooo")==4
  try:
    python_objs.seq_len(1)
  except TypeError:
    pass
  except:
    assert 0,'wrong exception type'
  else:
    assert 0,'no exception'


def sumtest():
  assert python_objs.sum_first2([1,2])==3
  assert python_objs.sum_first2((1,3))==4
  try:
    python_objs.sum_first2((1.,3.))
  except TypeError:
    pass
  except:
    assert 0,'wrong exception type'
  else:
    assert 0,'no exception'
  try:
    python_objs.sum_first2('foo')
  except TypeError:
    pass
  except:
    assert 0,'wrong exception type'
  else:
    assert 0,'no exception'
  try:
    python_objs.sum_first2(1)
  except TypeError:
    pass
  except:
    assert 0,'wrong exception type'
  else:
    assert 0,'no exception'


lentest()
sumtest()
