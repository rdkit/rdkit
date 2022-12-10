/* -----------------------------------------------------------------------------
 * extend_std_vector.i
 * ----------------------------------------------------------------------------- */


%extend std::vector {
  bool equals(const vector<T> &o){
    if(self->size()==o.size()){
      std::vector< T >::const_iterator sIt=self->begin();
      std::vector< T >::const_iterator oIt=o.begin();
      while(sIt != self->end()){
        if(*sIt != *oIt) return false;
        ++sIt;
        ++oIt;
      }
      return true;
    } else {
      return false;
    }
  }
 };

// suggested in the CHANGELOG for SWIG 4.0
%extend std::vector {
    vector(size_type count) { return new std::vector< T >(count); }
}

%include "std_vector.i"




