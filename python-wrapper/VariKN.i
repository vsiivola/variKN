%include "exception.i"
%include "std_string.i"
%include "std_vector.i"
%module varikn

%exception {
  try {
    $action
  }
  catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    SWIG_exception(SWIG_RuntimeError, "Exception");
  }
  catch (std::string &e) {
    std::cerr << "Exception: ";
    std::cerr << e << std::endl;
    SWIG_exception(SWIG_RuntimeError, "Exception");
  }
	catch (...) {
    SWIG_exception(SWIG_RuntimeError, "Unknown Exception");
	}
}

#if defined(SWIGPYTHON)
%typemap(in) std::string& {
  if (!PyString_Check($input)) {
    PyErr_SetString(PyExc_TypeError, "not a string");
    return NULL;
  }
  $1 = new std::string(PyString_AsString($input),
           PyString_Size($input));
}

%typemap(freearg) std::string& {
  delete $1;
}

%typemap(out) std::string& {
  $result = Py_BuildValue("s#",$1->c_str(),$1->size());
}

#ifndef PYTHON3 
## Python 3 no longer supports FILE* typemaps, this can be enabled for py2
%typemap(in) FILE* {
        if (!(PyFile_Check($input))) {
                PyErr_SetString(PyExc_TypeError, "not a file pointer");
                return NULL;
        }
        $1=PyFile_AsFile($input);}
#endif
#endif

%inline %{
  typedef int indextype;
%}

%template(stringvector) std::vector<std::string>;
%template(floatvector) std::vector<float>;

%feature("notabstract") InterTreeGram;
%include Vocabulary.hh
%include NGram.hh
%include InterTreeGram.hh
%include PerplexityFuncs.hh
%include VarigramFuncs.hh
%{
#include "Vocabulary.hh"
#include "NGram.hh"
#include "InterTreeGram.hh"
#include "PerplexityFuncs.hh"
#include "VarigramFuncs.hh"
%}

%template(VarigramTrainer) Varigram_t<int, int>;

