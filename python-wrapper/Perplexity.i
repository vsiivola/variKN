# -*- tab-width: 2 -*-

%include exception.i

%module Perplexity

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
 
%typemap(in) FILE* {
        if (!(PyFile_Check($input))) {
                PyErr_SetString(PyExc_TypeError, "not a file pointer");
                return NULL;
        }
        $1=PyFile_AsFile($input);}
#endif

%include PerplexityFuncs.hh
%include VarigramFuncs.hh
%{
#include "PerplexityFuncs.hh"
#include "VarigramFuncs.hh"
%}
