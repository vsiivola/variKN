###################################
# Set OS and compiler specific flags
ARCH = $(shell uname -m)

# Opteron cluster
ifeq ($(ARCH),x86_64)
OPT ?=-O2 -fPIC
endif

ifeq ($(ARCH),i686)
OPT ?=-O2 -fPIC
endif

ifeq ($(ARCH),i386)
OPT ?=-O2 -fPIC
endif

# IRIX
ifeq ($(ARCH),IP27)
CXX = CC
ifeq ($(CXX),CC)
OPT ?= -Ofast -r10000
CFLAGS = $(OPT) -Dno_inline_funcs 
CXXFLAGS = -LANG:std -DDMEDIA -woff 1183
endif
endif

# Tru64
ifeq ($(ARCH),alpha)
CXX = cxx
ifeq ($(CXX),cxx)
OPT ?= -fast
CXXFLAGS += -D__USE_STD_IOSTREAM -arch host -msg_display_number -msg_disable 767 -msg_disable 186
endif
endif
##############################

##############################
# GNU compiler
ifeq ($(CXX),g++)
CXXFLAGS += -Wall -Wno-sign-compare
CFLAGS = $(CXXFLAGS)
endif

LIBS=-lm
#LIBS += -lparse_input
CXXFLAGS += $(INCLUDES) $(OPT)
#CFLAGS += $(INCLUDES) $(OPT) Why is this not needed?
CC=cc
##############################

##############################
# Define sources
PROGS_INDEP = varigram_kn counts2kn bin2arpa arpa2bin arpa2arpa arpasize 
#find_ngram_counts merge_ngram_counts take_n_most_common_words 
PROGS_DEP = perplexity check_model add_zeroprob_grams 
#moc_writetest prune varicluster  
CLASS_SRCS_INDEP = QFit.cc Vocabulary.cc Endian.cc GramSorter.cc conf.cc io.cc str.cc
CLASS_SRCS_DEP = TreeGram.cc TreeGramArpaReader.cc PerplexityFuncs.cc 
TEMPLATE_HEADERS =Storage.hh NgramCounts.hh MultiOrderCounts.hh InterKn.hh VarigramFuncs.hh HashGram.hh
C_OBJS = sparse_matrix.o
#PROGS_DEP += qfit_test clwritetest hashgramtest
HEADERS= def.hh sikMatrix.hh NGram.hh

CLASS_SRCS = $(CLASS_SRCS_INDEP) $(CLASS_SRCS_DEP)
TEMPLATE_SRCS = $(TEMPLATE_HEADERS:.hh=_tmpl.hh) $(TEMPLATE_HEADERS)
CLASS_HEADERS_INDEP = $(CLASS_SRCS_INDEP:.cc=.hh) $(HEADERS)
CLASS_HEADERS = $(CLASS_SRCS:.cc=.hh) $(TEMPLATE_SRCS) $(HEADERS)

PROGS = $(PROGS_INDEP) $(PROGS_DEP)
CLASS_OBJS_INDEP = $(CLASS_SRCS_INDEP:.cc=.o)
CLASS_OBJS_PARTDEP = $(CLASS_SRCS_PARTDEP:.cc=.o) $(CLASS_SRCS_DEP:.cc=.o)
CLASS_OBJS = $(CLASS_OBJS_INDEP) $(CLASS_OBJS_PARTDEP) 
C_SRCS = $(C_OBJS:.o=.c)
C_HEADERS = $(C_OBJS:.o=.h) matrix_common.h

################################################################
# Compilation commands
all: $(PROGS) $(C_OBJS) $(PROGS_FMAP) libvarikn.a swig

libvarikn.a: $(C_OBJS) $(CLASS_OBJS) 
	ar r libvarikn.a $(C_OBJS) $(CLASS_OBJS)

swig:
	(cd swig; make)

$(PROGS): %: %.cc $(CLASS_OBJS) $(C_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@ $(CLASS_OBJS) $(C_OBJS) $(LIBS)

$(CLASS_OBJS_INDEP): %.o: %.cc $(CLASS_HEADERS_INDEP)
	$(CXX) -c $(CXXFLAGS) $(LDFLAGS) $< -o $@

$(CLASS_OBJS_PARTDEP): %.o: %.cc $(CLASS_HEADERS) $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $(LDFLAGS) $< -o $@

$(C_OBJS): $(C_SRCS) $(C_HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean swig

clean:
	rm -f $(CLASS_OBJS) $(PROGS) $(C_OBJS)
	rm -rf ii_files