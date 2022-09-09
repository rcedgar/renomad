BINDIR := ../bin
OBJDIR := o
BINPATH := $(BINDIR)/renomad

CPPFLAGS := $(CPPFLAGS) -DNDEBUG -pthread

CXX = g++
CXXFLAGS := $(CXXFLAGS) -O3 -fopenmp -ffast-math -msse -mfpmath=sse

UNAME_S := $(shell uname -s)
LDFLAGS := $(LDFLAGS) -O3 -fopenmp -pthread -lpthread
ifeq ($(UNAME_S),Linux)
    LDFLAGS += -static
endif

HDRS = \
  alpha.h \
  atab.h \
  cmds.h \
  countsort.h \
  fastaseqsource.h \
  fastq.h \
  fastqseqsource.h \
  fileseqsource.h \
  filetype.h \
  gobuff.h \
  joiner.h \
  joinera.h \
  joinercv.h \
  kmc_file.h \
  kmer_api.h \
  kmer_defs.h \
  linereader.h \
  lockobj.h \
  lockobjs.h \
  mmer.h \
  mxreader.h \
  myopts.h \
  myutils.h \
  nomad.h \
  nomadsim.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  omplock.h \
  params.h \
  seqdb.h \
  seqinfo.h \
  seqsource.h \
  sort.h \

OBJS = \
  $(OBJDIR)/alpha.o \
  $(OBJDIR)/atab.o \
  $(OBJDIR)/dump.o \
  $(OBJDIR)/fastaseqsource.o \
  $(OBJDIR)/fastq.o \
  $(OBJDIR)/fastqseqsource.o \
  $(OBJDIR)/fileseqsource.o \
  $(OBJDIR)/filetype.o \
  $(OBJDIR)/join.o \
  $(OBJDIR)/joiner.o \
  $(OBJDIR)/joinera.o \
  $(OBJDIR)/joinercv.o \
  $(OBJDIR)/joinp.o \
  $(OBJDIR)/joinpvc.o \
  $(OBJDIR)/kmc_file.o \
  $(OBJDIR)/kmer_api.o \
  $(OBJDIR)/linereader.o \
  $(OBJDIR)/lockobj.o \
  $(OBJDIR)/mmer.o \
  $(OBJDIR)/mxreader.o \
  $(OBJDIR)/nomad.o \
  $(OBJDIR)/nomad_bench.o \
  $(OBJDIR)/nomadsim.o \
  $(OBJDIR)/nomad_sim_bench.o \
  $(OBJDIR)/objmgr.o \
  $(OBJDIR)/paper_example.o \
  $(OBJDIR)/params.o \
  $(OBJDIR)/pval.o \
  $(OBJDIR)/pvaldoanchor.o \
  $(OBJDIR)/pvalue_brute.o \
  $(OBJDIR)/rekmc_main.o \
  $(OBJDIR)/myutils.o \
  $(OBJDIR)/repval.o \
  $(OBJDIR)/seqdb.o \
  $(OBJDIR)/seqdbfromfasta.o \
  $(OBJDIR)/seqinfo.o \
  $(OBJDIR)/seqsource.o \

.PHONY: clean

$(BINPATH) : $(BINDIR)/ $(OBJDIR)/ $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(BINPATH)
	strip -d $(BINPATH)

$(OBJDIR)/ :
	mkdir -p $(OBJDIR)/

$(BINDIR)/ :
	mkdir -p $(BINDIR)/

$(OBJDIR)/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR)/ $(BINPATH)
