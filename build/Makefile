TOPDIR := $(shell cd .. && pwd)
INCDIR := $(TOPDIR)/include
SRC_DIR := $(TOPDIR)/src
BUILDDIR := $(TOPDIR)/build
LIBDIR := $(BUILDDIR)/lib
BINDIR := $(BUILDDIR)/bin

SCRIPTDIR   := $(TOPDIR)/scripts
VERSION     ?= $(shell git describe --tags --always --dirty 2>/dev/null || echo 0.0.0)
GIT_REV     := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
CXX_STD     ?= c++17
USE_ROOT    ?= yes

CONFIG_IN   := $(SCRIPTDIR)/rarexsec-config.in
CONFIG_OUT  := $(BINDIR)/rarexsec-config
ROOT_WRAPPER:= $(SCRIPTDIR)/rarexsec-root.sh
SETUP_MACRO := $(SCRIPTDIR)/setup_rarexsec.C

MKDIR := mkdir -p
PREFIX ?= /usr/local
INSTALL_LIBDIR := $(PREFIX)/lib
INSTALL_INCDIR := $(PREFIX)/include

UNKNOWN_REV := unknown version

ifeq ($(MAKECMDGOALS),debug)
  override CXXFLAGS := -std=$(CXX_STD) $(CXXFLAGS) -O0 -g
else
  override CXXFLAGS := -O3 -std=$(CXX_STD) $(CXXFLAGS)
endif

CXX_STD_STR := $(or $(patsubst -std=%,%,$(filter -std=%,$(CXXFLAGS))),$(CXX_STD))

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SOEXT := dylib
else ifeq ($(UNAME_S),Linux)
  SOEXT := so
else
  $(warning Unrecognized operating system encountered.)
  SOEXT := so
endif

GIT := $(shell command -v git 2> /dev/null)

ifndef GIT
  GIT_REVISION := $(UNKNOWN_REV)
else
  GIT_REVPARSE_CODE := $(shell git rev-parse 2> /dev/null && echo "$$?")
  ifeq ($(GIT_REVPARSE_CODE),0)
    GIT_REVISION := $(shell git rev-parse --short HEAD)
    GIT_DIFF_INDEX_CODE := $(shell git diff-index --quiet HEAD 2> /dev/null && echo "$$?")
    ifneq ($(GIT_DIFF_INDEX_CODE),0)
      GIT_REVISION := $(GIT_REVISION)-dirty
    endif
  else
    GIT_REVISION := $(UNKNOWN_REV)
  endif
endif

ifneq (,$(wildcard ../.VERSION))
  PROJECT_VERSION := $(shell cat ../.VERSION)
else
  PROJECT_VERSION := $(GIT_REVISION)
endif

VERSION := $(PROJECT_VERSION)
GIT_REV := $(GIT_REVISION)

override CXXFLAGS += -I$(INCDIR) -Wall -Wextra -Wpedantic
override CXXFLAGS += -DRAREXSEC_VERSION=\"$(PROJECT_VERSION)\"

LIBNAME := rarexsec
ROOT_SHARED_LIB_NAME := rarexsec_root

CXX ?= g++
AR ?= ar
LDFLAGS ?=

SOURCES := $(wildcard $(SRC_DIR)/*.cc)
OBJECTS := $(notdir $(SOURCES:.cc=.o))
OBJECTS := $(filter-out $(ROOT_SHARED_LIB_NAME)_dict.o,$(OBJECTS))

ROOTCONFIG := $(shell command -v root-config 2> /dev/null)
ROOT := $(shell command -v root 2> /dev/null)
ROOTCLING := $(shell command -v rootcling 2> /dev/null)
ifndef ROOTCLING
  ROOTCLING := $(shell command -v rootcint 2> /dev/null)
endif

ifeq ($(ROOTCONFIG),)
  $(info ROOT not found. Building without dictionary support.)
  USE_ROOT := no
else
  ROOT_VERSION := $(shell $(ROOTCONFIG) --version)
  $(info Found ROOT version $(ROOT_VERSION) in $(ROOT))
  USE_ROOT := yes
  ROOT_CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)
  ROOT_LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
  ROOT_LIBS := $(shell $(ROOTCONFIG) --libs)
endif

SHARED_LIB := $(LIBDIR)/lib$(LIBNAME).$(SOEXT)
STATIC_LIB := $(LIBDIR)/lib$(LIBNAME).a
ROOT_SHARED_LIB := $(LIBDIR)/lib$(ROOT_SHARED_LIB_NAME).$(SOEXT)

ALL_TARGETS := $(SHARED_LIB) $(STATIC_LIB)
ifeq ($(USE_ROOT),yes)
  ALL_TARGETS += $(ROOT_SHARED_LIB)
endif

all: $(ALL_TARGETS) $(CONFIG_OUT)

debug: $(ALL_TARGETS) $(CONFIG_OUT)

.INTERMEDIATE: $(OBJECTS)

$(LIBDIR):
	@$(MKDIR) $@

$(BINDIR):
	@$(MKDIR) $@

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -fPIC -o $@ -c $<

$(SHARED_LIB): $(OBJECTS) | $(LIBDIR)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -shared -o $@ $^ $(LDFLAGS)

$(STATIC_LIB): $(OBJECTS) | $(LIBDIR)
	$(AR) rcs $@ $^

ifeq ($(USE_ROOT),yes)
ROOT_DICT_SRC := $(ROOT_SHARED_LIB_NAME)_dict.cxx
ROOT_DICT_OBJ := $(ROOT_SHARED_LIB_NAME)_dict.o
ROOT_DICT_PCM := $(ROOT_SHARED_LIB_NAME)_dict_rdict.pcm

$(ROOT_DICT_OBJ): $(SRC_DIR)/rarexsecLinkDef.h \
        $(INCDIR)/rarexsec/EventProcessor.h \
        $(INCDIR)/rarexsec/PreSelection.h \
        $(INCDIR)/rarexsec/MuonSelector.h \
        $(INCDIR)/rarexsec/TruthClassifier.h
	$(RM) $(ROOT_DICT_SRC) $(ROOT_DICT_PCM)
	$(ROOTCLING) -f $(ROOT_DICT_SRC) -c -I$(INCDIR) \
                rarexsec/EventProcessor.h \
                rarexsec/PreSelection.h \
                rarexsec/MuonSelector.h \
                rarexsec/TruthClassifier.h \
                $(SRC_DIR)/rarexsecLinkDef.h
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -I$(INCDIR) -fPIC -o $@ -c $(ROOT_DICT_SRC)

$(ROOT_SHARED_LIB): $(SHARED_LIB) $(ROOT_DICT_OBJ) | $(LIBDIR)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -shared -o $@ $(ROOT_DICT_OBJ) -L$(LIBDIR) -l$(LIBNAME) $(ROOT_LDFLAGS) $(ROOT_LIBS)
endif

$(CONFIG_OUT): $(CONFIG_IN) | $(BINDIR)
	sed -e 's|@@VERSION@@|$(VERSION)|g' \
    -e 's|@@GIT_REVISION@@|$(GIT_REV)|g' \
    -e 's|@@CXX_STD@@|$(CXX_STD_STR)|g' \
    -e 's|@@USE_ROOT@@|$(USE_ROOT)|g' \
    $< > $@
	chmod +x $@

install: all
	@$(MKDIR) $(INSTALL_LIBDIR) $(INSTALL_INCDIR) $(PREFIX)/bin $(PREFIX)/scripts
	@cp -a $(SHARED_LIB) $(INSTALL_LIBDIR)/
	@cp -a $(STATIC_LIB) $(INSTALL_LIBDIR)/
ifeq ($(USE_ROOT),yes)
	@cp -a $(ROOT_SHARED_LIB) $(INSTALL_LIBDIR)/
endif
	@rsync -a --delete $(INCDIR)/ $(INSTALL_INCDIR)/
	@cp -a $(CONFIG_OUT) $(PREFIX)/bin/
	@cp -a $(ROOT_WRAPPER) $(PREFIX)/bin/rarexsec-root
	@cp -a $(SETUP_MACRO) $(PREFIX)/scripts/
	@echo "Installed helper scripts to $(PREFIX)/bin and macro to $(PREFIX)/scripts"

clean:
	$(RM) $(OBJECTS) $(STATIC_LIB) $(SHARED_LIB) $(CONFIG_OUT)
ifeq ($(USE_ROOT),yes)
	$(RM) $(ROOT_DICT_SRC) $(ROOT_DICT_OBJ) $(ROOT_DICT_PCM) $(ROOT_SHARED_LIB)
else
	$(RM) $(ROOT_SHARED_LIB)
endif

.PHONY: all debug clean install
