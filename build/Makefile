TOP_DIR := $(shell cd .. && pwd)
INCLUDE_DIR := $(TOP_DIR)/include
SRC_DIR := $(TOP_DIR)/src
BUILD_DIR := $(TOP_DIR)/build

CXX_STD ?= c++17
UNKNOWN_REV := unknown version

ifeq ($(MAKECMDGOALS),debug)
  override CXXFLAGS := -std=$(CXX_STD) $(CXXFLAGS) -O0 -g
else
  override CXXFLAGS := -O3 -std=$(CXX_STD) $(CXXFLAGS)
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHARED_LIB_SUFFIX := dylib
else ifeq ($(UNAME_S),Linux)
  SHARED_LIB_SUFFIX := so
else
  $(warning Unrecognized operating system encountered.)
  SHARED_LIB_SUFFIX := so
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

override CXXFLAGS += -DRAREXSEC_VERSION=\"$(PROJECT_VERSION)\"

PROJECT_SHARED_LIB_NAME := rarexsec
SHARED_LIB := lib$(PROJECT_SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)
ROOT_SHARED_LIB_NAME := rarexsec_root
ROOT_SHARED_LIB := lib$(ROOT_SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)

CXX ?= g++
override CXXFLAGS += -I$(INCLUDE_DIR) -Wall -Wextra -Wpedantic
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

ALL_TARGETS := $(SHARED_LIB)
ifeq ($(USE_ROOT),yes)
  ALL_TARGETS += $(ROOT_SHARED_LIB)
endif

all: $(ALL_TARGETS)

debug: $(ALL_TARGETS)

.INTERMEDIATE: $(OBJECTS)

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -fPIC -o $@ -c $<

$(SHARED_LIB): $(OBJECTS)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -shared -o $@ $^ $(LDFLAGS)

ifeq ($(USE_ROOT),yes)
ROOT_DICT_SRC := $(ROOT_SHARED_LIB_NAME)_dict.cxx
ROOT_DICT_OBJ := $(ROOT_SHARED_LIB_NAME)_dict.o
ROOT_DICT_PCM := $(ROOT_SHARED_LIB_NAME)_dict_rdict.pcm


$(ROOT_DICT_OBJ): $(SRC_DIR)/rarexsecLinkDef.h \
	$(INCLUDE_DIR)/rarexsec/EventProcessor.h \
	$(INCLUDE_DIR)/rarexsec/PreSelection.h \
	$(INCLUDE_DIR)/rarexsec/MuonSelector.h \
	$(INCLUDE_DIR)/rarexsec/TruthClassifier.h
	$(RM) $(ROOT_DICT_SRC) $(ROOT_DICT_PCM)
        $(ROOTCLING) -f $(ROOT_DICT_SRC) -c -I$(INCLUDE_DIR) \
                rarexsec/EventProcessor.h \
                rarexsec/PreSelection.h \
                rarexsec/MuonSelector.h \
                rarexsec/TruthClassifier.h \
                $(SRC_DIR)/rarexsecLinkDef.h
        $(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -I$(INCLUDE_DIR) -fPIC -o $@ -c $(ROOT_DICT_SRC)

$(ROOT_SHARED_LIB): $(SHARED_LIB) $(ROOT_DICT_OBJ)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -shared -o $@ $(ROOT_DICT_OBJ) -L. -l$(PROJECT_SHARED_LIB_NAME) $(ROOT_LDFLAGS) $(ROOT_LIBS)
endif

clean:
	$(RM) $(SHARED_LIB) $(ROOT_SHARED_LIB) $(OBJECTS)
ifeq ($(USE_ROOT),yes)
	$(RM) $(ROOT_DICT_SRC) $(ROOT_DICT_OBJ) $(ROOT_DICT_PCM)
endif

.PHONY: all debug clean
