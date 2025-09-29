TOP_DIR := $(abspath ..)
INCLUDE_DIR := $(TOP_DIR)/include
SRC_DIR := $(TOP_DIR)/src
BUILD_DIR := $(TOP_DIR)/build
OBJ_DIR := $(BUILD_DIR)/obj
LIB_DIR := $(BUILD_DIR)/lib

SCRIPTDIR   := $(TOP_DIR)/scripts
VERSION     ?= $(shell cd $(TOP_DIR) && git describe --tags --always --dirty 2>/dev/null || echo 0.0.0)
GIT_REV     := $(shell cd $(TOP_DIR) && git rev-parse --short HEAD 2>/dev/null || echo unknown)
# Extract c++ standard from CXXFLAGS or default:
CXX_STD_STR ?= $(or $(patsubst -std=%,%,$(filter -std=%,$(CXXFLAGS))),c++17)
USE_ROOT    ?= yes   # your project uses ROOT; set 'no' if ever optional

CONFIG_IN   := $(SCRIPTDIR)/rarexsec-config.in
CONFIG_OUT  := $(BUILD_DIR)/bin/rarexsec-config
ROOT_WRAPPER:= $(SCRIPTDIR)/rarexsec-root.sh
SETUP_MACRO := $(SCRIPTDIR)/setup_rarexsec.C

PROJECT_SHARED_LIB_NAME := rarexsec
ROOT_SHARED_LIB_NAME := rarexsec_root
SHARED_LIB := $(LIB_DIR)/lib$(PROJECT_SHARED_LIB_NAME)
ROOT_SHARED_LIB := $(LIB_DIR)/lib$(ROOT_SHARED_LIB_NAME)

MKDIR := mkdir -p
RM := rm -f
RMDIR := rm -rf
AR ?= ar

CXX_STD ?= c++17
UNKNOWN_REV := unknown version

ROOTCXX := $(shell root-config --cxx 2> /dev/null)
CXX ?= $(if $(ROOTCXX),$(ROOTCXX),g++)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SOEXT := dylib
  SHARED_LDFLAGS := -dynamiclib
else ifeq ($(UNAME_S),Linux)
  SOEXT := so
  SHARED_LDFLAGS := -shared
else
  $(warning Unrecognized operating system encountered.)
  SOEXT := so
  SHARED_LDFLAGS := -shared
endif

GIT := $(shell command -v git 2> /dev/null)

ifndef GIT
  GIT_REVISION := $(UNKNOWN_REV)
else
  GIT_REVPARSE_CODE := $(shell $(GIT) -C $(TOP_DIR) rev-parse 2> /dev/null && echo "$$?")
  ifeq ($(GIT_REVPARSE_CODE),0)
    GIT_REVISION := $(shell $(GIT) -C $(TOP_DIR) rev-parse --short HEAD)
    GIT_DIFF_INDEX_CODE := $(shell $(GIT) -C $(TOP_DIR) diff-index --quiet HEAD 2> /dev/null; echo "$$?")
    ifneq ($(GIT_DIFF_INDEX_CODE),0)
      GIT_REVISION := $(GIT_REVISION)-dirty
    endif
  else
    GIT_REVISION := $(UNKNOWN_REV)
  endif
endif

ifneq (,$(wildcard $(TOP_DIR)/.VERSION))
  PROJECT_VERSION := $(shell cat $(TOP_DIR)/.VERSION)
else
  PROJECT_VERSION := $(VERSION)
endif

GIT_REV := $(GIT_REVISION)

ROOTCONFIG := $(shell command -v root-config 2> /dev/null)
ROOT := $(shell command -v root 2> /dev/null)
ROOTCLING := $(shell command -v rootcling 2> /dev/null)
ifndef ROOTCLING
  ROOTCLING := $(shell command -v rootcint 2> /dev/null)
endif

ROOT_PREFIX :=
ifeq ($(ROOTCONFIG),)
  $(warning ROOT not found. Dictionary library will not be built.)
  USE_ROOT := no
else ifeq ($(ROOTCLING),)
  $(warning ROOT dictionary tool (rootcling/rootcint) not found. Dictionary library will not be built.)
  USE_ROOT := no
else
  ROOT_VERSION := $(shell $(ROOTCONFIG) --version)
  $(info Found ROOT version $(ROOT_VERSION) in $(ROOT))
  USE_ROOT := yes
  ROOT_PREFIX := $(shell $(ROOTCONFIG) --prefix)
  ROOT_CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)
  ROOT_LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
  ROOT_LIBS := $(shell $(ROOTCONFIG) --libs)
endif

EXTRA_INC ?=
EXTRA_CXXFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=

NLOHMANN_JSON_INC ?=
ifeq ($(strip $(NLOHMANN_JSON_INC)),)
  ROOT_BIN_INCLUDE := $(if $(ROOTCONFIG),$(abspath $(dir $(ROOTCONFIG))../include))
  NLOHMANN_JSON_SEARCH_DIRS := \
    $(if $(ROOT_PREFIX),$(ROOT_PREFIX)/include) \
    $(ROOT_BIN_INCLUDE) \
    /usr/include \
    /usr/local/include
  NLOHMANN_JSON_HEADER := $(firstword $(foreach d,$(filter-out ,$(NLOHMANN_JSON_SEARCH_DIRS)),$(wildcard $(d)/nlohmann/json.hpp)))
  ifneq ($(NLOHMANN_JSON_HEADER),)
    EXTRA_INC += -I$(patsubst %/nlohmann/json.hpp,%,$(NLOHMANN_JSON_HEADER))
  endif
else
  EXTRA_INC += -I$(strip $(NLOHMANN_JSON_INC))
endif

OPTFLAGS := -O3
ifeq ($(MAKECMDGOALS),debug)
  OPTFLAGS := -O0 -g
endif

USER_CXXFLAGS := $(CXXFLAGS)
USER_LDFLAGS := $(LDFLAGS)

BASE_CXXFLAGS := -std=$(CXX_STD) -Wall -Wextra -Wpedantic -fPIC -I$(INCLUDE_DIR) \
                  $(ROOT_CXXFLAGS) $(EXTRA_INC) $(EXTRA_CXXFLAGS) \
                  -DRAREXSEC_VERSION=\"$(PROJECT_VERSION)\"
CXXFLAGS := $(OPTFLAGS) $(BASE_CXXFLAGS) $(USER_CXXFLAGS)

LINK_FLAGS := $(ROOT_LDFLAGS) $(EXTRA_LDFLAGS) $(USER_LDFLAGS)
LDLIBS := $(ROOT_LIBS) $(EXTRA_LIBS)

SRCS_CC := $(shell find $(SRC_DIR) -type f -name '*.cc' 2> /dev/null)
SRCS_CPP := $(shell find $(SRC_DIR) -type f -name '*.cpp' 2> /dev/null)
SRCS := $(SRCS_CC) $(SRCS_CPP)

OBJS := $(SRCS_CC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o) \
        $(SRCS_CPP:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OBJS := $(filter-out $(OBJ_DIR)/$(ROOT_SHARED_LIB_NAME)_dict.o,$(OBJS))

ifeq ($(strip $(SRCS)),)
  $(warning No sources found under '$(SRC_DIR)'.)
endif

ROOT_DICT_HEADERS_REL := rarexsec/EventProcessor.h rarexsec/PreSelection.h rarexsec/MuonSelector.h rarexsec/TruthClassifier.h
ROOT_DICT_HEADERS := $(addprefix $(INCLUDE_DIR)/,$(ROOT_DICT_HEADERS_REL))
ROOT_DICT_SRC := $(OBJ_DIR)/$(ROOT_SHARED_LIB_NAME)_dict.cxx
ROOT_DICT_OBJ := $(OBJ_DIR)/$(ROOT_SHARED_LIB_NAME)_dict.o
ROOT_DICT_PCM := $(OBJ_DIR)/$(ROOT_SHARED_LIB_NAME)_dict_rdict.pcm

ALL_TARGETS := $(SHARED_LIB).$(SOEXT) $(LIB_DIR)/lib$(PROJECT_SHARED_LIB_NAME).a
ifeq ($(USE_ROOT),yes)
  ALL_TARGETS += $(ROOT_SHARED_LIB).$(SOEXT)
endif
ALL_TARGETS += $(CONFIG_OUT)

.PHONY: all debug shared static rootlib install run clean distclean print

all: $(ALL_TARGETS)

shared: $(SHARED_LIB).$(SOEXT)

static: $(LIB_DIR)/lib$(PROJECT_SHARED_LIB_NAME).a

rootlib: $(ROOT_SHARED_LIB).$(SOEXT)

debug: $(ALL_TARGETS)

$(LIB_DIR) $(OBJ_DIR):
	@$(MKDIR) $@

$(BUILD_DIR)/bin:
	@$(MKDIR) $@

$(CONFIG_OUT): $(CONFIG_IN) | $(BUILD_DIR)/bin
	sed -e 's|@@VERSION@@|$(VERSION)|g' \
	    -e 's|@@GIT_REVISION@@|$(GIT_REV)|g' \
	    -e 's|@@CXX_STD@@|$(CXX_STD_STR)|g' \
	    -e 's|@@USE_ROOT@@|$(USE_ROOT)|g' \
	    $< > $@
	chmod +x $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	@$(MKDIR) $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@$(MKDIR) $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(SHARED_LIB).$(SOEXT): $(OBJS) | $(LIB_DIR)
	$(CXX) $(SHARED_LDFLAGS) -o $@ $^ $(LINK_FLAGS) $(LDLIBS)
	@echo "Built $@"

$(LIB_DIR)/lib$(PROJECT_SHARED_LIB_NAME).a: $(OBJS) | $(LIB_DIR)
	$(AR) rcs $@ $^
	@echo "Built $@"

ifeq ($(USE_ROOT),yes)
$(ROOT_DICT_OBJ): $(SRC_DIR)/rarexsecLinkDef.h $(ROOT_DICT_HEADERS) | $(OBJ_DIR)
	@$(MKDIR) $(dir $@)
	$(RM) $(ROOT_DICT_SRC) $(ROOT_DICT_PCM)
	$(ROOTCLING) -f $(ROOT_DICT_SRC) -c -I$(INCLUDE_DIR) \
	$(ROOT_DICT_HEADERS_REL) \
	$(SRC_DIR)/rarexsecLinkDef.h
	$(CXX) $(CXXFLAGS) -MMD -MP -c $(ROOT_DICT_SRC) -o $@

$(ROOT_SHARED_LIB).$(SOEXT): $(SHARED_LIB).$(SOEXT) $(ROOT_DICT_OBJ) | $(LIB_DIR)
	$(CXX) $(SHARED_LDFLAGS) -o $@ $(ROOT_DICT_OBJ) -L$(LIB_DIR) -l$(PROJECT_SHARED_LIB_NAME) $(LINK_FLAGS) $(LDLIBS)
	@echo "Built $@"
endif

PREFIX ?= $(TOP_DIR)/install
INSTALL_LIBDIR := $(PREFIX)/lib
INSTALL_INCDIR := $(PREFIX)/include

install: all
	@$(MKDIR) $(INSTALL_LIBDIR) $(INSTALL_INCDIR) $(PREFIX)/bin $(PREFIX)/scripts
	@cp -a $(SHARED_LIB).$(SOEXT) $(INSTALL_LIBDIR)/
	@cp -a $(LIB_DIR)/lib$(PROJECT_SHARED_LIB_NAME).a $(INSTALL_LIBDIR)/
ifeq ($(USE_ROOT),yes)
	@cp -a $(ROOT_SHARED_LIB).$(SOEXT) $(INSTALL_LIBDIR)/
endif
	@rsync -a --delete $(INCLUDE_DIR)/ $(INSTALL_INCDIR)/
	@cp -a $(CONFIG_OUT) $(PREFIX)/bin/
	@cp -a $(ROOT_WRAPPER) $(PREFIX)/bin/rarexsec-root
	@cp -a $(SETUP_MACRO) $(PREFIX)/scripts/
	@echo "Installed to $(PREFIX)"
	@echo "  - libs:    $(INSTALL_LIBDIR)"
	@echo "  - headers: $(INSTALL_INCDIR)"
	@echo "  - scripts: $(PREFIX)/bin"
	@echo "  - macros:  $(PREFIX)/scripts"

run: $(SHARED_LIB).$(SOEXT)
	@test -n "$(ROOT)" || (echo "ROOT executable not found" && exit 1)
	@root -l -q -e 'gSystem->Load("$(abspath $(SHARED_LIB)).$(SOEXT)"); \
                         gSystem->AddIncludePath("-I$(abspath $(INCLUDE_DIR))"); \
                         gROOT->ProcessLine(".x $(SETUP_MACRO)$(ARGS)");'

print:
	@echo "CXX           = $(CXX)"
	@echo "CXXFLAGS      = $(CXXFLAGS)"
	@echo "LINK_FLAGS    = $(LINK_FLAGS)"
	@echo "LDLIBS        = $(LDLIBS)"
	@echo "SRCS          = $(SRCS)"
	@echo "OBJS          = $(OBJS)"
	@echo "LIB_DIR       = $(LIB_DIR)"
	@echo "SOEXT         = $(SOEXT)"

clean:
	$(RMDIR) $(OBJ_DIR) $(LIB_DIR) $(BUILD_DIR)/bin

distclean: clean
	$(RMDIR) $(PREFIX)

DEPS := $(OBJS:.o=.d)
ifeq ($(USE_ROOT),yes)
DEPS += $(ROOT_DICT_OBJ:.o=.d)
endif

-include $(DEPS)
