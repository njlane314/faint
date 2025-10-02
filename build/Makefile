TOP:=$(abspath ..)
INC:=$(TOP)/include
SRC:=$(TOP)/src
BUILD:=$(TOP)/build
OBJ:=$(BUILD)/obj
LIB:=$(BUILD)/lib
BIN:=$(BUILD)/bin
NAME:=rarexsec

ROOTCONFIG:=$(shell command -v root-config 2>/dev/null)
ifeq ($(ROOTCONFIG),)
$(error root-config not found in PATH)
endif

UNAME:=$(shell uname -s)
ifeq ($(UNAME),Darwin)
SOEXT:=dylib
SHAREDFLAGS:=-dynamiclib
else
SOEXT:=so
SHAREDFLAGS:=-shared
endif

CXX:=$(shell root-config --cxx)
AR?=ar

ifneq (,$(wildcard $(TOP)/.VERSION))
VERSION:=$(shell cat $(TOP)/.VERSION)
else
VERSION:=$(shell cd $(TOP) && (git describe --tags --always --dirty 2>/dev/null || echo 0.0.0))
endif
GIT_REV:=$(shell cd $(TOP) && (git rev-parse --short HEAD 2>/dev/null || echo unknown))

CXXSTD?=c++17
OPT?=-O3

CPPFLAGS+=-I$(INC) -DRAREXSEC_VERSION=\"$(VERSION)\"
CXXFLAGS+=-Wall -Wextra -Wpedantic -fPIC $(OPT)
ifneq ($(findstring -std=,$(CXXFLAGS)),)
else
CXXFLAGS+=-std=$(CXXSTD)
endif
CPPFLAGS+=$(shell $(ROOTCONFIG) --cflags)

define _split_paths
$(strip $(foreach path,$(subst :, ,$(1)),$(if $(path),$(path))))
endef

JSON_INCLUDE_FLAGS:=$(strip $(NLOHMANN_JSON_CFLAGS))
ifeq ($(JSON_INCLUDE_FLAGS),)
JSON_INCLUDE_DIRS:=$(call _split_paths,$(NLOHMANN_JSON_INC) $(NLOHMANN_JSON_INCLUDE_DIRS))
JSON_INCLUDE_DIRS+=$(addsuffix /include,$(call _split_paths,$(NLOHMANN_JSON_FQ_DIR) $(NLOHMANN_JSON_DIR)))
JSON_INCLUDE_FLAGS:=$(strip $(foreach dir,$(JSON_INCLUDE_DIRS),-I$(dir)))
endif

CPPFLAGS+=$(JSON_INCLUDE_FLAGS)
LDFLAGS+=$(shell $(ROOTCONFIG) --ldflags)
LDLIBS+=$(shell $(ROOTCONFIG) --libs)

SRCS:=$(shell find $(SRC) -type f \( -name '*.cc' -o -name '*.cpp' \) 2>/dev/null)
OBJS:=$(SRCS:$(SRC)/%.cc=$(OBJ)/%.o)
OBJS:=$(OBJS:$(SRC)/%.cpp=$(OBJ)/%.o)

SHARED:=$(LIB)/lib$(NAME).$(SOEXT)
STATIC:=$(LIB)/lib$(NAME).a

CONFIG_IN:=$(TOP)/rarexsec-config.in
CONFIG_OUT:=$(BIN)/rarexsec-config

.DEFAULT_GOAL:=all
.PHONY:all debug clean distclean install print run

all: $(SHARED) $(STATIC) $(CONFIG_OUT)

debug: OPT:=-O0 -g
debug: all

$(OBJ) $(LIB) $(BIN):
	@mkdir -p $@

$(OBJ)/%.o: $(SRC)/%.cc | $(OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(OBJ)/%.o: $(SRC)/%.cpp | $(OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(SHARED): $(OBJS) | $(LIB)
	$(CXX) $(SHAREDFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)
	@echo Built $@

$(STATIC): $(OBJS) | $(LIB)
	$(AR) rcs $@ $^
	@echo Built $@

$(CONFIG_OUT): $(CONFIG_IN) | $(BIN)
	sed -e 's|@@VERSION@@|$(VERSION)|g' \
	-e 's|@@GIT_REVISION@@|$(GIT_REV)|g' \
	-e 's|@@CXX_STD@@|$(CXXSTD)|g' \
	-e 's|@@USE_ROOT@@|yes|g' \
	$< > $@
	chmod +x $@

PREFIX?=$(TOP)/install
INSTALL_LIB:=$(PREFIX)/lib
INSTALL_INC:=$(PREFIX)/include
INSTALL_BIN:=$(PREFIX)/bin
INSTALL_SCRIPTS:=$(PREFIX)/scripts

install: all
	@mkdir -p $(INSTALL_LIB) $(INSTALL_INC) $(INSTALL_BIN) $(INSTALL_SCRIPTS)
	@cp -a $(LIB)/* $(INSTALL_LIB)/
	@rsync -a --delete $(INC)/ $(INSTALL_INC)/
	@cp -a $(CONFIG_OUT) $(INSTALL_BIN)/
        @[ -f $(TOP)/rarexsec-root.sh ] && chmod +x $(TOP)/rarexsec-root.sh || true
        @[ -f $(TOP)/rarexsec-root.sh ] && cp -a $(TOP)/rarexsec-root.sh $(INSTALL_BIN)/rarexsec-root || true
	@[ -f $(TOP)/setup_rarexsec.C ] && cp -a $(TOP)/setup_rarexsec.C $(INSTALL_SCRIPTS)/ || true
	@echo Installed to $(PREFIX)

run: $(SHARED)
	@command -v root >/dev/null 2>&1 || (echo ROOT executable not found && exit 1)
	@root -l -q -e 'gSystem->Load("$(abspath $(SHARED))"); gSystem->AddIncludePath("-I$(abspath $(INC))"); gROOT->ProcessLine(".x $(TOP)/setup_rarexsec.C$(ARGS)");'

print:
	@echo CXX=$(CXX)
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo LDFLAGS=$(LDFLAGS)
	@echo LDLIBS=$(LDLIBS)
	@echo SRCS=$(SRCS)
	@echo OBJS=$(OBJS)
	@echo SHARED=$(SHARED)
	@echo STATIC=$(STATIC)
	@echo SOEXT=$(SOEXT)
	@echo JSON_INCLUDE_FLAGS=$(JSON_INCLUDE_FLAGS)

clean:
	@rm -rf $(OBJ) $(LIB) $(BIN)

distclean: clean
	@rm -rf $(PREFIX)

DEPS:=$(OBJS:.o=.d)
-include $(DEPS)
