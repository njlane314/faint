TOP := $(abspath .)
SRC := $(TOP)/src
INC := $(TOP)/include
BUILD := $(TOP)/build
OBJ := $(BUILD)/obj
LIB := $(BUILD)/lib
NAME := rarexsec

UNAME := $(shell uname -s)
SOEXT := $(if $(filter Darwin,$(UNAME)),dylib,so)
SHAREDFLAGS := $(if $(filter Darwin,$(UNAME)),-dynamiclib,-shared)

CXX ?= $(shell root-config --cxx)
CPPFLAGS += -I$(INC) $(shell root-config --cflags) $(NLOHMANN_JSON_CFLAGS)
CXXFLAGS += -O3 -std=c++17 -Wall -Wextra -Wpedantic -fPIC
LDFLAGS  += $(shell root-config --ldflags)
LDLIBS   += $(shell root-config --libs)

SRCS := $(shell find $(SRC) -type f \( -name '*.cc' -o -name '*.cpp' \) 2>/dev/null)
OBJS := $(patsubst $(SRC)/%.cc,$(OBJ)/%.o,$(filter %.cc,$(SRCS))) \
        $(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(filter %.cpp,$(SRCS)))

SHARED := $(LIB)/lib$(NAME).$(SOEXT)

all: $(SHARED)

$(OBJ)/%.o: $(SRC)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
$(OBJ)/%.o: $(SRC)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(SHARED): $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(SHAREDFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -rf $(BUILD)

DEPS := $(OBJS:.o=.d)
-include $(DEPS)
