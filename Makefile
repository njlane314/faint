TOP := $(abspath .)
SRC := $(TOP)/src
BUILD := $(TOP)/build
OBJ := $(BUILD)/obj
LIB := $(BUILD)/lib
NAME := rarexsec

SOEXT := so
SHAREDFLAGS := -shared

CXX ?= $(shell root-config --cxx)
NLOHMANN_JSON_CFLAGS ?= -isystem $(NLOHMANN_JSON_INC)
CPPFLAGS += -I$(SRC) $(shell root-config --cflags) $(NLOHMANN_JSON_CFLAGS)
CXXFLAGS += -O3 -std=c++17 -Wall -Wextra -Wpedantic -fPIC
LDFLAGS  += $(shell root-config --ldflags)
LDLIBS   += $(shell root-config --libs)

SRCS := $(shell find $(SRC) -type f -name '*.cxx' 2>/dev/null)
OBJS := $(patsubst $(SRC)/%.cxx,$(OBJ)/%.o,$(SRCS))

SHARED := $(LIB)/lib$(NAME).$(SOEXT)

all: $(SHARED)

$(OBJ)/%.o: $(SRC)/%.cxx
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(SHARED): $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(SHAREDFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -rf $(BUILD)

DEPS := $(OBJS:.o=.d)
-include $(DEPS)
