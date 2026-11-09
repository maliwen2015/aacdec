CC ?= gcc
AR ?= ar
CFLAGS ?= -O2 -std=c99 -Wall -Wextra
CPPFLAGS += -I. -Isrc
LDLIBS += -lm

ifeq ($(OS),Windows_NT)
EXEEXT := .exe
endif

# Native Windows make (SHELL defaults to a bare "sh.exe" that does not exist)
# needs cmd.exe syntax, while MSYS2 / Git Bash / Linux / macOS make uses a
# POSIX shell and the POSIX tools.
ifeq ($(OS),Windows_NT)
ifneq (,$(findstring /,$(SHELL)))
WIN_CMD := 0
else
WIN_CMD := 1
endif
else
WIN_CMD := 0
endif

ifeq ($(WIN_CMD),1)
define MAKE_DIR
@if not exist "$(subst /,\,$(1))" mkdir "$(subst /,\,$(1))"
endef
else
define MAKE_DIR
@mkdir -p "$(1)"
endef
endif

BUILD := build
LIB := libaacdec.a
SRC := src/aacdec.c src/aac_core.c src/aac_usac.c src/aac_usac_ac.c src/aac_tables.c src/aac_huffman.c src/aac_ps.c src/aac_sbr.c
OBJ := $(patsubst %.c,$(BUILD)/obj/%.o,$(SRC))
EXAMPLE := aacdec_example$(EXEEXT)
TEST := tests/api_check$(EXEEXT)
TEST_OBJ := $(BUILD)/tests/api_check.o
COMPARE := tests/compare_check$(EXEEXT)
COMPARE_OBJ := $(BUILD)/tests/compare_check.o
HE_CHECK := tests/he_check$(EXEEXT)
HE_CHECK_OBJ := $(BUILD)/tests/he_check.o
USAC_CHECK := tests/usac_check$(EXEEXT)
USAC_CHECK_OBJ := $(BUILD)/tests/usac_check.o
PS_CHECK := tests/ps_check$(EXEEXT)
PS_CHECK_OBJ := $(BUILD)/tests/ps_check.o
PS_UNIT := tests/ps_unit$(EXEEXT)
LD_CHECK := tests/ld_check$(EXEEXT)
LD_CHECK_OBJ := $(BUILD)/tests/ld_check.o
ELD_CHECK := tests/eld_check$(EXEEXT)
ELD_CHECK_OBJ := $(BUILD)/tests/eld_check.o
DRC_CHECK := tests/drc_check$(EXEEXT)
DRC_CHECK_OBJ := $(BUILD)/tests/drc_check.o
CONCEAL_CHECK := tests/conceal_check$(EXEEXT)
CONCEAL_CHECK_OBJ := $(BUILD)/tests/conceal_check.o
ER_CHECK := tests/er_check$(EXEEXT)
ER_CHECK_OBJ := $(BUILD)/tests/er_check.o
CRC_CHECK := tests/crc_check$(EXEEXT)
CRC_CHECK_OBJ := $(BUILD)/tests/crc_check.o

.PHONY: all clean check

all: $(LIB) $(EXAMPLE)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

$(BUILD)/obj/%.o: %.c
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(EXAMPLE): example_decoder.c $(LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/api_check.o: tests/api_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TEST): $(TEST_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/compare_check.o: tests/compare_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(COMPARE): $(COMPARE_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/he_check.o: tests/he_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(HE_CHECK): $(HE_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/usac_check.o: tests/usac_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(USAC_CHECK): $(USAC_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/ps_check.o: tests/ps_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(PS_CHECK): $(PS_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(PS_UNIT): tests/ps_unit.c src/aac_ps.c src/aac_ps.h src/aac_ps_tables.h
	$(CC) $(CPPFLAGS) $(CFLAGS) tests/ps_unit.c $(LDLIBS) -o $@

$(BUILD)/tests/ld_check.o: tests/ld_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(LD_CHECK): $(LD_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/eld_check.o: tests/eld_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(ELD_CHECK): $(ELD_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/drc_check.o: tests/drc_check.c aacdec.h src/aac_core.h src/aac_bits.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(DRC_CHECK): $(DRC_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/conceal_check.o: tests/conceal_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(CONCEAL_CHECK): $(CONCEAL_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/er_check.o: tests/er_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(ER_CHECK): $(ER_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/tests/crc_check.o: tests/crc_check.c aacdec.h
	$(call MAKE_DIR,$(dir $@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(CRC_CHECK): $(CRC_CHECK_OBJ) $(LIB)
	$(CC) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

check: $(TEST) $(COMPARE) $(HE_CHECK) $(USAC_CHECK) $(PS_CHECK) $(PS_UNIT) $(LD_CHECK) $(ELD_CHECK) $(DRC_CHECK) $(CONCEAL_CHECK) $(ER_CHECK) $(CRC_CHECK)
	$(call RUN_EXE,$(TEST)) tests/data/tiny_stereo.aac
	$(call RUN_EXE,$(TEST)) tests/data/tiny_51.aac
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_stereo.aac tests/data/tiny_stereo.s16
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_51.aac tests/data/tiny_51.s16
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_stereo.latm tests/data/tiny_stereo.s16 32
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_stereo.adif tests/data/tiny_stereo.s16 32
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_stereo.raw tests/data/tiny_stereo.s16 32 raw tests/data/tiny_stereo.asc
	$(call RUN_EXE,$(HE_CHECK)) tests/data/he_stereo.aac
	$(call RUN_EXE,$(USAC_CHECK)) tests/data/usac_stereo.loas
	$(call RUN_EXE,$(PS_CHECK)) tests/data/he_v2.aac
	$(call RUN_EXE,$(PS_UNIT))
	$(call RUN_EXE,$(LD_CHECK)) tests/data/tiny_ld.latm tests/data/tiny_ld_ref.pcm
	$(call RUN_EXE,$(ELD_CHECK)) tests/data/tiny_eld.latm tests/data/tiny_eld_ref.pcm
	$(call RUN_EXE,$(DRC_CHECK))
	$(call RUN_EXE,$(CONCEAL_CHECK)) tests/data/tiny_stereo.aac tests/data/tiny_stereo.s16 tests/data/tiny_stereo.asc
	$(call RUN_EXE,$(ER_CHECK))
	$(call RUN_EXE,$(CRC_CHECK)) tests/data/tiny_adts_crc4.aac tests/data/tiny_adts_crc_ref.s16
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_latm_amv1.loas tests/data/tiny_latm_ref.s16 32
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_latm_2sub.loas tests/data/tiny_latm_ref.s16 32
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_adts_crc1.aac tests/data/tiny_adts_crc_ref.s16 32
	$(call RUN_EXE,$(COMPARE)) tests/data/tiny_adts_crc4.aac tests/data/tiny_adts_crc_ref.s16 32

ifeq ($(WIN_CMD),1)
define RUN_EXE
$(subst /,\,$(1))
endef
else
define RUN_EXE
./$(1)
endef
endif

clean:
ifeq ($(WIN_CMD),1)
	@if exist "$(subst /,\,$(BUILD))" rmdir /s /q "$(subst /,\,$(BUILD))"
	@if exist "$(subst /,\,$(LIB))" del /q /f "$(subst /,\,$(LIB))"
	@if exist "$(subst /,\,$(EXAMPLE))" del /q /f "$(subst /,\,$(EXAMPLE))"
	@if exist "$(subst /,\,$(TEST))" del /q /f "$(subst /,\,$(TEST))"
	@if exist "$(subst /,\,$(COMPARE))" del /q /f "$(subst /,\,$(COMPARE))"
	@if exist "$(subst /,\,$(HE_CHECK))" del /q /f "$(subst /,\,$(HE_CHECK))"
	@if exist "$(subst /,\,$(USAC_CHECK))" del /q /f "$(subst /,\,$(USAC_CHECK))"
	@if exist "$(subst /,\,$(PS_CHECK))" del /q /f "$(subst /,\,$(PS_CHECK))"
	@if exist "$(subst /,\,$(PS_UNIT))" del /q /f "$(subst /,\,$(PS_UNIT))"
	@if exist "$(subst /,\,$(LD_CHECK))" del /q /f "$(subst /,\,$(LD_CHECK))"
	@if exist "$(subst /,\,$(ELD_CHECK))" del /q /f "$(subst /,\,$(ELD_CHECK))"
	@if exist "$(subst /,\,$(DRC_CHECK))" del /q /f "$(subst /,\,$(DRC_CHECK))"
	@if exist "$(subst /,\,$(CONCEAL_CHECK))" del /q /f "$(subst /,\,$(CONCEAL_CHECK))"
	@if exist "$(subst /,\,$(ER_CHECK))" del /q /f "$(subst /,\,$(ER_CHECK))"
	@if exist "$(subst /,\,$(CRC_CHECK))" del /q /f "$(subst /,\,$(CRC_CHECK))"
else
	rm -rf $(BUILD) $(LIB) $(EXAMPLE) $(TEST) $(COMPARE) $(HE_CHECK) \
	       $(USAC_CHECK) $(PS_CHECK) $(PS_UNIT) $(LD_CHECK) $(ELD_CHECK) \
	       $(DRC_CHECK) $(CONCEAL_CHECK) $(ER_CHECK) $(CRC_CHECK)
endif
