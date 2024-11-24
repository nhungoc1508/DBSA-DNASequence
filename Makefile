EXTENSION	= dnasequence
MODULES		= dnasequence
DATA		= dnasequence--1.0.sql dnasequence.control

PG_CPPFLAGS = -I$(shell brew --prefix icu4c)/include
SHLIB_LINK = -L$(shell brew --prefix icu4c)/lib -licui18n -licuuc -licudata

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)