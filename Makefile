EXTENSION	= dnasequence
MODULES		= dnasequence
DATA		= dnasequence--1.0.sql dnasequence.control

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)