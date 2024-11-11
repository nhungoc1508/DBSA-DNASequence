# Makefile tells WSL how to compile the extension code
# Extension defines extension name
EXTENSION = dna_extension
# Modules tells postgreSQL the c file to compile
MODULES = dna_extension
# Data specifies the sql scripts
DATA = dna_extension--1.0.sql

# we need to always provide PG_CONFIG when building or testing the extension
# it tells Makefile where to find the configurations for postgreSQL
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

