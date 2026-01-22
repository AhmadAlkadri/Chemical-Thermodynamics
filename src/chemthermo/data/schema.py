"""Schema details for the packaged component databank."""

SCHEMA_VERSION = 1

TOP_LEVEL_KEYS = ("schema_version", "components")
COMPONENT_KEYS = ("name", "formula", "properties")
PROPERTY_KEYS = ("MW_kg_per_mol", "Tc_K", "Pc_Pa", "omega")
ANTOINE_KEYS = ("A", "B", "C", "Tmin_K", "Tmax_K")
