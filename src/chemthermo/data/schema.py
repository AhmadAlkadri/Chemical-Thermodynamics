"""Schema details for the packaged component databank."""

SCHEMA_VERSION = 1

TOP_LEVEL_KEYS = ("schema_version", "components")
COMPONENT_KEYS = ("name", "formula", "MW", "Tc", "Pc", "omega")
# Property keys inside the parameter objects
PARAMETER_KEYS = ("value", "units", "source_key")
ANTOINE_KEYS = ("A", "B", "C", "Tmin_K", "Tmax_K", "units", "source_key")
