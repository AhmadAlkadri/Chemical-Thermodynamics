# Component database

The `flash` module looks for a consolidated HDF5 store at `database/database.h5` by default. The `DEFAULT_DATABASE_PATH` constant exported by `flash` points to this location so that client code and notebooks do not need to hard-code paths.

The `database.h5` file is generated from the tab-delimited text sources in this directory (for example `organics.txt` and `inorganics.txt`). The `sort_database.py` helper can be used to rebuild the HDF5 store from those raw tables.
