# Copyright 2006-2025 Mark Diekhans

# FIXME: danger of dump, etc, methods conflicting with columns.  maybe
# a better convention to avoid collisions or make these functions rather
# than methods
# FIXME: need accessor functions for columns
# FIXME: need way to get raw row with Nones for sql
# FIXME: this could actually be a dict-like object since py3

def tsvRowToDict(row):
    """convert a TSV row to a dict"""
    return {col: getattr(row, col) for col in row._columnSpecs_.columns}


class TsvRow:
    "Row of a TSV where columns are fields."
    # n.b.: doesn't inherit from list, as this results in columns in two
    # places when they are stored as fields

    def __init__(self, reader, row):
        # disconnect from reader to allow object to be pickled
        self._columnSpecs_ = reader.columnSpecs
        for i in range(len(self._columnSpecs_.columns)):
            setattr(self, self._columnSpecs_.columns[i], row[i])

    def __getitem__(self, key):
        "access a column by string key or numeric index"
        if isinstance(key, int):
            return getattr(self, self._columnSpecs_.columns[key])
        else:
            return getattr(self, key)

    def __setitem__(self, key, val):
        "set a column by string key or numeric index"
        if isinstance(key, int):
            setattr(self, self._columnSpecs.columns[key], val)
        else:
            setattr(self, key, val)

    def get(self, key, default=None):
        """Return the value for string key or numeric index, if key is in the
        dictionary, else default."""
        if isinstance(key, int):
            return getattr(self, self._columnSpecs_.columns[key]) if key < len(self._columnSpecs_.columns) else default
        else:
            return getattr(self, key, default)

    def __len__(self):
        return len(self._columnSpecs_.columns)

    def __iter__(self):
        for col in self._columnSpecs_.columns:
            yield getattr(self, col)

    def __contains__(self, key):
        return hasattr(self, key)

    def __str__(self):
        return "\t".join(self.getRow())

    def __repr__(self):
        return str(self)

    def getRow(self):
        return self._columnSpecs_.formatRow(self)

    def getColumns(self, colNames):
        """get a subset of the columns in the row as a list"""
        subRow = []
        for col in colNames:
            subRow.append(self[col])
        return subRow

    def getHeader(self):
        return "\t".join(self._columnSpecs_.extColumns)

    def writeHeader(self, fh):
        fh.write(self.getHeader())
        fh.write("\n")

    def write(self, fh):
        fh.write(str(self))
        fh.write("\n")

    def dump(self, fh):
        i = 0
        for col in self:
            if i > 0:
                fh.write("\t")
            fh.write(self._columnSpecs_.columns[i])
            fh.write(": ")
            fh.write(str(col))
            i += 1
        fh.write("\n")
