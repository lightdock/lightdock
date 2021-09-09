"""PyDock_Table: processes table-like files"""

import os


class Table(object):
    def __init__(self, arrays=None, header=None):
        """Create a new Table instance"""
        self.__numCols = 0
        self.__numRows = 0
        self.__content = {}
        self.__header = []

        if arrays is not None and len(arrays) >= 1:

            # Check for homogeneous array length
            length = len(arrays[0])
            for array in arrays[1:]:
                if len(array) != length:
                    raise TableError("Heterogeneous array length")

            if header is not None:
                if len(header) != len(arrays):
                    raise TableError("Header length different from array size")

            else:
                # Build header
                header = ["Col" + str(i) for i in range(len(arrays))]

            self.__numCols = len(arrays)
            self.__numRows = length
            self.__content = {}
            self.__header = header

            for headerName, array in zip(header, arrays):
                self.__content[headerName] = array

    def is_null(self):
        """Returns true if this table has no info"""
        return self.__numRows == 0 and self.__numCols == 0

    def is_empty(self):
        """Returns true if this table is empty (has column names, but not info"""
        return self.__numRows == 0 and self.__numCols != 0

    def get_column_names(self):
        """Get header"""
        return self.__header

    def get_content(self):
        """Get content"""
        return self.__content

    def get_num_rows(self):
        """Return number of rows"""
        return self.__numRows

    def get_num_cols(self):
        """Return number of columns"""
        return self.__numCols

    def clone(self):
        """Returns a copy of this table"""
        new_table = Table()
        new_table.__content = self.__content.copy()
        new_table.__header = self.__header[:]
        new_table.__numCols = self.__numCols
        new_table.__numRows = self.__numRows
        return new_table

    def __str__(self):
        return self._to_string()

    def show(self, extra_info=True):
        """Show a Table"""
        if extra_info:
            print("Table: %sx%s" % (self.__numCols, self.__numRows))
        print(self.__str__())

    def _to_string(self):
        """Show a Table instance"""
        if self.__numCols == 0:
            return "Empty Table"
        table_out = ""
        line1 = ""
        line2 = ""
        for key in self.__header:
            line1 += "%12s" % key
            line2 += "------------"

        table_out += line1 + os.linesep + line2 + os.linesep

        for i in range(self.__numRows):
            line = ""
            for key in self.__header:
                if type(self.__content[key][i]) == int:
                    line += "%12i" % self.__content[key][i]
                elif type(self.__content[key][i]) == float:
                    line += "%12.3f" % self.__content[key][i]
                else:
                    line += "%12s" % self.__content[key][i]

            table_out += line + os.linesep

        return table_out

    @staticmethod
    def read(table_file):
        """Read a Table from a given file"""
        try:
            f = open(table_file)
            tab = f.readlines()
            f.close()
        except IOError as e:
            raise TableError(str(e))

        if tab[0].split()[0] == "#>T":  # for icm table format
            tab[0] = tab[1].replace("#>", " ").replace("-", " ")
            tab[1] = "-----------"

        arrays = []
        header = tab[0].split()
        for i in range(len(header)):
            arrays.append([])

        # in case header has two lines (first for columns, second just a separator)
        if len(tab[1].split()) == len(tab[0].split()):
            i_first = 1
        else:
            i_first = 2

        for i in range(i_first, len(tab[i_first:]) + i_first):
            line = tab[i].split()
            for k in range(len(header)):
                try:
                    arrays[k].append(int(line[k]))
                except:
                    try:
                        arrays[k].append(float(line[k]))
                    except:
                        try:
                            arrays[k].append(str(line[k]))
                        except:
                            arrays[k].append(" ")

        return Table(arrays, header)

    def write(self, table_file, table_format=None):
        """Write a Table to a given file"""
        try:
            f = open(table_file, "w")
            if table_format == "icm":
                f.write("#>T TABLE" + os.linesep)
                f.write("#>-")
                for key in self.__header:
                    col = "%12s" % key
                    f.write(col.replace(" ", "-"))
                f.write(os.linesep)
            else:
                line1 = ""
                line2 = ""
                for key in self.__header:
                    line1 += "%12s" % key
                    line2 += "------------"
                f.write(line1 + os.linesep)
                f.write(line2 + os.linesep)

            for i in range(self.__numRows):
                for key in self.__header:
                    content_type = type(self.__content[key][i])
                    if content_type == int:
                        f.write("%12i" % self.__content[key][i])
                    elif content_type == float:
                        f.write("%12.3f" % self.__content[key][i])
                    elif content_type.__name__ == "float64":
                        f.write("%12.3f" % self.__content[key][i])
                    else:
                        f.write("%12s" % self.__content[key][i])
                f.write(os.linesep)

            f.close()

        except Exception as e:
            raise TableError(str(e))

    def sort(self, column_to_sort, reverse=False):
        """Sort a table by columnToSort"""
        if type(column_to_sort) is not str or column_to_sort not in self.__header:
            raise TableError("Column does not exist: %s" % column_to_sort)

        indices = list(range(self.__numRows))
        indices.sort(key=self.__content[column_to_sort].__getitem__, reverse=reverse)

        new_content = dict([(i, list()) for i in self.__header])

        for key in self.__header:
            for i in range(self.__numRows):
                new_content[key].append(self.__content[key][indices[i]])

        self.__content = new_content

    def get_tab_by_columns(self, col_names):
        """Creates a new Table with selected columns"""
        content = []
        try:
            for key in col_names:
                if key not in self.__header:
                    raise TableError("Column name not found: %s" % key)
                content.append(self.__content[key])
        except Exception as e:
            raise TableError(str(e))

        return Table(content, col_names)

    def get_tab_by_rows(self, row_from, row_to):
        """Creates a new Table with selected rows (first row is 0, last is N-1)"""
        if (
            row_from < 0
            or row_from > self.__numRows
            or row_to < 0
            or row_to > self.__numRows
        ):
            raise TableError(
                "Error in selection: index must be in [0:%s]" % self.__numRows
            )

        content = []
        try:
            for key in self.__header:
                content.append(self.__content[key][row_from : row_to + 1])
        except Exception as e:
            raise TableError(str(e))

        return Table(content, self.__header)

    def delete_columns(self, col_names):
        """Delete given columns from Table"""
        if type(col_names) is str:
            col_names = [col_names]

        if not set(col_names).issubset(set(self.__header)):
            raise TableError("One or more columns does not exist")

        for colName in col_names:
            self.__header.remove(colName)
            del self.__content[colName]
            self.__numCols -= 1

    def append_columns_from_table(self, orig_table):
        """Append columns from origTable"""
        if set(orig_table.__header).issubset(set(self.__header)):
            raise TableError("One or more columns override")

        if self.__numRows != orig_table.__numRows:
            raise TableError("Different number of rows")

        self.__header.extend(orig_table.__header)
        self.__numCols += orig_table.__numCols
        self.__content.update(orig_table.__content)

    def append_rows_from_table(self, orig_table):
        """Append to this table origTable rows"""
        if not isinstance(orig_table, Table):
            raise TableError("Argument not a Table")

        if not set(self.__header).issubset(set(orig_table.__header)):
            raise TableError("Columns does not correspond")

        for key in self.__header:
            self.__content[key].extend(orig_table.__content[key])

        self.__numRows += orig_table.__numRows

    def append_array(self, col_array, col_name=None):
        """Append a given array to a Table"""
        if not self.is_null() and len(col_array) != self.__numRows:
            raise TableError("Array size and number of rows does not correspond")

        if col_name in self.__header:
            raise TableError("Column name already exist")

        if col_name is None:
            col_name = "Col" + str(self.__numCols)

        self.__header.append(col_name)
        self.__numCols += 1
        if self.is_empty:
            self.__numRows = len(col_array)
        self.__content[col_name] = col_array[:]

    def get_column(self, col_name):
        """Get a column identified by colName"""
        if col_name not in self.__header:
            raise TableError("colName is not a valid column name")

        return self.__content[col_name]

    def __getitem__(self, key):
        return self.get_column(key)

    def rename(self, old_header, new_header):
        """Rename Table columns"""
        try:
            i = self.__header.index(old_header)
        except ValueError:
            raise TableError("%s is not a valid column name" % old_header)

        if self.__header.count(new_header) > 0:
            raise TableError("%s already exists" % new_header)

        self.__header[i] = new_header
        self.__content[new_header] = self.__content[old_header]
        del self.__content[old_header]

    def add_index_column(self, name="RANK"):
        """Adds a new column with numeric indexes"""
        indexes = [i for i in range(1, self.get_num_rows() + 1)]
        self.append_array(indexes, name)


class TableError(Exception):
    """Table Exception class"""

    def __init__(self, value):
        self.parameter = value

    def __str__(self):
        return self.parameter
