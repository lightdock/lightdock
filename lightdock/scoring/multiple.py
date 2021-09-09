import os


class ScoringConfiguration(object):
    """Manages scoring configuration files"""

    @staticmethod
    def parse_file(file_name):
        functions = {}
        with open(file_name) as input_file:
            lines = input_file.readlines()
            for line in lines:
                try:
                    line = line.rstrip(os.linesep)
                    if not line.startswith("#"):
                        fields = line.split()
                        function_name = fields[0]
                        weight = float(fields[1])
                        functions[function_name] = weight
                except ValueError:
                    pass
        return functions
