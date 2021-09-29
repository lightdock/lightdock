"""Logging utilities"""

import os


class Logger(object):

    (ERROR, WARNING, INFO, DEBUG) = (0, 1, 2, 4)

    def __init__(self, tag, file_name="", level=INFO):
        """Creates a Logger object.

        If file_name is set, output will be printed to a file called file_name.
        Default level is set to 'info', debug messages will be ignored.

        """
        if file_name:
            self._output = open(file_name, "a")
            self._write = self._write_to_file
        else:
            self._write = self._write_to_std_out

        self._level = level
        self._tag = tag

    def set_level(self, level):
        """Sets logging level"""
        try:
            level = int(level)
            if level >= Logger.ERROR and level <= Logger.DEBUG:
                self._level = level
        except:
            pass

    def get_level(self):
        """Gets the current level of logging"""
        return self._level

    def get_tag(self):
        """Gets the logging name"""
        return self._tag

    def debug(self, message=""):
        """Prints to selected output a debug message"""
        if self._level >= self.DEBUG:
            self._write(message, "DEBUG")

    def info(self, message=""):
        """Prints to selected output an informative message"""
        if self._level >= self.INFO:
            self._write(message, "INFO")

    def warning(self, message=""):
        """Prints to selected output a warning message"""
        if self._level >= self.WARNING:
            self._write(message, "WARNING")

    def error(self, message=""):
        """Prints to selected output an error message"""
        if self._level >= self.ERROR:
            self._write(message, "ERROR")

    def _write_to_file(self, message, level):
        """Outputs log info to the selected file"""
        out_message = "[%s] %s: %s" % (self._tag, level, message)
        self._output.write(out_message + os.linesep)
        self._output.flush()

    def _write_to_std_out(self, message, level):
        """Outputs log info to the standard output"""
        out_message = "[%s] %s: %s" % (self._tag, level, message)
        print(out_message)

    def __del__(self):
        """Frees file descriptor"""
        try:
            self._output.close()
        except:
            pass


class LoggingManager(object):
    """Logger facility"""

    _loggers = {}

    @staticmethod
    def get_logger(tag, file_name="", level=Logger.INFO):
        """Gets logger object identified by tag if exists or creates a new one"""
        try:
            return LoggingManager._loggers[tag]
        except KeyError:
            try:
                level = int(os.environ["LIGHTDOCK_LOGGING"])
            except KeyError:
                pass
            log = Logger(tag, file_name, level)
            LoggingManager._loggers[tag] = log
            return log
