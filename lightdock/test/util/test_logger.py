"""Tests for Logger module"""
import shutil
import os
import filecmp
import sys
from pathlib import Path
from io import StringIO
from lightdock.util.logger import Logger


class TestLogger:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_logger"
        self.golden_data_path = self.path / "golden_data"

    def setup(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def teardown(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass

    def test_create_default_logger(self):
        log = Logger("test")

        assert Logger.INFO == log.get_level()
        assert log.get_tag() == "test"

    def test_set_level(self):
        log = Logger("test")

        log.set_level(Logger.DEBUG)
        assert Logger.DEBUG == log.get_level()

        log.set_level("a")
        assert Logger.DEBUG == log.get_level()

        log.set_level("-1")
        assert Logger.DEBUG == log.get_level()

    def test_create_logger_with_file(self):
        log = Logger("test", file_name=self.test_path / "test_file.log")
        log.info("Testing")

        assert (self.test_path / "test_file.log").exists()

    def test_output_debug(self):
        log = Logger("test", file_name=self.test_path / "test_file.log")
        log.debug("Testing DEBUG")
        log.set_level(Logger.DEBUG)
        log.debug("Testing DEBUG")

        assert filecmp.cmp(
            self.golden_data_path / "output_debug.log", self.test_path / "test_file.log"
        )

    def test_output_warning(self):
        log = Logger("test", file_name=self.test_path / "test_file.log")
        log.warning("Testing WARNING")

        assert filecmp.cmp(
            self.golden_data_path / "output_warning.log",
            self.test_path / "test_file.log",
        )

    def test_output_error(self):
        log = Logger("test", file_name=self.test_path / "test_file.log")
        log.error("Testing ERROR")

        assert filecmp.cmp(
            self.golden_data_path / "output_error.log", self.test_path / "test_file.log"
        )

    def test_output_std_out(self):
        output = StringIO()
        log = Logger("test")
        saved_stdout = sys.stdout
        sys.stdout = output
        log.info("Testing")
        sys.stdout = saved_stdout

        assert output.getvalue().rstrip() == "[test] INFO: Testing"
