"""Tests for Logger module"""
import shutil
import os
import filecmp
import sys

from lightdock.util.logger import Logger
from io import StringIO


class TestLogger:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch/'
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
    
    def test_create_default_logger(self):
        log = Logger('test')
        
        assert Logger.INFO == log.get_level()
        assert "test" == log.get_tag()
        
    def test_set_level(self):
        log = Logger('test')
        
        log.set_level(Logger.DEBUG)
        assert Logger.DEBUG == log.get_level()
        
        log.set_level('a')
        assert Logger.DEBUG == log.get_level()
        
        log.set_level('-1')
        assert Logger.DEBUG == log.get_level()
        
    def test_create_logger_with_file(self):
        log = Logger('test', file_name=self.test_path+'test_file.log')
        log.info('Testing')
        
        assert os.path.exists(self.test_path+'test_file.log')
        
    def test_output_debug(self):
        log = Logger('test', file_name=self.test_path+'test_file.log')
        log.debug('Testing DEBUG')
        log.set_level(Logger.DEBUG)
        log.debug('Testing DEBUG')
        
        assert filecmp.cmp(self.golden_data_path+'output_debug.log', self.test_path+'test_file.log')
        
    def test_output_warning(self):
        log = Logger('test', file_name=self.test_path+'test_file.log')
        log.warning('Testing WARNING')
        
        assert filecmp.cmp(self.golden_data_path+'output_warning.log', self.test_path+'test_file.log')
        
    def test_output_error(self):
        log = Logger('test', file_name=self.test_path+'test_file.log')
        log.error('Testing ERROR')
        
        assert filecmp.cmp(self.golden_data_path+'output_error.log', self.test_path+'test_file.log')
        
    def test_output_std_out(self):
        output = StringIO()        
        log = Logger('test')
        saved_stdout = sys.stdout
        sys.stdout = output
        log.info("Testing")        
        sys.stdout = saved_stdout

        assert "[test] INFO: Testing" == output.getvalue().rstrip()
