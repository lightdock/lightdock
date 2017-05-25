import shutil
import os


class RegressionTest(object):
    def __init__(self):
        self.path = ''
        self.test_path = ''
        self.golden_data_path = ''

    def ini_test_path(self):
        try:
            if self.test_path:
                shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def clean_test_path(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
