import shutil
import os


class RegressionTest(object):
    def __init__(self):
        self.path = None
        self.test_path = None
        self.golden_data_path = None

    def ini_path(self):
        try:
            if self.test_path:
                shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def clean_path(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
