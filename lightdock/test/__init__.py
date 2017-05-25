import os

# Set global path variables
test_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__)))
configuration_path = "%s%s%s%s" % (test_path, os.sep, 'etc', os.sep)

os.environ['LIGHTDOCK_CONF_PATH'] = configuration_path
