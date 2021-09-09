import os
from pathlib import Path

# Set global path variables
test_path = Path(__file__).absolute().parent
configuration_path = test_path / "etc"

os.environ["LIGHTDOCK_CONF_PATH"] = str(configuration_path)
