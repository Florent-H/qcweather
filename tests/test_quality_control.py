import pytest
import os
import shutil
from pathlib import Path
"""
~~~ PYTEST FIXTURES ~~~
"""

@pytest.fixture()
def epw_path():
    from pathlib import Path
    # CAN-QC - Montreal YUL 716270 - ISD 2015
    # yield Path("datafiles/epw_files/tampered.epw")
    yield Path("datafiles/epw_files/CAN_QC_Montreal-McTavish.716120_CWEC2016.epw")
    # yield Path("datafiles/epw_files/CAN-QC - Montreal YUL 716270 - ISD 2015.epw")

@pytest.fixture()
def config():
    import shutil
    if os.path.isdir(Path("datafiles/raw_tbl")):
        shutil.rmtree(Path("datafiles/raw_tbl"))
    Path("datafiles/raw_tbl").mkdir()

class TestsQualityControl:
    """A class that encapsulates the different tests for the functions found in
    the quality_control.py module
    """

    """
    ~~~ FUNCTION TESTS ~~~
    """

    def test_run_quality_control(self, config, capsys, epw_path):

        with capsys.disabled():
            from qcweather import run_quality_control
            run_quality_control(epw_path)

    def test_get_month_rad_graphs(self, config, capsys):

        with capsys.disabled():
            from qcweather import get_month_rad_graphs
            get_month_rad_graphs()