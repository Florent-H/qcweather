import pytest
import os
from pathlib import Path
from qcweather import Weather
"""
~~~ PYTEST FIXTURES ~~~
"""

@pytest.fixture()
def epw_path():
    from pathlib import Path
    # CAN-QC - Montreal YUL 716270 - ISD 2015
    yield "datafiles\\epw_files\\tampered.epw"
    # yield "datafiles\\epw_files\\CAN_QC_Montreal-McTavish.716120_CWEC2016.epw"
    # yield "datafiles\\epw_files\\CAN-QC - Montreal YUL 716270 - ISD 2015.epw"

    # yield Path("datafiles/epw_files/tampered.epw")
    # yield Path("datafiles/epw_files/CAN_QC_Montreal-McTavish.716120_CWEC2016.epw")
    # yield Path("datafiles/epw_files/CAN-QC - Montreal YUL 716270 - ISD 2015.epw")

@pytest.fixture()
def drive_letter():
    yield "E"

@pytest.fixture()
def ashrae_path():
    yield "datafiles\\ashrae\\2017DesignConditions_s.xlsx"

@pytest.fixture()
def output_dir():
    yield "datafiles"

@pytest.fixture()
def config(output_dir):
    import shutil
    if os.path.isdir(Path(output_dir) / "raw_tbl"):
        shutil.rmtree(Path(output_dir) / "raw_tbl")
    (Path(output_dir) / "raw_tbl").mkdir()

class TestsQualityControl:
    """A class that encapsulates the different tests for the functions found in
    the quality_control.py module
    """

    """
    ~~~ FUNCTION TESTS ~~~
    """

    def test_run_quality_control(self, config, capsys, epw_path, drive_letter, ashrae_path, output_dir):

        with capsys.disabled():
            weather = Weather.get_weather(epw_path, drive_letter, ashrae_path, output_dir)
            weather.run_quality_control()

    def test_get_month_rad_graphs(self, config, capsys):

        with capsys.disabled():
            from qcweather import get_month_rad_graphs
            get_month_rad_graphs(drive_letter, ashrae_path, output_dir)

    def test_get_month_rad_extreme_graphs(self, config, capsys, drive_letter, ashrae_path, output_dir):
        with capsys.disabled():
            from qcweather import get_month_rad_extreme_graphs
            get_month_rad_extreme_graphs(drive_letter, ashrae_path, output_dir)

    def test_get_graph(self, config, capsys, epw_path):
        with capsys.disabled():
            weather = Weather.get_weather(epw_path, drive_letter, ashrae_path, output_dir)
            weather.get_graph("hourly", "dew_point")

