import pytest


class TestsQualityAssurance:
    """A class that encapsulates the different tests for the functions found in
    the sensitivity_analysis.py module
    """

    """
    ~~~ PYTEST FIXTURES ~~~
    """

    @pytest.fixture()
    def epw_path(self):
        from pathlib import Path
        yield Path("datafiles/epw_files/CAN_QC_Montreal-McTavish.716120_CWEC2016.epw")

    """
    ~~~ FUNCTION TESTS ~~~
    """

    def test_run_quality_assurance(self, capsys, epw_path):

        with capsys.disabled():
            from qcweather import run_quality_assurance
            run_quality_assurance(epw_path)

    def test_do_multi_lin_reg(self, capsys, a):
        with capsys.disabled():
            from qcweather import y

    def test_write_csv_table(self, a):
        with capsys.disabled():
            from qcweather import z