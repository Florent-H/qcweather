# qcweather

qcweather is a Python package for vetting Energyplus .epw weather files. It runs physics- and statistics-based quality control procedures to help check that the inputted weather file is reliable.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install qcweather.

```bash
pip install qcweather
```

NOTE: qcweather is not published yet, so it is not possible to install it via pip. The above is a placeholder.

## Usage

First import the ```Weather``` class from ```qcweather```.

```python
from qcweather import Weather
```

Then create an object from the ```Weather``` class, which requires the following arguments:

1. The path for the .epw weather file; 
2. The drive letter where ASHRAE's Weather Data Viewer (2017) iso file is mounted (not provided);
3. The path for ASHRAE's 2017 Design Conditions Excel file (not provided); and
4. The output directory for the intermediate and final results files.

```python
weather_path = "datafiles\\epw_files\\CAN-QC - Montreal YUL 716270 - 2015.epw"  # Windows path (Mac, Linux paths possible too)
drive_letter = "E"
ashrae_path = "datafiles\\ashrae\\2017DesignConditions_s.xlsx"
output_dir = "datafiles\\results"

weather = Weather.get_weather(weather_path, drive_letter, ashrae_path, output_dir)
```

Finally, run the quality control procedures.

```python
weather.run_quality_control()
```
A csv file with flagged hours, days, and months is created in the ```output_dir``` along with graphs showing these flagged data points in context with neighbouring data.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)