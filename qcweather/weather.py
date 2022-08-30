import csv
import math
import os
import re
import subprocess
import pandas as pd
import numpy as np
import calendar
import datetime
from pathlib import Path
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
from scipy.stats import norm


class Weather(object):
    def __init__(self, weather_path, drive_letter, ashrae_path, output_dir):
        """
        @type weather_path: basestring, pathlib Path
        @param weather_path: the file path (either Windows or Mac/Linux) for the .epw weather file
        @type drive_letter: basestring
        @param drive_letter: the drive letter (without the colon) where the ASHRAE Weather Viewer (2017) iso file is
        mounted
        @type ashrae_path: basestring, pathlib Path
        @param ashrae_path: the file path (either Windows or Mac/Linux) for the ASHRAE Design Conditions (2017)
        Excel file
        @type output_dir: basestring, pathlib Path
        @param output_dir: the output folder/directory where intermediate and final results files will be saved
        """
        from pathlib import Path

        self.weather_path = Path(weather_path)
        self.drive_letter = drive_letter
        self.ashrae_path = Path(ashrae_path)
        self.output_dir = Path(output_dir)
        self.wmo = None  # World Meteorological Organization (WMO) weather station (type: int)
        self.folder = None  # folder for the raw cumulative distribution function (CDF) tables (type: pathlib Path)
        self.location = None  # tuple of (latitude, longitude) for weather station in radians (type: tuple of floats)
        self.time_zone = None  # offset from Coordinated Universal Time (UTC) in hours; negative to the west
        # (type: float)
        self.meteo_vars = None  # table of meteorological variables found in .epw weather file (type: pandas DataFrame)
        self.solar_angles = None  # four different solar angles or solar attributes ("declination", "hour_angle",
        # "zenith", "day") (type: dictionary)
        self.historical = None  # monthly probability distribution functions as well as monthly averages and standard
        # deviations for some meteorological variables (type: dictionary)

    def get_meteo_params(self):
        """
        Assign meteorological parameters (World Meteorological Organization (WMO) station number, latitide + longitude,
        and time zone) and meteorological variables (year, month, dry bulb temperature, pressure,
        etc.)) from .epw weather file to attributes of Weather() object
        """
        # read in header information of .epw file as well as first datetime index
        if self.weather_path.suffix == ".epw":
            # open energy plus weather file
            with open(self.weather_path, "r") as fp:
                reader = csv.reader(fp)
                # get relevant rows
                rel_rows = [row for idx, row in enumerate(reader) if idx in (0, 8)]

                # read in header information
                header_line = rel_rows[0]
                self.wmo = int(header_line[5])
                # (latitude, longitude) in degrees -> convert to radians
                # using convention that to the west is positive
                self.location = (
                    math.radians(float(header_line[6])),
                    math.radians(-float(header_line[7])),
                )
                self.time_zone = -float(header_line[8])

                # read in first datetime index
                dti_line = rel_rows[1]
                first_date = datetime.datetime(
                    int(dti_line[0]),
                    int(dti_line[1]),
                    int(dti_line[2]),
                    int(dti_line[3]),
                )

            # all header names
            all_cols = [
                "year",
                "month",
                "day",
                "hour",
                "minute",
                "source_uncertainty",
                "dry_bulb",
                "dew_point",
                "rel_hum",
                "pressure",
                "ext_hor_rad",
                "ext_dir_norm",
                "hor_inf_rad_int",
                "glob_hor_rad",
                "dir_norm_rad",
                "diff_hor_rad",
                "glob_hor_ill",
                "dir_norm_ill",
                "diff_hor_ill",
                "zen_lum",
                "wind_dir",
                "wind_speed",
                "tot_sky_cov",
                "op_sky_cov",
                "visibility",
                "ceil_height",
                "pres_wth_obs",
                "pres_wth_cod",
                "prec_water",
                "aero_opt_dep",
                "snow_depth",
                "days_lst_snow",
                "albedo",
                "liq_prec_dep",
                "liq_prec_qty",
            ]

            # read epw as csv to create pandas DataFrame for meteo_vars attribute
            meteo_vars = pd.read_csv(
                self.weather_path,
                skiprows=8,
                header=None,
                names=all_cols,
                usecols=all_cols[6:],
                index_col=False,
                dtype=np.float64,
            )

            # get datetime index for dataframe using first datetime and number of rows in dataframe
            meteo_vars.index = pd.date_range(
                first_date, periods=len(meteo_vars.index), freq="H"
            )

            # add columns for quality control procedure flags
            meteo_vars["hourly_flags"] = [[] for _ in meteo_vars.index]
            meteo_vars["step_flags"] = [[] for _ in meteo_vars.index]
            meteo_vars["daily_flags"] = [[] for _ in meteo_vars.index]
            meteo_vars["monthly_flags"] = [[] for _ in meteo_vars.index]

            # add flag for missing data points
            rows, cols = np.where(meteo_vars == 9999)
            # get integer position of hourly_flags column
            hf_col = list(meteo_vars.columns).index("hourly_flags")
            for row, col in zip(rows, cols):
                meteo_vars.iloc[row, hf_col].append(
                    f"missing data for {list(meteo_vars.columns)[col]}"
                )

        print("done reading epw file")
        self.meteo_vars = meteo_vars

    def get_solar_angles(self):
        """
        Assign solar angles and solar attributes (declination angle, hour angle, day/night boolean, and zenith angle)
        to solar_angles attribute of Weather() object
        """
        # get location of weather station in radians
        latitude = self.location[0]
        longitude = self.location[1]

        # solar_angles = pd.DataFrame(None, columns=["declination", "hour_angle", "zenith", "altitude", "azimuth", "day"],
        #     index=self.meteo_vars.index,
        # )
        solar_angles = pd.DataFrame(
            None,
            columns=["declination", "hour_angle", "zenith", "day"],
            index=self.meteo_vars.index,
        )

        for time, _ in self.meteo_vars.iterrows():
            # fractional day of the year
            n = time.dayofyear - 0.5 + (time.hour - 0.5) / 24
            # equation of time
            B = math.radians((n - 1) * 360 / 365)
            E = 2.2918 * (
                0.0075
                + 0.1868 * math.cos(B)
                - 3.2077 * math.sin(B)
                - 1.4615 * math.cos(2 * B)
                - 4.089 * math.sin(2 * B)
            )

            # solar declination angle
            declin = (
                0.006918
                - 0.399912 * math.cos(B)
                + 0.070257 * math.sin(B)
                - 0.006758 * math.cos(2 * B)
                + 0.000907 * math.sin(2 * B)
                - 0.002697 * math.cos(3 * B)
                + 0.00148 * math.sin(3 * B)
            )
            solar_angles["declination"][time] = declin

            # solar time
            t_s = (
                time.hour
                - 0.5
                + (4 * (self.time_zone * 15 - math.degrees(longitude)) + E) / 60
            ) % 24
            # solar hour angle
            hour_angle = math.radians((t_s - 12) * 15)
            solar_angles["hour_angle"][time] = hour_angle

            # calculate sunrise and sunset angles
            sunset = math.acos(-math.sin(latitude) * math.tan(declin))
            sunrise = -sunset
            if hour_angle < sunrise:
                solar_angles["day"][time] = False

            elif hour_angle > sunset:
                solar_angles["day"][time] = False

            else:
                solar_angles["day"][time] = True

            # solar zenith angle
            zenith = math.acos(
                math.cos(latitude) * math.cos(declin) * math.cos(hour_angle)
                + math.sin(latitude) * math.sin(declin)
            )

            solar_angles["zenith"][time] = zenith

            # altitude = math.asin(math.cos(zenith))
            # solar_angles["altitude"][time] = altitude
            #
            # # todo: find equivalent of sign() in python (not math.copysign)
            #
            # azimuth = math.copysign(
            #     abs(
            #         math.acos(
            #             (math.cos(zenith) * math.sin(latitude) - math.sin(declin))
            #             / (math.sin(zenith) * math.cos(latitude)),
            #         )
            #     ),
            #     hour_angle,
            # )
            # solar_angles["azimuth"][time] = azimuth

        print("done calculating solar angles")
        self.solar_angles = solar_angles

    def get_pdf(self, meteo_var):
        """
        get monthly historical probability distribution functions (PDFs) for either dry bulb temperature, dew point
        temperature, relative humidity, wind direction, wind speed, and dry bulb temperature crossed with time of day
        as dictionaries of list of 12 lists (months) or list of 12 lists (months) of 24 lists (hours) for crossed PDF
        @type meteo_var: basestring
        @param meteo_var: the meteorological variable for the PDF creation
        @return: probability distribution function (PDF) for the passed meteorological variable
        """
        # get the meteorological variable ID found in the cumulative distribution text files, or warn the user if the
        # wrong meteo_var string is passed as an argument
        if meteo_var == "dry_bulb":
            meteo_id = "DB"

        elif meteo_var == "dew_point":
            meteo_id = "DP"

        elif meteo_var == "rel_hum":
            meteo_id = "DBDP"

        elif meteo_var == "wind_dir":
            meteo_id = "WD"

        elif meteo_var == "wind_speed":
            meteo_id = "WS"

        elif meteo_var == "dry_bulb_time":
            meteo_id = "DBTD"
        else:
            raise ValueError(
                "meteo_var must be either 'dry_bulb', 'dew_point', 'rel_hum', 'wind_dir', 'wind_speed', or 'dry_bulb_time'"
            )

        # create probability distribution function dictionaries with either list of 12 lists (months) for single
        # meteorological variable cumulative distribution function (CDF) text tables or list of 12 lists (months) of
        # 24 lists (hours) for crossed CDF text tables
        if meteo_id == "DBTD":
            pdf = {
                "bin": [[[] for i in range(24)] for i in range(12)],
                "freq": [[[] for i in range(24)] for i in range(12)],
                "cdf": [[[] for i in range(24)] for i in range(12)],
                "pdf": [[[] for i in range(24)] for i in range(12)],
            }

        else:
            pdf = {
                "bin": [[] for i in range(12)],
                "freq": [[] for i in range(12)],
                "cdf": [[] for i in range(12)],
                "pdf": [[] for i in range(12)],
            }

        # for each month of the year, read in the data from the cumulative distribution function (CDF) text tables,
        # which were generated with ASHRAE Weather Viewer's tblxpand executable
        for m in range(12):
            # add a zero in front of the month string if the month number is less than 10 to match the text file path
            # names
            if m < 9:
                month = f"0{m + 1}"
            else:
                month = f"{m + 1}"
            month_path = self.folder / f"{self.wmo}_{meteo_id}_{month}.txt"

            # open CDF text file for the month
            with open(month_path, "r") as fp:
                # skip the five-line or four-line header
                if meteo_id in ("DBDP", "DBTD"):
                    skip = 4
                else:
                    skip = 5

                for n in range(skip):
                    next(fp)

                # for each line in the CDF text table, read in the bins, frequencies, and probabilities
                for j, line in enumerate(fp):
                    # to get the relative humidity PDF, the crossed dry bulb and dew point temperature CDF
                    # must be used; first read the dry bulb and dew point temperatures as lists of floats
                    if meteo_id == "DBDP":
                        if j == 0:
                            dew_point = [float(s) for s in line.split()[1:]]
                            dry_bulb = []
                            freq = []

                        else:
                            dry_bulb.append(float(line.split()[0]))
                            freq.append([float(s) for s in line.split()[1:]])

                    # to get the dry bulb temperature crossed with the time of day, read in the time of day and dry bulb
                    # temperature as lists of floats
                    elif meteo_id == "DBTD":
                        if j == 0:
                            time_day = [float(s) for s in line.split()[1:]]
                            dry_bulb_day = []
                            freq_day = [[] for s in line.split()[1:]]

                        else:
                            dry_bulb_day.append(float(line.split()[0]))
                            for s, _ in enumerate(freq_day):
                                freq_day[s].append(float(line.split()[s + 1]))

                    # for the other PDFs (dry bulb temperature, dew point temperature, wind direction, and wind speed,
                    # directly read in the bin, frequency, and probability
                    else:
                        pdf["bin"][m].append(float(line.split()[0]))
                        pdf["freq"][m].append(float(line.split()[1]))
                        pdf["cdf"][m].append(float(line.split()[2]))

            # use the lists of dry bulb and dew point temperature to calculate the relative humidity
            if meteo_id == "DBDP":
                bin_unsorted = []
                freq_unsorted = []
                run_freq_sum = [0]

                for x, db in enumerate(dry_bulb):
                    for y, dp in enumerate(dew_point):
                        bin_unsorted.append(
                            100
                            * (17.625 * dp / (243.04 + dp))
                            / (
                                17.625
                                * (db + 0.00000000000000000000001)
                                / (243.04 + db)
                            )
                        )
                        freq_unsorted.append(freq[x][y])
                        run_freq_sum.append(run_freq_sum[-1] + freq[x][y])

                # sort bins and corresponding frequencies
                pdf["bin"][m] = sorted(bin_unsorted)
                pdf["freq"][m] = [
                    f for _, f in sorted(zip(pdf["bin"][m], freq_unsorted))
                ]
                pdf["cdf"][m] = [
                    c / run_freq_sum[-1]
                    for _, c in sorted(zip(pdf["bin"][m], run_freq_sum))
                ]

            # if meteorological variable is dry bulb temperature crossed with time of day, calculate PDF for each month
            # and hour of the day
            if meteo_id == "DBTD":
                for h, _ in enumerate(time_day):
                    run_freq_sum = 0
                    for x, _ in enumerate(dry_bulb_day):
                        pdf["bin"][m][h].append(dry_bulb_day[x])
                        pdf["freq"][m][h].append(freq_day[h][x])
                        pdf["pdf"][m][h].append(freq_day[h][x] / sum(freq_day[h]))
                        run_freq_sum += freq_day[h][x]
                        pdf["cdf"][m][h].append(run_freq_sum / sum(freq_day[h]))

            # otherwise, calculate probability distribution function for each month
            else:
                total_obs = sum(pdf["freq"][m])
                pdf["pdf"][m] = [f / total_obs for f in pdf["freq"][m]]

        return pdf

    def get_historical(self):
        # create folder for collecting cumulative distribution function text tables
        self.folder = Path(self.output_dir) / f"raw_tbl/{self.wmo}"
        os.mkdir(self.folder)

        # run tblexpand to gather raw historical weather data
        run_string = (
            str(Path(f"{self.drive_letter}:/tblxpand/tblxpand.exe"))
            + f" {self.folder} "
            + str(Path(f"{self.drive_letter}:/data/{self.wmo}.wdv"))
            + " ALL 99 SI"
        )
        subprocess.run(run_string, cwd=Path.cwd(), capture_output=True)

        # get monthly probability distribution functions
        db_pdf = self.get_pdf("dry_bulb")
        dp_pdf = self.get_pdf("dew_point")
        rel_pdf = self.get_pdf("rel_hum")
        wd_pdf = self.get_pdf("wind_dir")
        ws_pdf = self.get_pdf("wind_speed")
        dbtd_pdf = self.get_pdf("dry_bulb_time")
        print("done getting monthly PDFs")

        # get monthly averages and standard deviations for irradiation and precipitation depth from ASHRAE Design
        # Conditions (2017) Excel file
        design_conditions = pd.read_excel(
            Path(self.ashrae_path),
            skiprows=4,
            header=None,
        )
        historical_avgs = pd.DataFrame(
            None,
            columns=[
                "clr_sky_dir_norm",
                "clr_sky_diff_hor",
                "dir_norm_taub",
                "diff_hor_taud",
                "all_sky_glob_avg",
                "all_sky_glob_std",
                "prec_dep_avg",
                "prec_dep_max",
                "prec_dep_min",
                "prec_dep_std",
            ],
            index=calendar.month_name[1:],
        )
        idx = design_conditions.index[design_conditions.iloc[:, 4] == self.wmo]
        historical_avgs["clr_sky_dir_norm"] = design_conditions.loc[
            idx, 533:544
        ].values.flatten()
        historical_avgs["clr_sky_dir_norm"] = design_conditions.loc[
            idx, 545:556
        ].values.flatten()
        historical_avgs["dir_norm_taub"] = design_conditions.loc[
            idx, 509:520
        ].values.flatten()
        historical_avgs["diff_hor_taud"] = design_conditions.loc[
            idx, 521:532
        ].values.flatten()
        historical_avgs["all_sky_glob_avg"] = design_conditions.loc[
            idx, 557:568
        ].values.flatten()
        historical_avgs["all_sky_glob_std"] = design_conditions.loc[
            idx, 569:580
        ].values.flatten()
        historical_avgs["prec_dep_avg"] = design_conditions.loc[
            idx, 206:217
        ].values.flatten()
        historical_avgs["prec_dep_min"] = design_conditions.loc[
            idx, 219:230
        ].values.flatten()
        historical_avgs["prec_dep_min"] = design_conditions.loc[
            idx, 232:243
        ].values.flatten()
        historical_avgs["prec_dep_std"] = design_conditions.loc[
            idx, 245:256
        ].values.flatten()
        print("done reading historical averages and std")
        self.historical = OrderedDict(
            (
                ("dry_bulb_pdf", db_pdf),
                ("dew_point_pdf", dp_pdf),
                ("rel_hum_pdf", rel_pdf),
                ("wind_dir_pdf", wd_pdf),
                ("wind_speed_pdf", ws_pdf),
                ("dry_bulb_time_pdf", dbtd_pdf),
            )
        )
        for col, avg in historical_avgs.iteritems():
            self.historical[col] = avg.values

    def run_quality_control(self):
        # run quality assurance on the hourly magnitudes of the meteorological variables
        self.hour_qc()

        # run quality assurance on the hourly steps of the meteorological variables
        self.step_qc()

        # run quality assurance on the daily profiles of the meteorological variables
        self.day_qc()

        # run quality assurance on the monthly magnitudes of the meteorological variables
        self.month_qc()

        # get graphs for different time spans and meteorological variables
        self.get_graph("hourly")
        self.get_graph("daily")
        self.get_graph("monthly")

        # write results csv file
        csv_path = Path(self.output_dir) / f"{self.weather_path.stem}_results.csv"
        self.meteo_vars.to_csv(csv_path)

        return self.meteo_vars

    def hour_qc(self):
        """
        check that meteorological variable values do not exceed hourly thresholds
        """

        for time, col in self.meteo_vars.iterrows():

            self.rad_hour_checks(time)
            self.db_hour_checks(time)
            self.dp_hour_check(time)
            self.rh_hour_check(time)
            self.ws_hour_check(time)

        print("done hourly value checks")

    def rad_hour_checks(self, time):

        # -----------------------
        # RADIATION / IRRADIANCE
        # -----------------------
        # (since radiation is in units of Wh/m^2 in a time series with hour time steps, both
        # are equivalent)

        # Solar constant
        G_sc = 1361

        # todo: add smaller than 0.03 * GHItoa (G_o)

        # if it is day and the solar zenith angle is less than 90 degrees
        if (
            self.solar_angles["day"][time]
            and self.solar_angles["zenith"][time] < math.pi / 2
        ):

            # calculate extraterrestrial solar radiation
            n = time.dayofyear - 0.5 + (time.hour - 0.5) / 24
            G_on = G_sc * (1 + 0.0334 * math.cos(math.radians((360 * (n - 3)) / 365)))
            G_o = G_on * math.cos(self.solar_angles["zenith"][time])

            # get month number
            month = time.month - 1

            # ---------------------------------------
            # global horizontal radiation/illuminance
            # ---------------------------------------

            # impossible physics check
            if self.meteo_vars["glob_hor_rad"][time] > min(
                1.2 * G_on,
                1.5 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 100,
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal irradiance is physically impossible [10.1191/0143624402bt038oa]"
                )

            # rare check
            if (
                self.meteo_vars["glob_hor_rad"][time]
                > 1.2 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 50
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal irradiance is rare [10.1191/0143624402bt038oa]"
                )

            # if it's not the first hour before sunrise or the final hour before sunset
            if (
                self.meteo_vars["glob_hor_rad"][time] < 0.03 * G_o
                and self.solar_angles["day"][time + datetime.timedelta(hours=-1)]
                and self.solar_angles["day"][time + datetime.timedelta(hours=+1)]
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal irradiance is rare [10.1191/0143624402bt038oa]"
                )

            # -----------------------------------
            # direct normal radiation/illuminance
            # -----------------------------------

            # impossible physics check
            if self.meteo_vars["dir_norm_rad"][time] > G_on:
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is physically impossible [BSRN Global Network recommended QC tests, V2.0]"
                )

            # rare check
            if (
                self.meteo_vars["dir_norm_rad"][time]
                > 0.95 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 0.2 + 10
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is rare [BSRN Global Network recommended QC tests, V2.0]"
                )

            # if it's not the first hour before sunrise or the final hour before sunset
            if (
                self.meteo_vars["dir_norm_rad"][time] == 0
                and self.solar_angles["day"][time + datetime.timedelta(hours=-1)]
                and self.solar_angles["day"][time + datetime.timedelta(hours=+1)]
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is rare [BSRN Global Network recommended QC tests, V2.0]"
                )

            # ----------------------------------------
            # diffuse horizontal radiation/illuminance
            # ----------------------------------------

            # impossible physics check
            if self.meteo_vars["diff_hor_rad"][time] > min(
                0.8 * G_on,
                0.95 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 50,
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is physically impossible [10.1016/S0038-092X(02)00121-4]"
                )

            # rare check
            if (
                self.meteo_vars["diff_hor_rad"][time]
                > 0.75 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 30
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is rare [10.1016/S0038-092X(02)00121-4]"
                )

            # if it's not the first hour before sunrise
            if (
                self.meteo_vars["diff_hor_rad"][time] < 0.03 * G_o
                and self.solar_angles["day"][time + datetime.timedelta(hours=-1)]
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is rare [10.1016/S0038-092X(02)00121-4]"
                )

        # if it is night
        if not self.solar_angles["day"][time]:
            # non-zero radiation during night check
            if self.meteo_vars["glob_hor_rad"][time] != 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation during night"
                )

            # non-zero radiation during night check
            if self.meteo_vars["dir_norm_rad"][time] != 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation during night"
                )

            # non-zero radiation during night check
            if self.meteo_vars["diff_hor_rad"][time] != 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation during night"
                )

    def db_hour_checks(self, time):

        # --------------------
        # DRY BULB TEMPERATURE
        # --------------------

        # impossible physics check
        if (
            self.meteo_vars["dry_bulb"][time] < -90
            or self.meteo_vars["dry_bulb"][time] > 60
        ):

            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is physically impossible or exceedingly rare, breaking world records "
                "[Guide on the Global Data-Processing System (WMO, 1993)"
            )

        # rare check
        if (
            self.meteo_vars["dry_bulb"][time] < -80
            or self.meteo_vars["dry_bulb"][time] > 50
        ):
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is rare [Guide on the Global Data-Processing System (WMO, 1993)]"
            )

        # 1st and 99th percentile check
        month = int(time.month) - 1

        perc_1st_idx = min(
            range(len(self.historical["dry_bulb_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dry_bulb_pdf"]["cdf"][month][i] - 0.01),
        )
        perc_1st = self.historical["dry_bulb_pdf"]["bin"][month][perc_1st_idx]

        perc_99th_idx = min(
            range(len(self.historical["dry_bulb_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dry_bulb_pdf"]["cdf"][month][i] - 0.99),
        )
        perc_99th = self.historical["dry_bulb_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["dry_bulb"][time] < perc_1st:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the lower 1st percentile"
            )

        if self.meteo_vars["dry_bulb"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the upper 99th percentile"
            )

    def dp_hour_check(self, time):

        # --------------------
        # DEWPOINT TEMPERATURE
        # --------------------

        # impossible physics check
        # todo: check that it is inferior to dry bulb

        # todo: check that it is superior to the dew_point as calculated by the dry bulb at 0% relative humidity

        month = int(time.month) - 1

        # 1st and 99th percentile check
        perc_1st_idx = min(
            range(len(self.historical["dew_point_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dew_point_pdf"]["cdf"][month][i] - 0.01),
        )
        perc_1st = self.historical["dew_point_pdf"]["bin"][month][perc_1st_idx]

        perc_99th_idx = min(
            range(len(self.historical["dew_point_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dew_point_pdf"]["cdf"][month][i] - 0.99),
        )
        perc_99th = self.historical["dew_point_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["dew_point"][time] < perc_1st:
            self.meteo_vars["hourly_flags"][time].append(
                "dew point is in the lower 1st percentile"
            )

        if self.meteo_vars["dew_point"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "dew point is in the upper 99th percentile"
            )

    def rh_hour_check(self, time):

        # -----------------
        # RELATIVE HUMIDITY
        # -----------------

        # impossible physics check
        if (
            self.meteo_vars["rel_hum"][time] < 0
            or self.meteo_vars["rel_hum"][time] > 100
        ):

            self.meteo_vars["hourly_flags"][time].append(
                "relative humidity is physically impossible [Guidelines on quality control procedures for data from "
                "automatic weather stations (Zahumenskỳ, 2004)]"
            )

        month = int(time.month) - 1

        # 1st and 99th percentile check
        perc_1st_idx = min(
            range(len(self.historical["rel_hum_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["rel_hum_pdf"]["cdf"][month][i] - 0.01),
        )
        perc_1st = self.historical["rel_hum_pdf"]["bin"][month][perc_1st_idx]

        perc_99th_idx = min(
            range(len(self.historical["rel_hum_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["rel_hum_pdf"]["cdf"][month][i] - 0.99),
        )
        perc_99th = self.historical["rel_hum_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["rel_hum"][time] < perc_1st:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the lower 1st percentile"
            )

        if self.meteo_vars["rel_hum"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the upper 99th percentile"
            )

    def ws_hour_check(self, time):
        # 99th percentile check
        month = int(time.month) - 1

        perc_99th_idx = min(
            range(len(self.historical["wind_speed_pdf"]["cdf"][month])),
            key=lambda i: abs(
                self.historical["wind_speed_pdf"]["cdf"][month][i] - 0.99
            ),
        )
        perc_99th = self.historical["wind_speed_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["wind_speed"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "wind speed is in the upper 99th percentile"
            )

    def step_qc(self):

        for time, _ in list(self.meteo_vars.iterrows())[:-1]:
            self.db_step_check(time)
            self.rh_step_check(time)
            self.ws_step_check(time)

        print("done hourly step checks")

    def db_step_check(self, time):

        # --------------------
        # DRY BULB TEMPERATURE
        # --------------------

        idx = self.meteo_vars.index.get_loc(time)
        db_step = (
            self.meteo_vars["dry_bulb"][idx + 1] - self.meteo_vars["dry_bulb"][idx]
        )

        # Exceedingly rare check
        if db_step > 8:
            self.meteo_vars["step_flags"][idx].append(
                "dry bulb temperature step is rare [Handbook of Automated Data Quality Control Checks and Procedures "
                "(NOAA, 2009)]"
            )
            self.meteo_vars["step_flags"][idx + 1].append(
                "dry bulb temperature step is rare [Handbook of Automated Data Quality Control Checks and Procedures "
                "(NOAA, 2009)]"
            )

    def rh_step_check(self, time):

        # -----------------
        # RELATIVE HUMIDITY
        # -----------------

        idx = self.meteo_vars.index.get_loc(time)
        rh_step = self.meteo_vars["rel_hum"][idx + 1] - self.meteo_vars["rel_hum"][idx]

        # Exceedingly rare check
        if rh_step > 30:
            self.meteo_vars["step_flags"][idx].append(
                "relative humidity step is rare [Quality Control of Meteorological Observations-Automatic "
                "Methods Used in the Nordic Countries (Vejen, 2002)]"
            )
            self.meteo_vars["step_flags"][idx + 1].append(
                "relative humidity step is rare [Quality Control of Meteorological Observations-Automatic "
                "Methods Used in the Nordic Countries (Vejen, 2002)]"
            )

    def ws_step_check(self, time):

        # -----------------
        # RELATIVE HUMIDITY
        # -----------------

        idx = self.meteo_vars.index.get_loc(time)
        ws_step = (
            self.meteo_vars["wind_speed"][idx + 1] - self.meteo_vars["wind_speed"][idx]
        )

        # Exceedingly rare check
        if ws_step > 15:
            self.meteo_vars["step_flags"][idx].append("wind speed step is rare [8 9]")
            self.meteo_vars["step_flags"][idx + 1].append(
                "wind speed step is rare [Handbook of Automated Data Quality Control Checks and Procedures "
                "(NOAA, 2009)]"
            )

    def day_qc(self):
        for time in pd.date_range(
            start=self.meteo_vars.index[0], end=self.meteo_vars.index[-1], freq="D"
        ):
            self.db_day_check(time)

        print("done dry bulb daily profile checks")

    def db_day_check(self, time):
        # get all indices
        idx = self.meteo_vars.index.get_loc(time)
        day_slice = slice(idx, idx + 24, 1)
        month = int(time.month) - 1

        # start a counter with the number of hours that are in the 1st and 99th percentile
        suspect_flag_count = []

        for i, (hour, _) in enumerate(self.meteo_vars[day_slice].iterrows()):
            # 1st and 99th percentile check
            perc_1st_idx = min(
                range(len(self.historical["dry_bulb_time_pdf"]["cdf"][month][i])),
                key=lambda j: abs(
                    self.historical["dry_bulb_time_pdf"]["cdf"][month][i][j] - 0.01
                ),
            )
            perc_1st = self.historical["dry_bulb_time_pdf"]["bin"][month][i][
                perc_1st_idx
            ]

            perc_99th_idx = min(
                range(len(self.historical["dry_bulb_time_pdf"]["cdf"][month][i])),
                key=lambda j: abs(
                    self.historical["dry_bulb_time_pdf"]["cdf"][month][i][j] - 0.99
                ),
            )
            perc_99th = self.historical["dry_bulb_time_pdf"]["bin"][month][i][
                perc_99th_idx
            ]

            if self.meteo_vars["dry_bulb"][hour] < perc_1st:
                suspect_flag_count.append(hour)

            if self.meteo_vars["dry_bulb"][hour] > perc_99th:
                suspect_flag_count.append(hour)

        # if six or more hours during the day are suspect because their corresponding temperatures fall in either the
        # 1st and 99th percentile, flag as a suspect daily profile
        if len(suspect_flag_count) >= 6:
            for t in suspect_flag_count:
                self.meteo_vars["daily_flags"][t].append(
                    "dry bulb temperature daily profile is suspect"
                )

    def month_qc(self):
        for time in pd.date_range(
            start=self.meteo_vars.index[0], periods=12, freq="MS"
        ):
            # self.prec_month_check(time)
            self.rad_month_check(time)

        print("done monthly value checks")

    def prec_month_check(self, time):
        # -------------------
        # WATER PRECIPITATION
        # -------------------

        # get all indices
        idx = self.meteo_vars.index.get_loc(time)
        days_in_month = calendar.monthrange(time.year, time.month)[1]
        month_slice = slice(idx, idx + days_in_month * 24, 1)
        month = int(time.month) - 1

        # max and min check
        if (
            sum(self.meteo_vars["liq_prec_dep"][month_slice])
            < self.historical["prec_dep_min"][month]
        ):
            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "liquid precipitation depth is less than the minimum observed"
                )

        if (
            sum(self.meteo_vars["liq_prec_dep"][month_slice])
            > self.historical["prec_dep_max"][month]
        ):

            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "liquid precipitation depth is more than the maximum observed"
                )

        # 1st and 99th percentile check
        perc_1st = (
            self.historical["prec_dep_avg"][month]
            - 2.32635 * self.historical["prec_dep_std"][month]
        )
        perc_99th = (
            self.historical["prec_dep_avg"][month]
            + 2.32635 * self.historical["prec_dep_std"][month]
        )

        if sum(self.meteo_vars["liq_prec_dep"][month_slice]) < perc_1st:
            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "liquid precipitation depth is in the lower 1st percentile"
                )

        if sum(self.meteo_vars["liq_prec_dep"][month_slice]) > perc_99th:
            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "liquid precipitation depth is in the upper 99th percentile"
                )

    def rad_month_check(self, time):
        # -----------------------
        # RADIATION / IRRADIANCE
        # -----------------------
        # (since radiation is in units of Wh/m^2 in a time series with hour time steps, both
        # are equivalent)

        # todo: add smaller than 0.03 * GHItoa (G_o)

        # get all indices
        idx = self.meteo_vars.index.get_loc(time)
        days_in_month = calendar.monthrange(time.year, time.month)[1]
        days_sum = []
        for day in pd.date_range(start=time, periods=days_in_month, freq="D"):
            day_idx = self.meteo_vars.index.get_loc(day)
            day_slice = slice(day_idx, day_idx + 24, 1)
            days_sum.append(sum(self.meteo_vars["glob_hor_rad"][day_slice]))
        month = int(time.month) - 1

        # # calculate extraterrestrial solar radiation
        # G_sc = 1361
        # n = time.dayofyear + 20 + (11.5) / 24
        # G_on = G_sc * (1 + 0.0334 * math.cos(math.radians((360 * (n - 3)) / 365)))
        # G_o = G_on * math.cos(self.solar_angles["zenith"][time])

        # ---------------------------------------
        # global horizontal radiation/illuminance
        # ---------------------------------------

        # 1st and 99th percentile check
        perc_all_1st = (
            self.historical["all_sky_glob_avg"][month]
            - 2.32635 * self.historical["all_sky_glob_std"][month]
        )
        perc_all_99th = (
            self.historical["all_sky_glob_avg"][month]
            + 2.32635 * self.historical["all_sky_glob_std"][month]
        )

        if np.mean(days_sum) / 1000 < perc_all_1st:
            for hour in range(idx, idx + days_in_month * 24):
                self.meteo_vars["monthly_flags"][hour].append(
                    "global horizontal radiation daily mean in the month is in the lower 1st percentile"
                )

        if np.mean(days_sum) / 1000 > perc_all_99th:

            for hour in range(idx, idx + days_in_month * 24):
                self.meteo_vars["monthly_flags"][hour].append(
                    "global horizontal radiation daily mean in the month is in the upper 99th percentile"
                )

        # # -----------------------------------
        # # direct normal radiation/illuminance
        # # -----------------------------------
        #
        # daylight_hours = sum(self.solar_angles["day"][month_slice])
        # day_adj_perc_1st = perc_all_1st * (days_in_month * 24) / daylight_hours
        # day_adj_perc_99th = perc_all_99th * (days_in_month * 24) / daylight_hours
        #
        # # Use Erbs model to split global horizontal average in diffuse and direct
        # k_t_1st = day_adj_perc_1st / G_o
        # Id_over_I_1st = self.get_Id_over_I(k_t_1st)
        #
        # k_t_99th = day_adj_perc_1st / G_o
        # Id_over_I_99th = self.get_Id_over_I(k_t_99th)
        #
        # perc_dir_1st =
        #
        #
        # ab = (
        #     1.454
        #     - 0.406 * self.historical["dir_nor_taub"][month_idx]
        #     - 0.268 * self.historical["diff_hor_taud"][month_idx]
        #     + 0.021
        #     * self.historical["dir_nor_taub"][month_idx]
        #     * self.historical["diff_hor_taud"][month_idx]
        # )
        #
        # ad = (
        #     0.507
        #     + 0.205 * self.historical["dir_nor_taub"][month_idx]
        #     - 0.08 * self.historical["diff_hor_taud"][month_idx]
        #     - 0.19
        #     * self.historical["dir_nor_taub"][month_idx]
        #     * self.historical["diff_hor_taud"][month_idx]
        # )
        #
        # # 1st and 99th percentile check
        # perc_dir_1st = perc_all_1st * x
        # perc_dir_99th = perc_all_99th * x
        #
        # if sum(self.meteo_vars["glob_hor_rad"][month_slice]) < perc_dir_1st:
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )
        #
        # if sum(self.meteo_vars["glob_hor_rad"][month_slice]) > perc_dir_99th:
        #
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )
        #
        # # ----------------------------------------
        # # diffuse horizontal radiation/illuminance
        # # ----------------------------------------
        #
        # # 1st and 99th percentile check
        # # perc_diff_1st = perc_all_1st * x
        # # perc_diff_99th = perc_all_99th * x
        #
        # if sum(self.meteo_vars["glob_hor_rad"][month_slice]) < perc_diff_1st:
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )
        #
        # elif sum(self.meteo_vars["glob_hor_rad"][month_slice]) > perc_diff_99th:
        #
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )

        # def get_Id_over_I(self, k_t):
        #     if k_t <= 0.22:
        #         Id_over_I = 1.0 - 0.09 * k_t
        #     elif k_t > 0.22 and k_t <= 0.80:
        #         Id_over_I = 0.9511 - 0.1604 * k_t + 4.388 * k_t ** 2 - 16.638 * k_t ** 3 + 12.336 * k_t ** 4
        #     elif k_t > 0.8:
        #         Id_over_I = 0.165

    def get_graph(self, check):

        def get_xlabel(meteo_var):
            """
            Get string with full meteorological variable name and units for xlabel of graphs

            @type: basestring
            @param meteo_var: the meteorological variable for the xlabel
            @rtype: basestring
            @return: xlabel for graph
            """
            # get full xlabel names
            if meteo_var == "dry_bulb":
                xlabel = "Dry bulb temperature [°C]"

            if meteo_var == "dew_point":
                xlabel = "Dew point temperature [°C]"

            if meteo_var == "rel_hum":
                xlabel = "Relative humidity [%]"

            if meteo_var == "glob_hor_rad":
                xlabel = r"Global horizontal irradiation [$\frac{kWh}{m^2}$]"

            if meteo_var == "dir_norm_rad":
                xlabel = r"Direct normal irradiation [$\frac{kWh}{m^2}$]"

            if meteo_var == "diff_hor_rad":
                xlabel = r"Diffuse horizontal irradiation [$\frac{kWh}{m^2}$]"

            if meteo_var == "wind_speed":
                xlabel = r"Wind speed [$\frac{m}{s}$]"

            if meteo_var == "liq_prec_dep":
                xlabel = "Liquid precipitation depth [mm]"

            return xlabel

        def get_low_high_div(min_v, max_v, num_div, fixed_low=False, fixed_high=False):

            # cast integer into floats
            if isinstance(min_v, int):
                min_v = float(min_v)

            if isinstance(max_v, int):
                max_v = float(max_v)

            div = (max_v - min_v) / (num_div - 2)
            digits = "%e" % div

            lead_digit = float(digits.split("e")[0])
            exp = int(digits.split("e")[1])

            if lead_digit >= 6.5:
                lead_digit = 8
            elif 6.5 > lead_digit >= 4.5:
                lead_digit = 5
            elif 4.5 > lead_digit >= 3:
                lead_digit = 4
            elif 3 > lead_digit >= 1.5:
                lead_digit = 2
            else:
                lead_digit = 1

            # try moving up and then down the list of possible leading digits
            dig_list = [1, 2, 4, 5, 8]

            i_d = dig_list.index(lead_digit)

            # changes = [1, -1, 2, -2]
            changes = [-3, -2, -1, 0, 1, 2, 3]

            low_list = []
            high_list = []
            div_list = []

            for chg_1 in changes:
                for chg_2 in changes:
                    i_1 = (i_d + chg_1) % len(dig_list)
                    if i_d + chg_1 > len(dig_list) - 1:
                        e_1 = exp + 1
                    elif i_d + chg_1 < 0:
                        e_1 = exp - 1
                    else:
                        e_1 = exp

                    d_1 = float(f"{dig_list[i_1]}e{e_1}")

                    i_2 = (i_d + chg_2) % len(dig_list)
                    if i_d + chg_2 > len(dig_list) - 1:
                        e_2 = exp + 1
                    elif i_d + chg_2 < 0:
                        e_2 = exp - 1
                    else:
                        e_2 = exp

                    d_2 = float(f"{dig_list[i_2]}e{e_2}")

                    if d_1 == d_2:
                        low_try = [d_1]
                        high_try = [d_1]
                        ticks_try = [d_1]

                    else:
                        low_try = [d_2, d_2, d_2, d_2, d_1, d_1, d_1, d_1]
                        high_try = [d_2, d_2, d_1, d_1, d_2, d_2, d_1, d_1]
                        ticks_try = [d_2, d_1, d_2, d_1, d_2, d_1, d_2, d_1]

                    for k in range(len(low_try)):
                        digs_low = "%e" % low_try[k]
                        e_low = int(digs_low.split("e")[1])

                        if fixed_low:
                            low = min_v
                        else:
                            low = low_try[k] * math.floor(min_v / low_try[k])
                            if low == min_v:
                                low = min_v - low_try[k]
                            low = round(low, -min(0, e_low))

                        digs_high = "%e" % high_try[k]
                        e_high = int(digs_high.split("e")[1])

                        if fixed_high:
                            high = max_v
                        else:
                            high = high_try[k] * math.ceil(max_v / high_try[k])
                            if high == max_v:
                                high = max_v + high_try[k]
                            high = round(high, -min(0, e_high))

                        ticks = round((high - low) / ticks_try[k] + 1, -min(0, e_low, e_high))
                        if ticks == num_div:
                            low_list.append(low)
                            high_list.append(high)
                            div_list.append(ticks_try[k])

            # now select the triplet (low, high, div) where the division has the least number of decimal places
            # if tied, choose the first one in the list
            dcm_places = np.array(
                [
                    sum(
                        [
                            0 if f.is_integer() else str(f)[::-1].find(".")
                            for f in [low_list[i], high_list[i], div_list[i]]
                        ]
                    )
                    for i in range(len(low_list))
                ]
            )
            min_idx = np.where(dcm_places == dcm_places.min())
            low_list_min = np.array(low_list)[min_idx]
            high_list_min = np.array(high_list)[min_idx]
            div_list_min = np.array(div_list)[min_idx]

            # finally select the smallest range
            spread = high_list_min - low_list_min
            final_idx = np.argmin(spread)
            final_low = low_list_min[final_idx]
            final_high = high_list_min[final_idx]
            final_div = div_list_min[final_idx]

            return final_low, final_high, final_div

        # create graphs
        if check == "hourly":

            # get flagged values for the meteorological variable
            for meteo_var in ["dry_bulb", "dew_point", "rel_hum", "wind_speed"]:
                string_check = meteo_var.split("_")[0]
                value_list = []
                time_list = []

                for time, value in self.meteo_vars[meteo_var].iteritems():
                    for flag in self.meteo_vars["hourly_flags"][time]:
                        if re.match(rf".*{string_check}", flag):
                            if "missing" not in flag:
                                value_list.append(value)
                                time_list.append(time)

                # get all days that have hourly flags; days could repeat
                all_days = [time_idx.day for time_idx in time_list]

                # get indices for unique days
                unique_days_idx = np.unique(all_days, return_index=True)[1]

                # get unique datetimeindices based on unique days
                unique_times = [time_list[i] for i in unique_days_idx]

                for unq_time in unique_times:
                    # create a day slice of the indices to get the daily profile of the meteorological variable
                    day = unq_time.day
                    month = unq_time.month - 1
                    year = unq_time.year

                    midnight = datetime.datetime(
                        year,
                        month + 1,
                        day,
                    )
                    midnight_idx = self.meteo_vars.index.get_loc(midnight)
                    day_slice = slice(midnight_idx, midnight_idx + 25, 1)
                    daily_profile = self.meteo_vars[meteo_var][day_slice]

                    pdf_name = meteo_var + "_pdf"

                    # get 1st percentile and 99th percentile values
                    perc_1st_idx = min(
                        range(len(self.historical[pdf_name]["cdf"][month])),
                        key=lambda i: abs(self.historical[pdf_name]["cdf"][month][i] - 0.01),
                    )
                    perc_1st = self.historical[pdf_name]["bin"][month][perc_1st_idx]

                    perc_99th_idx = min(
                        range(len(self.historical[pdf_name]["cdf"][month])),
                        key=lambda i: abs(self.historical[pdf_name]["cdf"][month][i] - 0.99),
                    )
                    perc_99th = self.historical[pdf_name]["bin"][month][perc_99th_idx]

                    bins = self.historical[pdf_name]["bin"][month]
                    pdf = self.historical[pdf_name]["pdf"][month]
                    fig = plt.figure(figsize=(3.3, 2.1))
                    ax1 = fig.add_subplot(111)
                    ax2 = ax1.twinx()
                    ax1.grid(alpha=0.2)
                    plt.rcParams["font.family"] = "Arial"
                    plt.rcParams['mathtext.fontset'] = "custom"
                    plt.rcParams['mathtext.rm'] = "Arial"
                    csfont = {"fontname": "Arial"}
                    fontsz = {"fontsize": 10}

                    ax1.hist(bins, bins=len(bins), weights=pdf, color="gray")
                    ax2.plot(daily_profile, range(25), color="lightseagreen")
                    ax1.set_ylabel("Probability", **csfont, **fontsz)
                    ax2.set_ylabel("Hour of the day", color="lightseagreen", **csfont, **fontsz)

                    # set xlabel based on meteo_var
                    xlabel = get_xlabel(meteo_var)
                    ax1.set_xlabel(f"{xlabel}", **csfont, **fontsz)

                    x_min = min(min(daily_profile), perc_1st)
                    x_max = max(max(daily_profile), perc_99th)
                    # determine low, high, and divisions for x axis based on range of meteorological variable
                    x_low, x_high, x_div = get_low_high_div(x_min, x_max, 10)
                    if x_div >= 1:
                        ax1.set_xticks(np.arange(x_low, x_high + x_div, x_div).astype(int))
                        ax1.set_xlim([int(x_low), int(x_high)])
                    else:
                        ax1.set_xticks(np.arange(x_low, x_high + x_div, x_div))
                        ax1.set_xlim([x_low, x_high])

                    y_min = 0
                    y_max = max(pdf)
                    # determine low, high ,and divisions for y axis based on range of historical bins of meteorological
                    # variable
                    y_low, y_high, y_div = get_low_high_div(y_min, y_max, 8, fixed_low=True)
                    ax1.set_yticks(np.arange(y_low, y_high + y_div, y_div))
                    ax1.set_ylim([y_low, y_high])
                    ax2.set_yticks(np.arange(0, 27, 3).astype(int))
                    ax2.set_ylim([0, 24])

                    ax1.set_xticklabels(ax1.get_xticks(), **csfont, **fontsz)
                    ax1.set_yticklabels(ax1.get_yticks(), **csfont, **fontsz)
                    ax2.set_yticklabels(ax2.get_yticks(), color="lightseagreen", **csfont, **fontsz)
                    plt.vlines(perc_1st, 0, 24, colors="coral")
                    plt.vlines(perc_99th, 0, 24, colors="coral")

                    # if self.weather_path.stem == "CAN_QC_Montreal-McTavish.716120_CWEC2016":
                    #     plt.title(
                    #         "Typical Meteorological Year", fontweight="bold", **csfont, **fontsz
                    #     )
                    #
                    # if self.weather_path.stem == "tampered":
                    #     plt.title("Tampered Year", fontweight="bold", **csfont, **fontsz)

                    plt.title(f"{self.weather_path.stem}\nday: {day}, month: {month + 1}", fontweight="bold", **csfont, **fontsz)

                    plt.tight_layout(pad=0.05)
                    plt.savefig(
                        Path(f"datafiles/figures/hourly_flags_{meteo_var}_day{day}_month{month + 1}_{self.weather_path.stem}.svg")
                    )

            print("done making hourly value and hourly step graphs")

        elif check == "daily":

            # get flagged daily values for the meteorological variable
            for meteo_var in ["dry_bulb"]:
                string_check = meteo_var.split("_")[0]
                value_list = []
                time_list = []

                for time, value in self.meteo_vars[meteo_var].iteritems():
                    for flag in self.meteo_vars["daily_flags"][time]:
                        if re.match(rf".*{string_check}", flag):
                            if "missing" not in flag:
                                value_list.append(value)
                                time_list.append(time)

                # get all days that have hourly flags; days could repeat
                all_days = [time_idx.day for time_idx in time_list]

                # get indices for unique days
                unique_days_idx = np.unique(all_days, return_index=True)[1]

                # get unique datetimeindices based on unique days
                unique_times = [time_list[i] for i in unique_days_idx]

                for unq_time in unique_times:
                    # create a day slice of the indices to get the daily profile of the meteorological variable
                    day = unq_time.day
                    month = unq_time.month - 1
                    year = unq_time.year

                    midnight = datetime.datetime(
                        year,
                        month + 1,
                        day,
                    )
                    midnight_idx = self.meteo_vars.index.get_loc(midnight)
                    day_slice = slice(midnight_idx, midnight_idx + 25, 1)
                    daily_profile = self.meteo_vars[meteo_var][day_slice]

                    pdf_name = meteo_var + "_time_pdf"

                    plt.figure(figsize=(3.3, 5))
                    gs1 = gridspec.GridSpec(24, 1)
                    gs1.update(wspace=0, hspace=0)  # set the spacing between axes.
                    plt.rcParams["font.family"] = "Arial"
                    plt.rcParams['mathtext.fontset'] = 'custom'
                    plt.rcParams['mathtext.rm'] = "Arial"
                    csfont = {"fontname": "Arial"}
                    fontsz = {"fontsize": 10}

                    # initialize lists to determine x and y limits of plots
                    perc_99th_list = []
                    perc_1st_list = []
                    bins_list = []
                    pdf_list = []

                    for hour in range(24):
                        # get 1st percentile and 99th percentile values
                        perc_1st_idx = min(
                            range(len(self.historical[pdf_name]["cdf"][month][hour])),
                            key=lambda i: abs(
                                self.historical[pdf_name]["cdf"][month][hour][i] - 0.01
                            ),
                        )
                        perc_1st = self.historical[pdf_name]["bin"][month][hour][perc_1st_idx]

                        perc_99th_idx = min(
                            range(len(self.historical[pdf_name]["cdf"][month][hour])),
                            key=lambda i: abs(
                                self.historical[pdf_name]["cdf"][month][hour][i] - 0.99
                            ),
                        )
                        perc_99th = self.historical[pdf_name]["bin"][month][hour][perc_99th_idx]

                        bins = self.historical[pdf_name]["bin"][month][hour]
                        pdf = self.historical[pdf_name]["pdf"][month][hour]

                        # record value of 99th and 1st percentile to calculate limits of x axis later
                        perc_99th_list.append(perc_99th)
                        perc_1st_list.append(perc_1st)
                        bins_list.append(bins)
                        pdf_list.append(pdf)

                    # determine maximum and minimum limits along with divisions for y and x axes
                    x_min = min(min(perc_1st_list), min(daily_profile))
                    x_max = max(max(perc_99th_list), max(daily_profile))
                    x_low, x_high, x_div = get_low_high_div(x_min, x_max, 10)

                    y_min = 0
                    y_max = max((max(pdf) for pdf in pdf_list))
                    y_low, y_high, y_div = get_low_high_div(y_min, y_max, 8, fixed_low=True)

                    # make a subplot of the histogram for each hour of the month
                    for hour in range(24):

                        plt.axis("on")
                        ax1 = plt.subplot(gs1[23 - hour])
                        ax2 = ax1.twinx()
                        ax1.xaxis.grid(True, alpha=0.2)
                        ax1.hist(bins_list[hour], bins=len(bins_list[hour]), weights=pdf_list[hour], color="gray")
                        ax2.plot(
                            daily_profile[hour : hour + 2],
                            np.arange(hour, hour + 2, 1),
                            color="lightseagreen",
                        )
                        if hour == 12:
                            ax1.set_ylabel("Probability", **csfont, **fontsz)
                            ax1.set_yticks([y_low, y_high])
                            ml = MultipleLocator(y_high / 2)
                            ax1.yaxis.set_minor_locator(ml)
                            ax1.set_yticklabels([f"{y}" if y != 0 else "0" for y in ax1.get_yticks()], **csfont,
                                                **fontsz)
                            ax2.set_ylabel("Hour of the day", **csfont, **fontsz, color="lightseagreen")
                        else:
                            for tic in ax1.yaxis.get_major_ticks():
                                tic.tick1line.set_visible(False)
                                tic.tick2line.set_visible(False)
                                tic.label1.set_visible(False)
                                tic.label2.set_visible(False)

                        ax1.set_xlim([x_low, x_high])
                        ax1.set_ylim([y_low, y_high])

                        ax2.set_ylim([hour, hour + 1])
                        ax2.set_yticklabels(
                            [f"{int(y)}" for y in ax2.get_yticks()],
                            color="lightseagreen",
                            **csfont,
                            **fontsz,
                        )

                        if hour == 0:
                            ax2.set_yticks([int(hour)])
                            xlabel = get_xlabel(meteo_var)
                            if x_div >= 1:
                                ax1.set_xticks(np.arange(x_low, x_high + x_div, x_div).astype(int))
                            else:
                                ax1.set_xticks(np.arange(x_low, x_high + x_div, x_div))
                            ax1.set_xlabel(xlabel, **csfont, **fontsz)
                            ax1.set_xticklabels(ax1.get_xticks(), **csfont, **fontsz)
                            ml = MultipleLocator(x_div / 2)
                            ax1.xaxis.set_minor_locator(ml)
                        else:
                            for tic in ax1.xaxis.get_major_ticks():
                                tic.tick1line.set_visible(False)
                                tic.tick2line.set_visible(False)
                                tic.label1.set_visible(False)
                                tic.label2.set_visible(False)
                            for tic in ax1.xaxis.get_minor_ticks():
                                tic.tick1line.set_visible(False)
                                tic.tick2line.set_visible(False)
                                # tic.label1.set_visible(False)
                                # tic.label2.set_visible(False)

                            ax2.set_yticks([int(hour), int(hour + 1)])

                        ax2.vlines(perc_1st_list[hour], hour, hour + 1, colors="coral")
                        ax2.vlines(perc_99th_list[hour], hour, hour + 1, colors="coral")

                    # if self.weather_path.stem == "CAN_QC_Montreal-McTavish.716120_CWEC2016":
                    #     plt.title(
                    #         "Typical Meteorological Year", fontweight="bold", **csfont, **fontsz
                    #     )
                    #
                    # if self.weather_path.stem == "tampered":
                    #     plt.title("Tampered Year", fontweight="bold", **csfont, **fontsz)

                    plt.title(f"{self.weather_path.stem}\nday: {day}, month: {month + 1}", fontweight="bold", **csfont,
                              **fontsz)

                    plt.tight_layout(pad=0.05)
                    plt.savefig(
                        Path(f"datafiles/figures/daily_flags_{meteo_var}_day{day}_month{month + 1}_{self.weather_path.stem}.svg")
                    )

            print("done making daily profile graphs")

        elif check == "monthly":

            for meteo_var in ["glob_hor_rad", "liq_prec_dep"]:
                string_check = meteo_var.split("_")[0]
                make_graph = False

                # determine if there is a month that is extreme, rare, or suspect
                for time, value in self.meteo_vars[meteo_var].iteritems():
                    for flag in self.meteo_vars["monthly_flags"][time]:
                        if re.match(rf".*{string_check}", flag):
                            if "missing" not in flag:
                                make_graph = True

                if make_graph:
                    plt.figure(figsize=(3.3, 4))
                    gs1 = gridspec.GridSpec(12, 1)
                    gs1.update(wspace=0, hspace=0)  # set the spacing between axes.
                    plt.rcParams["font.family"] = "Arial"
                    plt.rcParams['mathtext.fontset'] = 'custom'
                    plt.rcParams['mathtext.rm'] = "Arial"

                    csfont = {"fontname": "Arial"}
                    fontsz = {"fontsize": 10}

                    past_days = 364

                    for month in range(11, -1, -1):
                        mu = self.historical["all_sky_glob_avg"][month]
                        sigma = self.historical["all_sky_glob_std"][month]

                        perc_1st = mu - 2.32635 * sigma
                        perc_99th = mu + 2.32635 * sigma

                        z1 = mu - 3 * sigma
                        z2 = mu + 3 * sigma

                        x = np.arange(z1, z2, 0.001)
                        y = norm.pdf(x, mu, sigma)

                        # get all indices
                        days_in_month = calendar.monthrange(
                            self.meteo_vars.index[0].year, month + 1
                        )[1]
                        days_sum = []
                        for day in range(past_days, past_days - days_in_month, -1):
                            day_slice = slice(day * 24, (day + 1) * 24, 1)
                            days_sum.append(sum(self.meteo_vars[meteo_var][day_slice]))
                        past_days -= days_in_month

                        month_average = np.mean(days_sum) / 1000

                        plt.axis("on")
                        ax1 = plt.subplot(gs1[month])
                        ax2 = ax1.twinx()
                        ax1.xaxis.grid(True, alpha=0.2, which="both")
                        ax1.barh(0, month_average, align="center", color="lightseagreen")
                        ax2.plot(x, y, color="gray")

                        x_min = 0
                        x_max = max(x)
                        x_low, x_high, x_div = get_low_high_div(x_min, x_max, 10, fixed_low=True)

                        if x_div >= 1:
                            ax1.set_xticks(np.arange(x_low, x_high + x_div, x_div).astype(int))
                            ax1.set_xlim([int(x_low), int(x_high)])
                        else:
                            ax1.set_xticks(np.arange(x_low, x_high + x_div, x_div))
                            ax1.set_xlim([x_low, x_high])

                        ml = MultipleLocator(0.5)
                        ax1.xaxis.set_minor_locator(ml)

                        ax1.set_yticks([-1, 0, 1])
                        ax1.set_yticklabels(
                            ["", calendar.month_name[month + 1], ""],
                            color="lightseagreen",
                            **csfont,
                            **fontsz,
                        )
                        for i, tic in enumerate(ax1.yaxis.get_major_ticks()):
                            if i != 1:
                                tic.label1.set_visible(False)
                                tic.label2.set_visible(False)
                            tic.tick1line.set_visible(False)
                            tic.tick2line.set_visible(False)

                        if month == 11:
                            xlabel = get_xlabel(meteo_var)
                            ax1.set_xlabel(
                                "Monthly average"
                                "\n"
                                f"{xlabel}",
                                **csfont,
                                **fontsz,
                            )
                            ax1.set_xticklabels(ax1.get_xticks(), **csfont, **fontsz)

                        else:
                            for tic in ax1.xaxis.get_major_ticks():
                                tic.tick1line.set_visible(False)
                                tic.tick2line.set_visible(False)
                                tic.label1.set_visible(False)
                                tic.label2.set_visible(False)
                            for tic in ax1.xaxis.get_minor_ticks():
                                tic.tick1line.set_visible(False)
                                tic.tick2line.set_visible(False)
                                tic.label1.set_visible(False)
                                tic.label2.set_visible(False)

                        ax2.get_yaxis().set_visible(False)

                        ax2.vlines(perc_1st, 0, max(y), colors="coral")
                        ax2.vlines(perc_99th, 0, max(y), colors="coral")

                    # if self.weather_path.stem == "CAN_QC_Montreal-McTavish.716120_CWEC2016":
                    #     plt.title(
                    #         "Typical Meteorological Year", fontweight="bold", **csfont, **fontsz
                    #     )
                    #
                    # if self.weather_path.stem == "CAN-QC - Montreal YUL 716270 - ISD 2015":
                    #     plt.title("2015", fontweight="bold", **csfont, **fontsz)

                    plt.title(self.weather_path.stem, fontweight="bold", **csfont,
                              **fontsz)

                    plt.tight_layout(pad=0.05)
                    plt.savefig(
                        Path(f"datafiles/figures/monthly_flags_{meteo_var}_{self.weather_path.stem}.svg")
                    )

            print("done making monthly value graphs")


def get_weather(weather_path, drive_letter, ashrae_path, output_dir):
    """
    @type weather_path: basestring, pathlib Path
    @param weather_path: the file path (either Windows or Mac/Linux) for the .epw weather file
    @type drive_letter: basestring
    @param drive_letter: the drive letter (without the colon) where the ASHRAE Weather Viewer (2017) iso file is
    mounted
    @type ashrae_path: basestring, pathlib Path
    @param ashrae_path: the file path (either Windows or Mac/Linux) for the ASHRAE Design Conditions (2017)
    Excel file
    @type output_dir: basestring, pathlib Path
    @param output_dir: the output folder/directory where intermediate and final results files will be saved
    @return: Weather() object which has meteorological parameters (longitude, latitude, wmo station, ...), meteorological
    variables (dry bulb temperature, solar irradiation, ...), solar angles (zenith, declination), and historical local
    weather data (probability distribution functions and monthly averages + standard deviations for the
    meteorological variables
    """
    from qcweather import Weather

    # instantiate Weather object
    weather = Weather(Path(weather_path), drive_letter, Path(ashrae_path), Path(output_dir))

    # get meteorological parameters (such as longitude and latitude of weather station, dry bulb temperature,
    # and wind speed)
    weather.get_meteo_params()

    # get solar angles
    weather.get_solar_angles()

    # gather historical local weather data
    weather.get_historical()

    return weather

def get_month_data(excel_sheet, wmo, var_type):

    # find the corresponding row of data for the weather station
    wmo_col = excel_sheet.col_values(4)[4:]
    row = wmo_col.index(str(wmo))

    # read the values from the corresponding columns in the row of the weather station
    if var_type == "clr_sky_dir_norm":
        cells = excel_sheet.row_slice(row, 533, 545)

    elif var_type == "clr_sky_diff_hor":
        cells = excel_sheet.row_slice(row, 545, 557)

    elif var_type == "dir_norm_taub":
        cells = excel_sheet.row_slice(row, 509, 521)

    elif var_type == "diff_hor_taud":
        cells = excel_sheet.row_slice(row, 521, 533)

    elif var_type == "all_sky_glob_avg":
        cells = excel_sheet.row_slice(row, 557, 569)

    elif var_type == "all_sky_glob_std":
        cells = excel_sheet.row_slice(row, 569, 581)

    elif var_type == "prec_dep_avg":
        cells = excel_sheet.row_slice(row, 206, 218)

    elif var_type == "prec_dep_max":
        cells = excel_sheet.row_slice(row, 219, 231)

    elif var_type == "prec_dep_min":
        cells = excel_sheet.row_slice(row, 232, 244)

    elif var_type == "prec_dep_std":
        cells = excel_sheet.row_slice(row, 245, 257)

    else:
        raise ValueError(
            "rad_type must be either 'clr_sky_dir_norm', 'clr_sky_diff_hor', 'dir_norm_taub', 'diff_hor_taud', "
            "'all_sky_glob_avg', 'all_sky_glob_std', 'prec_dep_avg', 'prec_dep_max', 'prec_dep_min', or "
            "'prec_dep_std'"
        )

    # turn N/A values into None and other numerical strings into floats
    month_data = [None if c.value == "N/A" else float(c.value) for c in cells]

    return month_data
