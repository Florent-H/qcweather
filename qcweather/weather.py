import csv
import itertools
import math
import os
import re
import subprocess

import pandas as pd
import numpy as np
import matplotlib as plt
import calendar
import datetime
from pathlib import Path
import epw
from collections import OrderedDict
import xlrd


class Weather(object):
    def __init__(self, file_path):
        self.file_path = file_path
        self.meteo_vars = None
        self.wmo = None
        self.folder = None
        self.location = None
        self.time_zone = None
        self.solar_angles = None
        self.historical = None

    @classmethod
    def get_weather(cls, file_path):

        # instantiate Weather object
        weather = cls(file_path)

        # if the weather file is an energy plus weather file
        if file_path.suffix == ".epw":
            # open energy plus weather file
            with open(file_path, "r") as fp:
                reader = csv.reader(fp)
                # iterate through each line of the energyplus file
                for i, line in enumerate(reader):

                    # for the first line of code, record the World Meteorological Organization station number and create
                    # the folder for the raw historical weather data
                    if re.match(r"(?i)location", line[0]):
                        weather.wmo = int(line[5])
                        weather.folder = Path(f"datafiles/raw_tbl/{weather.wmo}")
                        os.mkdir(weather.folder)
                        # (latitude, longitude) in degrees
                        # using convention that to the west is positive
                        weather.location = (
                            math.radians(float(line[6])),
                            math.radians(-float(line[7])),
                        )
                        weather.time_zone = -int(line[8])
                        weather.solar_angles = get_solar_angles(weather)

                    # if the line starts with 4 digits (the year), then we have passed through the header and can start
                    # recording the values of the meteorological variables to the class's meteo_vars attribute
                    elif re.match(r"\d{4}", line[0]):
                        if not weather.meteo_vars:
                            # get number of hours in the year as it could be a leap year (i.e. 8760 vs 8784 hours)
                            hours = len(list(reader)) + 1
                            first_date = datetime.datetime(
                                int(line[0]), int(line[1]), int(line[2]), int(line[3])
                            )
                            # create dataframe template for meteorological variables for the energyplus weather file
                            weather.meteo_vars = pd.DataFrame(
                                OrderedDict(
                                    (
                                        ("dry_bulb", np.empty(hours)),
                                        ("dew_point", np.empty(hours)),
                                        ("rel_hum", np.empty(hours)),
                                        ("pressure", np.empty(hours)),
                                        ("ext_hor_rad", np.empty(hours)),
                                        ("ext_dir_norm", np.empty(hours)),
                                        ("hor_inf_rad_int", np.empty(hours)),
                                        ("glob_hor_rad", np.empty(hours)),
                                        ("dir_norm_rad", np.empty(hours)),
                                        ("diff_hor_rad", np.empty(hours)),
                                        ("glob_hor_ill", np.empty(hours)),
                                        ("dir_norm_ill", np.empty(hours)),
                                        ("diff_hor_ill", np.empty(hours)),
                                        ("zen_lum", np.empty(hours)),
                                        ("wind_dir", np.empty(hours)),
                                        ("wind_speed", np.empty(hours)),
                                        ("tot_sky_cov", np.empty(hours)),
                                        ("op_sky_cov", np.empty(hours)),
                                        ("visibility", np.empty(hours)),
                                        ("ceil_height", np.empty(hours)),
                                        ("pres_wth_obs", np.empty(hours)),
                                        ("pres_wth_cod", np.empty(hours)),
                                        ("prec_water", np.empty(hours)),
                                        ("aero_opt_dep", np.empty(hours)),
                                        ("snow_depth", np.empty(hours)),
                                        ("days_lst_snow", np.empty(hours)),
                                        ("albedo", np.empty(hours)),
                                        ("liq_prec_dep", np.empty(hours)),
                                        ("liq_prec_qty", np.empty(hours)),
                                        ("hourly_flags", np.empty(hours, dtype=list)),
                                        ("daily_flags", np.empty(hours, dtype=list)),
                                        ("monthly_flags", np.empty(hours, dtype=list)),
                                    )
                                ),
                                index=pd.date_range(
                                    first_date, periods=hours, freq="H"
                                ),
                            )
                            # keep track of the hour number
                            hour = 0

                        # if the meteorological variable dataframe has been created
                        if isinstance(weather.meteo_vars, pd.DataFrame):
                            # write the hour's data to the corresponding columns
                            for j, (col, arr) in enumerate(
                                weather.meteo_vars.iteritems()
                            ):
                                if col not in (
                                    "hourly_flags",
                                    "daily_flags",
                                    "monthly_flags",
                                ):
                                    # missing value
                                    if line[j + 6] == "9999":
                                        weather.meteo_vars["hourly_flags"][hour].append(
                                            "missing value"
                                        )
                                        weather.meteo_vars[col][hour] = None

                                    weather.meteo_vars[col][hour] = float(line[j + 6])
                            hour += 1

        # get solar angles
        weather.solar_angles = get_solar_angles(weather)

        # gather historical data for the meteorological variables

        # run tblexpand to gather raw historical weather data
        run_string = (
            str(Path("E:/tblxpand/tblxpand.exe"))
            + f" {weather.folder} "
            + str(Path(f"E:/data/{weather.wmo}.wdv"))
            + " ALL 99 SI"
        )
        subprocess.run(run_string, cwd=Path.cwd(), capture_output=True)

        # get monthly CDFs
        db_pdf = get_pdf(weather, "dry_bulb")
        wb_pdf = get_pdf(weather, "wet_bulb")
        dp_pdf = get_pdf(weather, "dew_point")
        rel_pdf = get_pdf(weather, "rel_hum")
        wd_pdf = get_pdf(weather, "wind_dir")
        ws_pdf = get_pdf(weather, "wind_speed")

        # get irradiation data
        # todo: generalize workbook location
        w_b = xlrd.open_workbook(Path("datafiles/ashrae/2017DesignConditions_s.xlsx"))

        # read first and only sheet of the excel document
        # todo: generalize selection of sheet index
        sheet_index = 0
        excel_sheet = w_b.sheet_by_index(sheet_index)

        clr_sky_dir_norm = get_month_data(excel_sheet, weather.wmo, "clr_sky_dir_norm")
        clr_sky_diff_hor = get_month_data(excel_sheet, weather.wmo, "clr_sky_diff_hor")
        dir_norm_taub = get_month_data(excel_sheet, weather.wmo, "dir_norm_taub")
        diff_hor_taud = get_month_data(excel_sheet, weather.wmo, "diff_hor_taud")
        all_sky_glob_avg = get_month_data(excel_sheet, weather.wmo, "all_sky_glob_avg")
        all_sky_glob_std = get_month_data(excel_sheet, weather.wmo, "all_sky_glob_std")
        prec_water_avg = get_month_data(excel_sheet, weather.wmo, "prec_water_avg")
        prec_water_max = get_month_data(excel_sheet, weather.wmo, "prec_water_max")
        prec_water_min = get_month_data(excel_sheet, weather.wmo, "prec_water_min")
        prec_water_std = get_month_data(excel_sheet, weather.wmo, "prec_water_std")

        weather.historical = pd.DataFrame(
            OrderedDict(
                (
                    ("dry_bulb_pdf", db_pdf),
                    ("wet_bulb_pdf", wb_pdf),
                    ("dew_point_pdf", dp_pdf),
                    ("rel_hum_pdf", rel_pdf),
                    ("wind_dir_pdf", wd_pdf),
                    ("wind_speed_pdf", ws_pdf),
                    ("clr_sky_dir_norm", clr_sky_dir_norm),
                    ("clr_sky_diff_hor", clr_sky_diff_hor),
                    ("dir_norm_taub", dir_norm_taub),
                    ("diff_hor_taud", diff_hor_taud),
                    ("all_sky_glob_avg", all_sky_glob_avg),
                    ("all_sky_glob_std", all_sky_glob_std),
                    ("prec_water_avg", prec_water_avg),
                    ("prec_water_max", prec_water_max),
                    ("prec_water_min", prec_water_min),
                    ("prec_water_std", prec_water_std),
                )
            ),
        )

        return weather

    def database_qc(self):
        pass

    def magnitude_qc(self):
        # check that meteorological variable values do not exceed thresholds or is negative

        for time, col in self.meteo_vars.iterrows():

            self.rad_mag_checks(time)
            self.db_mag_checks(time)
            self.dp_mag_check(time)
            self.rh_mag_check(time)
            self.ws_mag_check(time)

        for time in pd.date_range(start=self.meteo_vars.index[0], periods=12, freq="M"):
            self.prec_mag_month_check(time)

    def rad_mag_checks(self, time):

        # -----------------------
        # RADIATION / IRRADIANCE
        # -----------------------
        # (since radiation is in units of Wh/m^2 in a time series with hour time steps, both
        # are equivalent)

        # Solar constant
        G_sc = 1361

        # negative radiation checks
        if self.meteo_vars["glob_hor_rad"][time] < 0:
            self.meteo_vars["hourly_flags"][time].append(
                "global horizontal radiation is negative"
            )
        if self.meteo_vars["dir_norm_rad"][time] < 0:
            self.meteo_vars["hourly_flags"][time].append(
                "direct normal radiation is negative"
            )
        if self.meteo_vars["diff_hor_rad"][time] < 0:
            self.meteo_vars["hourly_flags"][time].append(
                "diffuse horizontal radiation is negative"
            )
        # if it is day
        if self.meteo_vars["day"][time]:

            # calculate extraterrestrial solar radiation
            n = time.day - 0.5 + (time.hour - 0.5) / 24
            G_on = G_sc * (1 + 0.0334 * math.cos(math.radians((360 * (n - 3)) / 365)))

            # ---------------------------------------
            # global horizontal radiation/illuminance
            # ---------------------------------------

            # impossible physics check
            if self.meteo_vars["glob_hor_rad"][time] >= min(
                1.2 * G_on,
                1.5 * G_on * math.cos(self.solar_angles["zenith"]) ** 1.2 + 100,
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation is physically impossible [3 4 7]"
                )

            # rare check
            elif (
                self.meteo_vars["glob_hor_rad"][time]
                >= 1.2 * G_on * math.cos(self.solar_angles["zenith"]) ** 1.2 + 50
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation is rare [4 5 7]"
                )

            # 5th and 95th percentile check
            perc_5th = (
                self.historical["all_sky_glob_avg"]
                - 2 * self.historical["all_sky_glob_std"]
            )
            perc_95th = (
                self.historical["all_sky_glob_avg"]
                + 2 * self.historical["all_sky_glob_std"]
            )

            if self.meteo_vars["glob_hor_rad"][time] < perc_5th:
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation is in the lower 5th percentile"
                )

            elif self.meteo_vars["glob_hor_rad"][time] > perc_95th:
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation is in the upper 95th percentile"
                )

            # -----------------------------------
            # direct normal radiation/illuminance
            # -----------------------------------

            # impossible physics check
            if self.meteo_vars["dir_norm_rad"][time] >= G_on:
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is physically impossible [4]"
                )

            # rare check
            elif (
                self.meteo_vars["dir_norm_rad"][time]
                >= 0.95 * G_on * math.cos(self.solar_angles["zenith"]) ** 0.2 + 10
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is rare [4]"
                )

            # -----------------------------------
            # diffuse horizontal radiation/illuminance
            # -----------------------------------

            # impossible physics check
            if self.meteo_vars["diff_hor_rad"][time] >= min(
                0.8 * G_on,
                0.95 * G_on * math.cos(self.solar_angles["zenith"]) ** 1.2 + 50,
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is physically impossible [3 4 5 7]"
                )

            # rare check
            elif (
                self.meteo_vars["diff_hor_rad"][time]
                >= 0.75 * G_on * math.cos(self.solar_angles["zenith"]) ** 1.2 + 30
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is rare [3 4]"
                )

        # if it is night
        else:
            # non-zero radiation during night check
            if self.meteo_vars["glob_hor_rad"][time] > 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation during night"
                )

            # non-zero radiation during night check
            if self.meteo_vars["dir_norm_rad"][time] > 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation during night"
                )

            # non-zero radiation during night check
            if self.meteo_vars["diff_hor_rad"][time] > 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation during night"
                )

    def db_mag_checks(self, time):

        # --------------------
        # DRY BULB TEMPERATURE
        # --------------------

        # impossible (or exceedingly rare) physics check
        if (
            self.meteo_vars["dry_bulb"][time] <= -90
            or self.meteo_vars["dry_bulb"][time] >= 60
        ):

            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is physically impossible or exceedingly rare, breaking world records "
                "[10 14]"
            )

        # rare check
        elif (
            self.meteo_vars["dry_bulb"][time] <= -80
            or self.meteo_vars["dry_bulb"][time] >= 50
        ):
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is rare [10 14]"
            )

        # 5th and 95th percentile check
        month = int(time.month) - 1

        perc_5th_idx = min(
            range(len(self.historical["dry_bulb_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dry_bulb_pdf"]["cdf"][month][i] - 0.05),
        )
        perc_5th = self.historical["dry_bulb_pdf"]["bin"][month][perc_5th_idx]

        perc_95th_idx = min(
            range(len(self.historical["dry_bulb_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dry_bulb_pdf"]["cdf"][month][i] - 0.95),
        )
        perc_95th = self.historical["dry_bulb_pdf"]["bin"][month][perc_95th_idx]

        if self.meteo_vars["dry_bulb"][time] < perc_5th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the lower 5th percentile"
            )

        elif self.meteo_vars["dry_bulb"][time] > perc_95th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the upper 95th percentile"
            )

    def dp_mag_check(self, time):

        # --------------------
        # DEWPOINT TEMPERATURE
        # --------------------
        month = time.month - 1

        # 5th and 95th percentile check
        perc_5th_idx = min(
            range(len(self.historical["dew_point_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dew_point_pdf"]["cdf"][month][i] - 0.05),
        )
        perc_5th = self.historical["dew_point_pdf"]["bin"][month][perc_5th_idx]

        perc_95th_idx = min(
            range(len(self.historical["dew_point_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dew_point_pdf"]["cdf"][month][i] - 0.95),
        )
        perc_95th = self.historical["dew_point_pdf"]["bin"][month][perc_95th_idx]

        if self.meteo_vars["dew_point"][time] < perc_5th:
            self.meteo_vars["hourly_flags"][time].append(
                "dew point is in the lower 5th percentile"
            )

        elif self.meteo_vars["dew_point"][time] > perc_95th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the upper 95th percentile"
            )

    def rh_mag_check(self, time):

        # -----------------
        # RELATIVE HUMIDITY
        # -----------------

        # impossible physics check
        if (
            self.meteo_vars["rel_hum"][time] < 0
            or self.meteo_vars["rel_hum"][time] > 100
        ):

            self.meteo_vars["hourly_flags"][time].append(
                "relative humidity is physically impossible"
            )

        month = time.month - 1

        # 5th and 95th percentile check
        perc_5th_idx = min(
            range(len(self.historical["rel_hum_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["rel_hum_pdf"]["cdf"][month][i] - 0.05),
        )
        perc_5th = self.historical["rel_hum_pdf"]["bin"][month][perc_5th_idx]

        perc_95th_idx = min(
            range(len(self.historical["rel_hum_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["rel_hum_pdf"]["cdf"][month][i] - 0.95),
        )
        perc_95th = self.historical["rel_hum_pdf"]["bin"][month][perc_95th_idx]

        if self.meteo_vars["rel_hum"][time] < perc_5th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the lower 5th percentile"
            )

        elif self.meteo_vars["rel_hum"][time] > perc_95th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the upper 95th percentile"
            )

    def ws_mag_check(self, time):
        # 95th percentile check
        month = time.month - 1

        perc_95th_idx = min(
            range(len(self.historical["wind_speed_pdf"]["cdf"][month])),
            key=lambda i: abs(
                self.historical["wind_speed_pdf"]["cdf"][month][i] - 0.95
            ),
        )
        perc_95th = self.historical["wind_speed_pdf"]["bin"][month][perc_95th_idx]

        if self.meteo_vars["wind_speed"][time] > perc_95th:
            self.meteo_vars["hourly_flags"][time].append(
                "wind speed is in the upper 95th percentile"
            )

    def prec_mag_month_check(self, time):
        # -------------------
        # WATER PRECIPITATION
        # -------------------

        # get all indices
        days_in_month = calendar.monthrange(time.year, time.month)
        end_time = time + datetime.timedelta(days=+days_in_month[1])
        month_slice = pd.date_range(start=time, end=end_time, freq="D")
        month_idx = time.month

        # max and min check
        if (
            sum(self.meteo_vars["prec_water"][month_slice])
            < self.historical["prec_water_min"][month_idx]
        ):
            for hour in month_slice:
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is less than the minimum observed"
                )

        elif (
            sum(self.meteo_vars["prec_water"][month_slice])
            > self.historical["prec_water_max"][month_idx]
        ):

            for hour in month_slice:
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is more than the maximum observed"
                )

        # 5th and 95th percentile check
        perc_5th = (
            self.historical["prec_water_avg"] - 2 * self.historical["prec_water_std"]
        )
        perc_95th = (
            self.historical["prec_water_avg"] + 2 * self.historical["prec_water_std"]
        )

        if sum(self.meteo_vars["prec_water"][month_slice]) < perc_5th:
            for hour in month_slice:
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is in the lower 5th percentile"
                )

        elif sum(self.meteo_vars["prec_water"][month_slice]) > perc_95th:

            for hour in month_slice:
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is in the upper 95th percentile"
                )

    def step_qc(self):

        for i, (time, _) in enumerate(self.meteo_vars["dry_bulb"].iteritems()[:-1]):
            self.db_step_check()
            self.rh_step_check()

    def time_qc(self):
        pass

    def db_step_check(self):
        pass

    def rh_step_check(self):
        pass


def get_pdf(weather, meteo_var):

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

    else:
        raise ValueError(
            "meteo_var must be either 'dry_bulb', 'dew_point', 'rel_hum', 'wind_dir', or 'wind_speed'"
        )

    pdf = [{"bin": [], "freq": [], "cdf": [], "pdf": []} for i in range(12)]

    for i in range(12):
        if i < 9:
            month = f"0{i + 1}"
        else:
            month = f"{i + 1}"
        month_path = weather.folder / f"{weather.wmo}_{meteo_id}_{month}.txt"

        with open(month_path, "r") as fp:
            # skip the five-line header
            for n in range(5):
                next(fp)

            for j, line in enumerate(fp):

                if meteo_id == "DBDP":

                    if j == 0:
                        dew_point = [float(s) for s in line.split()[1:]]
                        dry_bulb = []
                        freq = []

                    else:
                        dry_bulb.append(float(line.split()[0]))
                        freq.append([float(s) for s in line.split()[1:]])

                else:
                    pdf[i]["bin"].append(float(line.split()[0]))
                    pdf[i]["freq"].append(float(line.split()[1]))
                    pdf[i]["cdf"].append(float(line.split()[2]))

        if meteo_id == "DBDP":
            bin_unsorted = []
            freq_unsorted = []
            run_freq_sum = [0]

            for j, (db, dp) in enumerate(itertools.zip_longest(dry_bulb, dew_point)):
                x = math.floor(j / len(dew_point))
                y = j - x * len(dew_point)
                bin_unsorted.append(
                    100 * (17.625 * dp / (243.04 + dp)) / (17.625 * db / (243.04 + db))
                )
                freq_unsorted.append(freq[x][y])
                run_freq_sum.append(run_freq_sum[-1] + freq[x][y])

            # sort bins and corresponding frequencies
            pdf[i]["bin"] = sorted(bin_unsorted)
            pdf[i]["freq"] = [f for _, f in sorted(zip(pdf[i]["bin"], freq_unsorted))]
            pdf[i]["cdf"] = [
                c / run_freq_sum[-1]
                for _, c in sorted(zip(pdf[i]["bin"], run_freq_sum))
            ]

        # calculate probability distribution function
        total_obs = sum(pdf[i]["freq"])
        pdf[i]["pdf"] = [b / total_obs for b in pdf[i]["bin"]]

    return pdf


def get_month_data(excel_sheet, wmo, rad_type):

    # find the corresponding row of data for the weather station
    wmo_col = excel_sheet.col_values(4)[2:]
    row = wmo_col.index(str(wmo))

    # read the values from the corresponding columns in the row of the weather station
    if rad_type == "clr_sky_dir_norm":
        cells = excel_sheet.row_slice(row, 533, 545)

    elif rad_type == "clr_sky_diff_hor":
        cells = excel_sheet.row_slice(row, 545, 557)

    elif rad_type == "dir_norm_taub":
        cells = excel_sheet.row_slice(row, 509, 521)

    elif rad_type == "diff_hor_taud":
        cells = excel_sheet.row_slice(row, 521, 533)

    elif rad_type == "all_sky_glob_avg":
        cells = excel_sheet.row_slice(row, 557, 569)

    elif rad_type == "all_sky_glob_std":
        cells = excel_sheet.row_slice(row, 569, 581)

    elif rad_type == "prec_water_avg":
        cells = excel_sheet.row_slice(row, 206, 218)

    elif rad_type == "prec_water_max":
        cells = excel_sheet.row_slice(row, 219, 231)

    elif rad_type == "prec_water_min":
        cells = excel_sheet.row_slice(row, 232, 244)

    elif rad_type == "prec_water_std":
        cells = excel_sheet.row_slice(row, 245, 257)

    else:
        raise ValueError(
            "rad_type must be either 'clr_sky_dir_norm', 'clr_sky_diff_hor', 'dir_norm_taub', 'diff_hor_taud', "
            "'all_sky_glob_avg', 'all_sky_glob_std', 'prec_water_avg', 'prec_water_max', 'prec_water_min', or "
            "'prec_water_std'"
        )

    # turn N/A values into None and other numerical strings into floats
    month_data = [None if c.value == "N/A" else float(c.value) for c in cells]

    return month_data


def get_solar_angles(weather):

    # get hours in meteorological file
    hours = len(weather.meteo_vars["dry_bulb"])
    latitude = weather.location[0]
    longitude = weather.location[1]

    solar_angles = pd.DataFrame(
        OrderedDict(
            (
                ("declination", np.empty(hours)),
                ("hour_angle", np.empty(hours)),
                ("zenith", np.empty(hours)),
                ("altitude", np.empty(hours)),
                ("azimuth", np.empty(hours)),
                ("night", np.empty(hours)),
            )
        ),
    )

    for i, (time, _) in enumerate(weather.meteo_vars["dry_bulb"].iteritems()):
        # fractional day of the year
        n = (i + 12) / 24
        # equation of time
        B = (n - 1) * 360 / 365
        E = 2.2918 * (
            0.0075
            + 0.1868 * math.cos(math.radians(B))
            - 3.2077 * math.sin(math.radians(B))
            - 1.4615 * math.cos(math.radians(2 * B))
            - 4.089 * math.sin(math.radians(2 * B))
        )

        # solar declination angle
        declin = math.radians(23.45 * math.sin(math.radians(360 * (284 + n) / 365)))
        solar_angles["declination"][i] = declin

        # solar time
        t_s = time.hour + (4 * (weather.time_zone * 15 - longitude) + E) / 60
        # solar hour angle
        hour_angle = math.radians(((t_s % 24) - 12) * 15)
        solar_angles["hour_angle"][i] = hour_angle

        # solar zenith angle
        zenith = math.acos(
            math.cos(latitude) * math.cos(declin) * math.cos(hour_angle)
            + math.sin(latitude) * math.sin(declin)
        )

        if zenith <= math.pi / 2:
            solar_angles["zenith"][i] = zenith
            solar_angles["day"][i] = True
        else:
            solar_angles["zenith"][i] = None
            solar_angles["day"][i] = False

        # if it is day
        if solar_angles["day"][i]:
            altitude = math.asin(math.cos(zenith))
            solar_angles["altitude"][i] = altitude
        else:
            solar_angles["altitude"][i] = None

        # if it is day
        if solar_angles["day"][i]:
            azimuth = math.copysign(
                abs(
                    math.acos(
                        math.copysign(
                            (math.cos(zenith) * math.sin(latitude) - math.sin(declin))
                            / (math.sin(zenith) * math.sin(latitude)),
                            latitude,
                        )
                    )
                ),
                hour_angle,
            )
            solar_angles["azimuth"][i] = azimuth
        else:
            solar_angles["azimuth"][i] = None

    return solar_angles
