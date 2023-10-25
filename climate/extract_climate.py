import os.path
import time

import pandas as pd
from emodpy_malaria.weather import (generate_weather, weather_to_csv, WeatherVariable, 
                                    csv_to_weather)

# ---| Request weather files |---
start_yr = 2010
length = 10
ds_name = "BK"
subfolder = f"{ds_name}_{length}yr"
coordcsv = f"simulation_inputs/coordinates/{ds_name}.csv"
extractdir = 'simulation_inputs/tmp/'
outrootdir = os.path.join('simulation_inputs', ds_name, subfolder)
shift_time = False
shift_temp = True

if not os.path.exists(extractdir):
    os.makedirs(extractdir)

weather_dir = extractdir
startdate = start_yr * 1000 + 1
enddate = (start_yr + length - 1) * 1000 + 365

wr = generate_weather(platform="Calculon",
                      site_file=os.path.join(coordcsv),
                      start_date=startdate,
                      end_date=enddate,
                      node_column="nodes",
                      id_reference=ds_name,
                      local_dir=weather_dir)
time.sleep(10)

df, wa = weather_to_csv(weather_dir)
weather_columns = {
    WeatherVariable.AIR_TEMPERATURE: "airtemp",
    WeatherVariable.RELATIVE_HUMIDITY: "humidity",
    WeatherVariable.RAINFALL: "rainfall"
}
weather_filenames = {
    WeatherVariable.AIR_TEMPERATURE: "air_temperature_daily.bin",
    WeatherVariable.RELATIVE_HUMIDITY: "relative_humidity_daily.bin",
    WeatherVariable.RAINFALL: "rainfall_daily.bin"
}

# Remove extra day in 2016
df = df[df.steps != 1096]
df.steps = [x if x < 1096 else x - 1 for x in df.steps]

# Shifting rainfall
if shift_time:
    for offset in [0, 30, 60, 120, -30, -60, -120]:  # fill in the time shift you want here...
        df1 = df.copy()
        df1['new_steps'] = df1.steps + offset
        df1.new_steps = df1.new_steps % (length * 365)
        cond = df1.new_steps <= 0
        df1.new_steps[cond] = df1.new_steps[cond] + (length * 365)
        df1.steps = df1.new_steps
        df1.drop('new_steps', axis=1, inplace=True)
        df1.sort_values('steps', inplace=True)
        df1.reset_index(inplace=True, drop=True)

        # Set constant air temperature to the mean
        airtempMean = df1["airtemp"].mean()
        df1.loc[:, "airtemp"] = airtempMean
        print(df1)

        outdir = os.path.join(outrootdir, f"timeshift_{offset}")
        csv_to_weather(df1, attributes=wa, weather_columns=weather_columns,
                       weather_dir=outdir,
                       weather_file_names=weather_filenames)


if shift_temp:
    for tempshift in [0]:
        df1 = df.copy()

        # Set constant air temperature to the mean
        airtempMean = df1["airtemp"].mean()
        df1.loc[:, "airtemp"] = airtempMean + tempshift

        outdir = os.path.join(outrootdir, f"tempshift_{tempshift}")
        csv_to_weather(df1, attributes=wa, weather_columns=weather_columns,
                       weather_dir=outdir,
                       weather_file_names=weather_filenames)
