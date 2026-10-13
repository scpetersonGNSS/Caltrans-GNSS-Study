# GNSS Static Occupation Data site

A public website for GNSS static occupation data. Each uploaded collector
export becomes its own dataset, with its mean coordinate, precision
statistics, a time series of deviations from the mean, and a horizontal
scatter plot with DRMS circles. Visitors pick an upload from the list or
filter the list by date. The original files are available for download.

## How it works

```
uploads/                 you drop collector exports here
scripts/process_uploads.py   turns each upload into a dataset
data/index.json          list of uploads with summary statistics
data/sessions/<id>.json  plot data for one upload (integer mm offsets)
data/raw/<id>_<file>     the original file, for download
scripts/satellite_geometry.py   predicts satellite counts and DOP per upload
data/geometry/<id>.json  predicted satellites and DOP for one upload
data/nav/                broadcast orbit files downloaded for those predictions
index.html               the website
.github/workflows/update-site.yml   runs the script and publishes the site
```

When a file lands in `uploads/`, GitHub runs the processing script, moves the
file into `data/raw/`, commits the results, and republishes the site. Each
upload is identified by its first timestamp, for example `20260831-072247`.
Uploading the exact same file twice is detected and skipped.

The script auto-detects the Northing, Easting, Elevation, and time columns, so
`.asc`, `.csv`, and `.txt` exports with similar headers work without changes.

After that, `scripts/satellite_geometry.py` predicts, for each tagged upload,
how many satellites of each constellation were above a 10° elevation mask and
the PDOP / HDOP / VDOP for the constellations in that test, every minute, like
Trimble's GNSS Planning charts. It downloads the IGS broadcast orbits for those
days (public, no login) into `data/nav/` and writes `data/geometry/<id>.json`.
The site shows these in the **Satellites and DOP over time** charts, with the
receiver's logged PDOP drawn over the prediction. The elevation mask, the time
zone of the uploaded timestamps, and each site's coordinate zone are set near
the top of that script. Untagged uploads are skipped, since their location is
not known.

## One-time setup

1. Create a free account at github.com if you don't have one.
2. Create a new repository. Make it **Public**. Name it something like
   `gnss-data`.
3. On the new repository page, click **uploading an existing file**. Drag in
   the *contents* of this folder (not the folder itself) and commit.
   The `.github` folder is hidden on some computers. If it doesn't upload,
   create the file by hand: **Add file > Create new file**, name it
   `.github/workflows/update-site.yml`, and paste in the contents.
4. Go to **Settings > Pages**. Under **Build and deployment**, set
   **Source** to **GitHub Actions**.
5. Go to the **Actions** tab and run **Process uploads and publish site**
   (click it, then **Run workflow**). When it finishes, the site is live at
   `https://<your-username>.github.io/<repository-name>/`.

## Adding data (the routine)

1. Download the export from the data collector.
2. In the repository on GitHub, open the `uploads` folder.
3. Click **Add file > Upload files**, drag in the export, and click
   **Commit changes**.
4. Wait a minute or two. The new upload appears at the top of the list.

You can follow progress in the **Actions** tab. A green check means the site
is updated.

## If an upload fails

If a file can't be read (missing columns, unrecognized time format), the run
is marked with a red X and the file stays in `uploads/`. Open the failed run
and expand **Process new uploads** to see the reason. Other files in the same
commit are still processed and published.

To add a new timestamp format, add it to `TIME_FORMATS` near the top of
`scripts/process_uploads.py`.

## Removing an upload

Delete its traces and commit: its entry in `data/index.json`, its file in
`data/sessions/`, its file in `data/raw/`, and its file in `data/geometry/`.

## Previewing on your own computer

```
python scripts/process_uploads.py
python scripts/satellite_geometry.py
python -m http.server 8000
```

Then open http://localhost:8000. Opening `index.html` directly by
double-clicking won't load the data, because browsers block that.

## Notes

* Times are shown exactly as recorded by the collector, with no time zone
  conversion.
* Coordinates are in the units of the export (meters). Statistics use
  population standard deviation. DRMS is √(σN² + σE²).
* GitHub warns about files over 50 MB and refuses files over 100 MB. A
  25-hour, 5-second export is about 1 MB, so this is only a concern for very
  long occupations.
