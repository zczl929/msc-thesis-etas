# Ridgecrest source catalogue

The frozen submission source is `ridgecrest_usgs_exact.csv`. It was retrieved
from the USGS ComCat FDSN Event Web Service at `2026-07-25T04:06:47Z` using:

```text
https://earthquake.usgs.gov/fdsnws/event/1/query.csv?starttime=2019-07-04T00%3A00%3A00Z&endtime=2019-07-16T03%3A19%3A53.040Z&latitude=35.7695&longitude=-117.5993&maxradiuskm=75&minmagnitude=2.95&orderby=time-asc
```

Frozen source checksums:

- MD5: `9d271affca4842c29237cee47b618bbe`
- SHA-256:
  `e4bec8199e7dfdf83856c5c8fabf7c79d9e7a45126f35a96d311b179af1e1923`

The query contains 1,035 unique events. Its centre uses the rounded USGS query
coordinate `(35.7695, -117.5993)`. Great-circle distances in the preparation
script are recomputed from the catalogue location of the M7.1 mainshock
`ci38457511`, `(35.7695, -117.5993333)`.

`scripts/submission/05_prepare_ridgecrest.R` validates the source checksum,
event IDs, query bounds, mainshock identity and analysis filters before
writing the frozen full and training catalogues.

The older `ridgecrest_raw.csv` is a broader source download retained only for
development history. It is not read by the submission workflow. The older
30-day, M2.5 diagnostic outputs are not submission results.
