from datetime import datetime
import pathlib

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.era5 import ERA5

if __name__ == '__main__':
    startdate=datetime(2023, 12, 1)
    rnday = 10
    hgrid=Hgrid.open('./input_files/hgrid.gr3',crs='EPSG:4326')
    bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')
    outdir = pathlib.Path('./atmospheric/ERA5/')#


    er=ERA5()
    er.write(outdir=outdir, start_date=startdate, rnday=rnday, air=True, rad=True, prc=True, bbox=bbox, overwrite=True)
