from datetime import datetime
import numpy as np

from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides

if __name__ == '__main__':
    start_date = datetime(2024, 1, 3)
    rnday = 14
    bctypes = [[3, 3, 0, 0], [3, 3, 0, 0]]
    constituents = 'major'
    database = 'tpxo'
    earth_tidal_potential = True
    sthconst = [np.nan, np.nan, 0]
    tobc = [0.5, 0.5, 1]
    sobc = [0.5, 0.5, 1]
    outdir = './'
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    bctides = Bctides(
        hgrid = hgrid,
        flags = bctypes,
        constituents = constituents,
        database = database,
        add_earth_tidal = earth_tidal_potential,
        sthconst = sthconst,
        tobc = tobc,
        sobc = sobc,
    )

    bctides.write(
        './tides/tpxo/', 
        start_date=start_date, 
        rnday=rnday, 
        overwrite=True,
    )
