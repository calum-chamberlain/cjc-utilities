"""
Convert Obspy inventories to NLL Station format:

GTSRCE {STATION:5s} LATLON {LAT:7.4f}  {LON:7.4f} {DEPTH:.1f} {ELEVATION:.3f}

"""

from obspy import read_inventory
from obspy.core.inventory import Inventory


def inv_to_nll(
    inv: Inventory
) -> str:
    """
    Convert an obspy inventory to a nll formatted string of station locations.
    """
    lats, lons, elevs, depths = dict(), dict(), dict(), dict()
    for net in inv:
        for sta in net:
            for chan in sta:
                depth = chan.depth / 1000.0
                elev = chan.elevation / 1000.0
                lat, lon = chan.latitude, chan.longitude
                
                # Check values are consistent
                assert depths.get(sta.code) in [None, depth]
                assert elevs.get(sta.code) in [None, elev]
                assert lats.get(sta.code) in [None, lat]
                assert lons.get(sta.code) in [None, lon]

                depths[sta.code] = depth
                elevs[sta.code] = elev
                lats[sta.code] = lat
                lons[sta.code] = lon

    lines = [
        f"GTSRCE {sta.code:5s} LATLON {lats[sta.code]:7.4f}  "
        f"{lons[sta.code]:7.4f} {depths[sta.code]:02.1f} "
        f"{elevs[sta.code]:.3f}"
        for net in inv for sta in net]

    return "\n".join(lines)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert an obspy readable Inventory to NonLinLoc format")

    parser.add_argument(
        "-i", "--infile", type=str, required=True)
    parser.add_argument(
        "-o", "--outfile", type=str, default="stations.nll")

    args = parser.parse_args()

    inv = read_inventory(args.infile)
    nll_lines = inv_to_nll(inv)
    with open(args.outfile, "w") as f:
        f.write(nll_lines)


if __name__ == "__main__":
    main()
