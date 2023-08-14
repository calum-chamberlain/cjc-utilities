"""
Tests for the magnitude_inversion functions
"""

import unittest

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

from magnitude_inversion import magnitude_inversion

class TestInversion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        client = Client("GEONET")
        cls.catalog = client.get_events(
            starttime=UTCDateTime(2019, 1, 1), endtime=UTCDateTime(2019, 1, 10),
            minlatitude=-42, maxlatitude=-40, minlongitude=172, 
            maxlongitude=178)

    def test_zero_residuals(self):
        """ 
        Check that giving the same catalog for both args returns similar
        magnitudes. 

        Note that GeoNet use some fucking wierd, almost totally undocumented
        method for computing magnitudes that does not conform to international
        standards and is therefore BOLLOCKS and cannot be compared directly.

        GEONET SUCK!
        """
        cat_out, gamma, station_corrections = magnitude_inversion(
            new_catalog=self.catalog.copy(), 
            callibration_catalog=self.catalog.copy(), magnitude_type="ML",
            frequency_dependent=False)
        mags_out = {
            ev.resource_id.id.split('/')[-1]: ev.preferred_magnitude().mag 
            for ev in cat_out}
        mags_in = {}
        for ev in self.catalog:
            mag = [mag for mag in ev.magnitudes if mag.magnitude_type == "ML"]
            if len(mag) == 0:
                continue
            mags_in.update({
                ev.resource_id.id.split('/')[-1]: mag[0].mag})

        residuals = {
            key: mags_in[key] - mags_out[key] for key in mags_out.keys()}
        self.assertAlmostEqual(
            sum(residuals.values()) / len(residuals), 0.0)


if __name__ == "__main__":
    unittest.main()