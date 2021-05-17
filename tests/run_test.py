import unittest
from datetime import datetime, timedelta
import cm4twc

from smart.surfacelayer import SMART as SMART_sl
from smart.subsurface import SMART as SMART_ss
from smart.openwater import SMART as SMART_ow


class TestContribution(unittest.TestCase):

    def test_smart(self):

        td = cm4twc.TimeDomain.from_start_end_step(
            start=datetime(2007, 1, 1, 0, 0, 0),
            end=datetime(2007, 3, 1, 0, 0, 0),
            step=timedelta(hours=1)
        )

        sd = cm4twc.LatLonGrid.from_extent_and_resolution(
            latitude_extent=(51, 52),
            latitude_resolution=1,
            longitude_extent=(0, 1),
            longitude_resolution=1
        )

        ds = cm4twc.DataSet(['in/rainfall_flux.nc',
                             'in/potential_water_evapotranspiration_flux.nc'])

        ds2 = cm4twc.DataSet('in/cell_area.nc')

        theta_t = (1.0, '1')
        theta_c = (1.0, '1')
        theta_h = (0.20845296027652363, '1')
        theta_d = (0.24606006380093334, '1')
        theta_s = (0.00012296588050682812, '1')
        theta_z = (105.25734595830215, 'kg m-2')
        theta_sk = (46.81961454361724 * 3600, 's')
        theta_fk = (315.5490902162102 * 3600, 's')
        theta_gk = (1066.7332319333473 * 3600, 's')
        theta_rk = (10.640277777777778 * 3600, 's')

        sl = SMART_sl(
            saving_directory='out',
            timedomain=td,
            spacedomain=sd,
            dataset=ds,
            parameters={
                'theta_t': theta_t,
                'theta_z': theta_z
            }
        )

        ss = SMART_ss(
            saving_directory='out',
            timedomain=td,
            spacedomain=sd,
            dataset=None,
            parameters={
                'theta_c': theta_c,
                'theta_h': theta_h,
                'theta_d': theta_d,
                'theta_s': theta_s,
                'theta_z': theta_z,
                'theta_sk': theta_sk,
                'theta_fk': theta_fk,
                'theta_gk': theta_gk
            }
        )

        ow = SMART_ow(
            saving_directory='out',
            timedomain=td,
            spacedomain=sd,
            dataset=ds2,
            parameters={
                'theta_rk': theta_rk
            },
            records={
                'outgoing_water_volume_transport_along_river_channel': {
                    timedelta(days=1): ['mean']
                }
            }
        )

        model = cm4twc.Model(
            identifier='test-smart',
            config_directory='out',
            saving_directory='out',
            surfacelayer=sl,
            subsurface=ss,
            openwater=ow
        )

        model.simulate()

        from_file = cm4twc.DataSet(
            'in/outgoing_water_volume_transport_along_river_channel.nc'
        )

        from_model = cm4twc.DataSet(
            'out/test-smart_openwater_run_records_daily.nc'
        )

        var_name = 'outgoing_water_volume_transport_along_river_channel'

        self.assertTrue(
            from_file[var_name].field.equals(from_model[var_name].field,
                                             verbose=3)
        )


if __name__ == '__main__':
    test_loader = unittest.TestLoader()
    test_suite = unittest.TestSuite()

    test_suite.addTests(
        test_loader.loadTestsFromTestCase(TestContribution))

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(test_suite)
