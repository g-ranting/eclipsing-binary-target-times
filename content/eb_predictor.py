from datetime import datetime
from pytz import timezone


import traitlets as tr
from astroplan import (
    EclipsingSystem,
    FixedTarget,
    AtNightConstraint,
    AltitudeConstraint,
    Observer,
    is_event_observable,
)

from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.table import join
import numpy as np


class EclipsingBinaryPredictor(tr.HasTraits):
    """
    All of the inputs for an eclipsing binary query in one place.
    """
    start_time_str = tr.Instance(datetime).tag(sync=True) #tr.Unicode(default_value="2024-07-30 12:00").tag(sync=True)
    days_in_future = tr.Integer(default_value=30).tag(sync=True)
    minimum_altitude = tr.Float(default_value=30).tag(sync=True)
    latitude = tr.Float(default_value=46.866767).tag(sync=True)
    longitude = tr.Float(default_value=-96.453299).tag(sync=True)
    elevation = tr.Float(default_value=311).tag(sync=True)

    def init(self):
        #super().__init__()
        # Observing time
        self.start_time = Time(self.start_time_str, scale="utc")  # observation time, Requires manual input
        self.end_time = self.start_time + self.days_in_future * u.day
        self.constraints = [
            AtNightConstraint.twilight_nautical(),
            AltitudeConstraint(min=self.minimum_altitude * u.deg),
        ]  # Altitude contraint is possible manual input

        # Observer location

        location = EarthLocation.from_geodetic(self.longitude, self.latitude, self.elevation * u.m)
        self.observer = Observer(
            location=location,
        )

    def can_we_observe(self, target_object, target_obs, n_transits):
        # Target Data

        midtransit_times = target_obs.next_primary_eclipse_time(
            self.start_time, n_eclipses=n_transits
        )

        is_observable = is_event_observable(
            self.constraints, self.observer, target_object, times=midtransit_times
        )

        return midtransit_times[is_observable.flatten()]

    def main(self):
        calc_new_name = Table.read("rich_table.ecsv")  # reads file and creates table
        can_ever_see = 90 - self.observer.location.lat.deg + calc_new_name["coordinates"].dec.deg
        can_ever_see = can_ever_see > self.minimum_altitude #.to_value("deg")
        calc_new_name = calc_new_name[can_ever_see]

        all_eclipses = []  # results in a full list of objects, regardless of whether we can observe it

        for stars in calc_new_name:
            target_object = FixedTarget(
                coord=stars["coordinates"]
            )
            target_obs = EclipsingSystem(
                primary_eclipse_time=stars["Epoch"],
                orbital_period=stars["Period"] * u.day,
                duration=stars["EclipseDuration"] * u.day,
            )
            #print(self.end_time - self.start_time)
            n_transits = (self.end_time - self.start_time)  / stars["Period"]
            n_transits = int(n_transits.jd) #+ 1
            if n_transits >0:
                all_eclipses.append(
                    self.can_we_observe(
                        target_object,
                        target_obs,
                        n_transits,
                    )
                )
            else:
                all_eclipses.append([])
        del calc_new_name["coordinates"]
        concise_eclipses = [
            stars for stars in all_eclipses if len(stars) > 0
        ]  # only returns objects with a duration.

        filtering = [len(stars) > 0 for stars in all_eclipses]

        calc_new_name = calc_new_name[
            filtering
        ]  # clips the table to fit only the objects with a duration

        names = []
        transits = []
        for star, eclipses in zip(calc_new_name, concise_eclipses):
            for eclipse in eclipses:
                names.append(star["Name"])
                transits.append(eclipse)

        eclipsetable = Table(data=[names, transits], names=["Name", "MidtransitTimes"])

        final_table = join(calc_new_name, eclipsetable, join_type="right")

        final_table["EclipseStart"] = (
            final_table["MidtransitTimes"] - final_table["EclipseDuration"] / 2
        )
        final_table["EclipseEnd"] = (
            final_table["MidtransitTimes"] + final_table["EclipseDuration"] / 2
        )

        # 👇👇👇👇👇 CDT is hard-coded here 👇👇👇👇👇
        final_table["EclipseDate"] = final_table["EclipseStart"] - 5 * u.hour

        final_table["EclipseDate"] = final_table["EclipseDate"].to_value(
            "iso", subfmt="date"
        )

        full_eclipse = []
        half1st_eclipse = []
        half2nd_eclipse = []

        for row in final_table:
            target_coord = SkyCoord(
                ra=row["RA2000"] * u.deg, dec=row["Declination2000"] * u.deg
            )
            target_object = FixedTarget(coord=target_coord)
            # 👇👇👇👇  add a padding of some sort here  👇👇👇👇👇👇
            # yeah, hear
            times_full = np.array([Time([row["EclipseStart"], row["EclipseEnd"]])])

            times_first_half = np.array(
                [Time([row["EclipseStart"], row["MidtransitTimes"]])]
            )
            times_second_half = np.array(
                [Time([row["MidtransitTimes"], row["EclipseEnd"]])]
            )
            full_eclipse.append(
                full_vis := is_event_observable(
                    self.constraints, self.observer, target_object, times_ingress_egress=times_full
                )[0, 0]
            )

            if full_vis:
                half1st_eclipse.append(True)
                half2nd_eclipse.append(True)
                continue

            half1st_eclipse.append(
                is_event_observable(
                    self.constraints,
                    self.observer,
                    target_object,
                    times_ingress_egress=times_first_half,
                )[0, 0]
            )
            half2nd_eclipse.append(
                is_event_observable(
                    self.constraints,
                    self.observer,
                    target_object,
                    times_ingress_egress=times_second_half,
                )[0, 0]
            )

        final_table["FullEclipse"] = full_eclipse
        final_table["FirstHalf"] = half1st_eclipse
        final_table["SecondHalf"] = half2nd_eclipse

        keep=final_table["FullEclipse"] | final_table["FirstHalf"] | final_table["SecondHalf"]
        

        self.predictions = final_table[keep]
        #final_table.write("final_table.csv", format="csv", overwrite=True)


if __name__ == "__main__":
    ebi = EclipsingSystemInputs()
    ebi.main()
