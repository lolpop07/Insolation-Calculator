import math

def calculate_insolation(day, coordinates, col_azimuth, col_tilt, solar_time, reflectance=0.2):
        """
    :param coordinates: tuple of (lattitude, longtitude, hight_above_sea_level) coordinates of the location where the plant is growing
    :param col_azimuth: float - angle of collector azimuth in degrees
    :param col_tilt: float - collector tilt angle in degrees
    :param solar_time: int - time in hours
    :param reflectance: reflectance value of the environment - default=0.2, None=0, snow=0.8
    :return: solar insolation intensity in W/hour (beam insolation, diffuse insolation, reflected insolation, total insolation)
    >>> calculate_insolation(141,(33.7,84.4,75.5), 20, 52, 12, reflectance=0.2)
    (696.5476928599205, 88.13843553845405, 37.89789582578626, 822.5840242241608)
        """
        L = coordinates[0]  # latitude in degrees
        phic = col_azimuth  # collector altitude in degrees
        Z = col_tilt  # collector tilt in degrees
        ST = solar_time  # time in hours
        theta = 23.45 * math.sin(math.radians(360 / 365 * (day - 81)))  # Declination of the sun in degrees
        H = 15 * (12 - ST)  # hour angle in degrees
        cosL = math.cos(math.radians(L))
        costheta = math.cos(math.radians(theta))
        cosH = math.cos(math.radians(H))
        sinL = math.sin(math.radians(L))
        sintheta = math.sin(math.radians(theta))
        beta = math.degrees(math.asin(cosL * costheta * cosH + sinL * sintheta))  # Altitude angle in degrees
        sinbeta = math.sin(math.radians(beta))
        sinH = math.sin(math.radians(H))
        cosbeta = math.cos(math.radians(beta))
        m = ((708 * sinbeta) ** 2 + 1417) ** 0.5 - 708 * sinbeta  # air mass ratio
        A = 1160 + 75 * math.sin(math.radians(360 / 365 * (day - 275)))  # Apparent extraterrestrial flux
        k = 0.174 + 0.035 * math.sin(math.radians(360 / 365 * (day - 100)))  # Optical depth
        Ib = A * math.exp(-k * m)  #
        tantheta = math.tan(math.radians(theta))
        tanL = math.tan(math.radians(L))
        if cosH > tantheta / tanL:
            phis = A * math.sin(costheta * sinH / cosbeta)
        else:
            phis = 180 - A * math.sin(costheta * sinH / cosbeta)
        sinz = math.sin(math.radians(Z))
        cosz = math.cos(math.radians(Z))
        Ibc = Ib * (math.cos(math.radians(phis - phic)) * cosbeta * sinz + sinbeta * cosz)  # direct beam radiation
        C = 0.095 + 0.04 * math.sin(math.radians(360 / 365 * (day - 100)))  # sky diffuse factor
        Idc = Ib*C*(1 + cosz) / 2  # diffuse radiation
        Irc = reflectance * Ib * (C + sinbeta)*(1 - cosz) / 2  # reflected radiation
        Ic = Ibc + Idc + Irc  # Total insolation
        return Ibc,Idc,Irc,Ic

res = calculate_insolation(141,(33.7,84.4,75.5), 20, 52, 12, reflectance=0.2)
print(res)

# to do: change formula, so it takes into account the hight above sea level. Ib formula
# to do: change formula, so it takes into account the weather at the plant site Ib and Idc formulas (k and m)
