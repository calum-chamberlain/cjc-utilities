"""
Utility functions for focal mechanism work.

"""

from math import (
    sin, cos, tan, acos, radians, degrees, pi, asin, atan2, sqrt, atan)
import numpy as np

from obspy.core.event import (
    NodalPlane, PrincipalAxes, Axis, MomentTensor, Tensor)


def aux_plane(nodal_plane):
    """
    Calculate the strike, dip and rake of the auxiliary nodal plane.

    :type nodal_plane: `obspy.core.event.NodalPlane`
    :param nodal_plane: The nodal-plane to calculate from.

    :return: `obspy.core.event.source.NodalPlane` auxiliary nodal plane.
    """
    # Following equations 18, 19, 16 and 17 P. 228 Stein & Wysession
    # slip angle = rake
    (strike_in, dip_in, rake_in) = _nodal_plane_radian(nodal_plane)
    dip = acos(sin(rake_in) * sin(dip_in))
    # sin_rake = math.cos(dip_rad_in) / math.sin(dip_rad_out)
    cos_rake = -1 * (sin(dip_in) * cos(rake_in) / sin(dip))
    rake = acos(cos_rake)

    sin_strike1_strike2 = cos(rake_in) / sin(dip)
    try:
        cos_strike1_strike2 = -1 / (tan(dip_in) * tan(dip))
    except ZeroDivisionError:
        cos_strike1_strike2 = -1
    strike1_strike2 = acos(cos_strike1_strike2)

    # Check the quadrant of strike
    if sin_strike1_strike2 < 0 and cos_strike1_strike2 < 0:
        strike1_strike2 *= -1
    elif sin_strike1_strike2 < 0 and cos_strike1_strike2 > 0:
        strike1_strike2 *= -1
    strike_rad_out = strike_in - strike1_strike2
    if rake_in < 0:
        rake = -(2 * pi - rake)
    strike = degrees(strike_rad_out)
    dip = degrees(dip)
    rake = degrees(rake)
    if 90 < dip < 180:
        strike += 180
        dip = 180 - dip
        rake = 360 - rake
    return NodalPlane(strike=round(strike % 360, 4), dip=round(dip % 360, 4),
                      rake=round(rake % 360, 4))


def np_to_principal_axes(nodal_plane):
    """
    Calculate the P, T and B axes of a focal mechanism given one focal plane.

    :type nodal_plane: `obspy.core.event.NodalPlane`
    :param nodal_plane: The nodal-plane to calculate from.

    :return: `obspy.core.event.source.PrincipalAxes` principal axes.

    .. Note::
        Length of axes is not calculated.
    """
    # Following Stein & Wysession equations 8, 9, 10 P. 228
    normal = normal_vector(nodal_plane=nodal_plane)
    slip = slip_vector(nodal_plane=nodal_plane)
    null = np.cross(normal, slip)
    tension = normal + slip
    # Looks like P is not great, but t is
    pressure = normal - slip
    principal_axes = PrincipalAxes(
        t_axis=trend_plunge(tension, "tension"),
        p_axis=trend_plunge(pressure, "pressure"),
        n_axis=trend_plunge(null, "null"))
    return principal_axes


def trend_plunge(vector, name=None):
    """
    Calculate the trend and plunge of a given vector

    :type vector: np.ndarray
    :param vector: Vector to calculate for, assumes [x, y, z]

    :return: `obspy.core.event.Axis`
    """
    if name:
        print("Calculating for {0}".format(name))
    magnitude = sqrt(np.sum(vector ** 2))
    print("\t{0}".format(vector))
    # Normalise to unit vector
    # vector /= magnitude
    plunge = degrees(asin(vector[2] / magnitude))
    plunge = plunge + 90 % 180 - 90
    azi = degrees(atan2(vector[1], vector[0]))
    if vector[1] < 0:
        azi += 180
    azi = azi % 360
    print("\tAzimuth: {0}\tPlunge: {1}\tLength:{2}".format(
        azi, plunge, magnitude))
    return Axis(plunge=plunge, azimuth=azi, length=magnitude)


def normal_vector(nodal_plane):
    """
    Calculate the normal to a nodal plane.

    :type nodal_plane: `obspy.core.event.NodalPlane`
    :param nodal_plane: The nodal-plane to calculate from.

    :returns: 1D np.ndarray of axes co-ordinates [N, E, Up]
    """
    (strike, dip, _) = _nodal_plane_radian(nodal_plane)
    normal = [-sin(dip) * sin(strike), sin(dip) * cos(strike), cos(dip)]
    return np.array(normal)


def slip_vector(nodal_plane):
    """
    Calculate the slip vector of a given nodal plane.

    :type nodal_plane: `obspy.core.event.NodalPlane`
    :param nodal_plane: The nodal-plane to calculate from.

    :returns: 1D np.ndarray of axes co-ordinates
    """
    (strike, dip, rake) = _nodal_plane_radian(nodal_plane)
    slip = [
        (cos(rake) * cos(strike)) + (sin(rake) * cos(dip) * sin(strike)),
        (cos(rake) * sin(strike)) - (sin(rake) * cos(dip) * cos(strike)),
        sin(rake) * sin(dip)]
    return np.array(slip)


def mt_from_dc(nodal_plane, scalar_moment=1):
    """
    Calculate a moment tensor from a nodal_plane - normalised

    :type nodal_plane: `obspy.core.event.NodalPlane`
    :param nodal_plane: The nodal-plane to calculate from.

    :returns: `obspy.core.event.MomentTensor`
    """
    n = normal_vector(nodal_plane)
    d = slip_vector(nodal_plane)
    mt = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            mt[i][j] = scalar_moment * (n[i] * d[j] + n[j] * d[i])
    mt = MomentTensor(scalar_moment=scalar_moment, tensor=Tensor(
        m_rr=mt[0][0], m_tt=mt[1][1], m_pp=mt[2][2], m_rt=mt[0][1],
        m_rp=mt[0][2], m_tp=mt[1][2]))
    return mt


def _nodal_plane_radian(nodal_plane):
    """
    Get strike, dip and rake in radians

    :type nodal_plane: `obspy.core.event.NodalPlane`
    :param nodal_plane: The nodal-plane to calculate from.

    :returns: tuple of strike, dip, rake
    """
    strike, dip, rake = (nodal_plane.strike, nodal_plane.dip, nodal_plane.rake)
    # strike = 360 - strike  # Strike is clockwise from N
    # rake = 90 - rake
    if dip > 90:
        dip = 180 - dip
        strike += + 180
    #print("Using S/D/R: {0}/{1}/{2}".format(strike, dip, rake))
    return radians(strike), radians(dip), radians(rake)


def cosd(angle):
    return cos(radians(angle))


def sind(angle):
    return sin(radians(angle))


def tand(angle):
    return tan(radians(angle))


def atan2d(value1, value2):
    return degrees(atan2(value1, value2))


def asind(value):
    return degrees(asin(value))


def sincosd(angle):
    return sin(radians(angle)), cos(radians(angle))


def hypot(adj, opp):
    return sqrt(adj ** 2 + opp ** 2)

def atand(value):
    return degrees(atan(value))

def meca_dc2axe(nodal_plane_1, nodal_plane_2):
    """
    From GMT mecautils.c
    
    From FORTRAN routines of Anne Deschamps :
    compute azimuth and plungement of P-T axis
    from nodal plane strikes, dips and rakes.
    """
    M_SQRT2 = sqrt(2)
    cd1 = cosd(nodal_plane_1.dip) * M_SQRT2
    sd1 = sind(nodal_plane_1.dip) * M_SQRT2
    cd2 = cosd(nodal_plane_2.dip) * M_SQRT2
    sd2 = sind(nodal_plane_2.dip) * M_SQRT2
    cp1 = - cosd(nodal_plane_1.strike) * sd1
    sp1 = sind(nodal_plane_1.strike) * sd1
    cp2 = - cosd(nodal_plane_2.strike) * sd2
    sp2 = sind(nodal_plane_2.strike) * sd2
    amz = - (cd1 + cd2)
    amx = - (sp1 + sp2)
    amy = cp1 + cp2
    dx = atan2d (hypot(amx, amy), amz) - 90.0
    px = atan2d (amy, -amx)
    if px < 0.0:
        px += 360.0
    if dx < 0.0001:
        if px > 90.0 and px < 180.0:
            px += 180.0
        if px >= 180.0 and px < 270.0:
            px -= 180.0
    amz = cd1 - cd2
    amx = sp1 - sp2
    amy = - cp1 + cp2
    dy = atan2d (hypot(amx, amy), -abs(amz)) - 90.0
    py = atan2d (amy, -amx)
    if amz > 0.0:
        py -= 180.0
    if py < 0.0:
        py += 360.0
    if dy < 0.0001:
        if py > 90.0 and py < 180.0:
            py += 180.0
        if py >= 180.0 and py < 270.0:
            py -= 180.0
    principal_axes = PrincipalAxes(
        t_axis=Axis(), p_axis=Axis(), n_axis=Axis())
    if nodal_plane_1.rake > 0.0:
        principal_axes.p_axis.plunge = dy
        principal_axes.p_axis.azimuth = py
        principal_axes.t_axis.plunge = dx
        principal_axes.t_axis.azimuth = px
    else:
        principal_axes.p_axis.plunge = dx
        principal_axes.p_axis.azimuth = px
        principal_axes.t_axis.plunge = dy
        principal_axes.t_axis.azimuth = py
    principal_axes.n_axis.azimuth = null_axis_strike(
        principal_axes.t_axis.azimuth, principal_axes.t_axis.plunge,
        principal_axes.p_axis.azimuth, principal_axes.p_axis.plunge)
    principal_axes.n_axis.plunge = null_axis_dip(
        principal_axes.t_axis.azimuth, principal_axes.t_axis.plunge,
        principal_axes.p_axis.azimuth, principal_axes.p_axis.plunge)
    return principal_axes


def null_axis_dip(str1, dip1, str2, dip2):
    """
    from gmt

    compute null axis dip when strike and dip are given
	for each nodal plane.  Angles are in degrees.

	   Genevieve Patau
	"""
    den = asind (sind (dip1) * sind (dip2) * sind (str1 - str2))
    if den < 0.:
        den = -den
    return den


def null_axis_strike(str1, dip1, str2, dip2):
    """
    from gmt

    Compute null axis strike when strike and dip are given
    for each nodal plane.   Angles are in degrees.

    Genevieve Patau
    """

    (sd1, cd1) = sincosd(dip1)
    (sd2, cd2) = sincosd(dip2)
    (ss1, cs1) = sincosd(str1)
    (ss2, cs2) = sincosd(str2)

    cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1
    sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1
    if sind(str1 - str2) < 0.0:
        cosphn = -cosphn
        sinphn = -sinphn
    phn = atan2d(sinphn, cosphn)
    if phn < 0.0:
        phn += 360.0
    return phn


if __name__ == "__main__":
    import doctest
    doctest.testmod()
