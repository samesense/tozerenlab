import utils_gmaps, PyMozilla, nose.tools

def test_drivingDirections_assumptions():
    """ Make sure the url form is right and the return html is in the format I think it is in.
    """

    moz_emu = PyMozilla.MozillaEmulator(cacher=None, trycount=0)
    web_page = moz_emu.download('http://maps.google.com/maps?f=d&source=s_d&hl=en&geocode=&saddr=mesquite%2C+tx&daddr=philadelphia%2C+pa&btnG=Get+Directions&output=js')
    distance = web_page.split('timedist')[1].split('mi')[0].split('\\')[-2].split('e')[1].replace(',' ,'')
    nose.tools.assert_equal('1460', distance, 'Malformed URL or HTML return')

def test_drivingDirections_names():
    """ Make sure drivingDirection's inputs are handled correctly for place names.
    """

    ans = utils_gmaps.drivingDistance('mesquite, tx',
                                      'philadelphia, pa')
    nose.tools.assert_equal('1460', ans)

    ans = utils_gmaps.drivingDistance('san diego, ca',
                                      'philadelphia, pa')
    nose.tools.assert_equal('2739', ans)

    ans = utils_gmaps.drivingDistance('hood canal, wa',
                                      'tampa, fl')
    nose.tools.assert_equal('3230', ans)

#def test_drivingDirections_longitude_latitude():
#    """ Make sure drivingDirection's inputs are handled
#        correctly for longitude, latitude inputs.
#    """
#    ans = utils_gmaps.drivingDistance('-96.466825, 32.755857',
#                                      '-75.34, 39.94')
#    nose.tools.assert_equal('1460', ans)
    
