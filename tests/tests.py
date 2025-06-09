from internal_ions.util import converter


def test_fragment_parse():
    c = converter.JSONConverter
    assert c._parse_fragment_code('c:zdot@2:4(+1)[]') == (2, 4, 'c', 'zdot', 1, '')
    assert c._parse_fragment_code('cdot:z@2:4(+1)[H2O]') == (2, 4, 'cdot', 'z', 1, 'H2O')
