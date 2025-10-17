import mocca

def test_symbol():
    assert mocca.symbol(1) == 'H'
    assert mocca.symbol(20) == 'Ca'
