from blimpy import Waterfall, LoTSSWaterfall

fb = LoTSSWaterfall('/Users/ascaife/SRC/GITHUB/blimpy/test_data/lotss_test.fits')
fb.info()

data = fb.data
