import pytest
import rasterio

@pytest.fixture
def original_data():
    out_data = rasterio.open('fixtures/coweeta.bil')
    return out_data.read(1)


@pytest.fixture
def new_data():
    out_data = rasterio.open('results/coweeta_output_SLOPE.bil')
    return out_data.read(1)


class TestingLSD():
    def test_slope(self, original_data, new_data):
        assert (original_data == new_data).all()
