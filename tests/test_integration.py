from integration import experiment

def test_shared_model_contract():
    result = experiment()
    assert result["rows"][-1]["buckling_relative_error"] < 1e-5
    assert abs(result["rows"][-1]["dimensionless_product"]/result["product_reference"]-1) < 0.001
    assert result["rows"][-1]["reduced_compliance_relative_error"] < result["rows"][-1]["full_compliance_relative_error"]
