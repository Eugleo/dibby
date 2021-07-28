def within_bounds(reference_mass, measured_mass, error_ppm: float = 10):
    return abs(reference_mass - measured_mass) <= err_margin(reference_mass, error_ppm)


def err_margin(reference_mass, error_ppm: float = 10):
    return (error_ppm / 1e6) * reference_mass


def compute_error(reference_mass, measured_mass):
    return 1e6 * abs(measured_mass - reference_mass) / reference_mass
