"""Compact-object photometric signatures implemented by Miletos."""


def compute_photometric_signatures(*args, **kwargs):
    """Calculate beaming, ellipsoidal, and self-lensing signal amplitudes."""

    from miletos.signatures import compute_photometric_signatures as calculate

    return calculate(*args, **kwargs)