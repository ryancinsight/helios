# Chapter 25 — Clinical Protocol Compliance

<!-- generated-figure-start -->
![Figure 37.1 — Clinical Protocol Compliance](figures/ch37/fig01_37_clinical_protocol_compliance.svg)
*Figure 37.1 — Clinical Protocol Compliance*
<!-- generated-figure-end -->

Helios validation targets the following published clinical benchmarks.

## TG-119 (IMRT Commissioning)

AAPM TG-119 defines test cases for IMRT commissioning with known dose distributions:
- Simple C-shape, head-and-neck, prostate plans
- Expected point-dose accuracy: ±3% / 3 mm gamma

## TRS-398 Absorbed Dose Protocol

IAEA TRS-398 calibration conditions (6 MV, 10×10 cm² field, 10 cm depth in water):
- Reference dose rate: 1 cGy/MU at calibration geometry
- Helios uses MassAttenuation::water() to reproduce this reference point

## TomoTherapy Workflow

The tomotherapy workflow example achieves:
- Water ROI μ reconstruction error: < 0.2% of μ_water
- DVH mean dose consistent with prescribed fluence
- 3%/2 mm gamma against a 2%-low comparison field: **100%** pass, with the
  6%-low negative control failing as required (see
  [Gamma Index and Plan Verification](planning_gamma.md))

A gamma of the dose against itself is not reported here: it is 100% by
construction and carries no information about the dose engine.

## Further Reading

- [Reference Phantoms](validation_phantoms.md)
- [Gamma Index and Plan Verification](planning_gamma.md)
