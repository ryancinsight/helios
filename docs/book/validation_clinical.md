# Chapter 25 — Clinical Protocol Compliance

Helios validation targets the following published clinical benchmarks.

## TG-119 (IMRT Commissioning)

AAPM TG-119 defines test cases for IMRT commissioning with known dose distributions:
- Simple C-shape, head-and-neck, prostate plans
- Expected point-dose accuracy: ±3% / 3 mm gamma

## TRS-398 Absorbed Dose Protocol

IAEA TRS-398 calibration conditions (6 MV, 10×10 cm² field, 10 cm depth in water):
- Reference dose rate: 1 cGy/MU at calibration geometry
- Helios uses MassAttenuation::water() to reproduce this reference point

## TomoTherapy Self-Gamma

The tomotherapy workflow example achieves:
- Water ROI μ reconstruction error: < 0.2% of μ_water
- DVH mean dose consistent with prescribed fluence
- 3%/2 mm self-gamma pass rate: **100%**

## Further Reading

- [Reference Phantoms](validation_phantoms.md)
- [Gamma Index and Plan Verification](planning_gamma.md)