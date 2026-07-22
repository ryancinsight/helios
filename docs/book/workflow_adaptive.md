# Chapter 19 — Adaptive Radiotherapy with MVCT

Adaptive radiotherapy (ART) modifies the treatment plan based on daily
MVCT imaging to account for anatomical changes during the course.

## Workflow

`
Daily MVCT → Rigid/deformable registration to planning CT
                            │
                    Contour propagation
                            │
                    Dose recalculation on current anatomy
                            │
                    Adaptation trigger? → Re-plan if needed
`

## Online vs Offline ART

| Mode | Timing | Latency |
|---|---|---|
| Offline | Between fractions | Days |
| Online | Same-day, pre-treatment | Hours |
| Real-time | During delivery | Seconds |

## Atlas Integration

Helios uses 
itk-registration for deformable registration and

itk-io for DICOM RT dose/structure set I/O.

## Further Reading

- [LINAC-Based Delivery](workflow_linac.md)
- [TomoTherapy Workflow](workflow_tomotherapy.md)