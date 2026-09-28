# problems to fix

* Electron mode: peak heights are multiplicity × Debye–Waller only — no
  electron scattering factor f_e(s), so relative intensities are not
  quantitative (low-g reflections should dominate once f_e(s) is added; it can
  be computed from `atomic_form_factor` by the Mott–Bethe relation)
