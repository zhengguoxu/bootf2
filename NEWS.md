# bootf2 0.3.0.9001

* Cosmetic changes of source files
* Changes in documentation.

# bootf2 0.3.0

* Change stop to warning for FDA CV rule in `calcf2()`.
* Round CV with precision `digits` then compared to 20%/10% rule.
* Add message for case with `both.TR.85 = TRUE` and `regulation = "FDA"`.
* Fix bug with missing `dp.cv` in `sim.dp()`.

# bootf2 0.2.0

* Added the main function `bootf2()` foe the bootstrap $f_2$.
* Added README.Rmd and vignettes for each main function

# bootf2 0.1.0

* Added 3 working functions:
    - `sim.dp()` for simulation of dissolution profile.
    - `calcf2()` for the calculation of similarity factor $f_2$.
    - `sim.dp.byf2()` for search a suitable model parameters such that the 
      simulated dissolution profile will have f~2~ equal to the predefined 
      target $f_2$. 
