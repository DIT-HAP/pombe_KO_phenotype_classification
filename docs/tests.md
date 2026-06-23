Test Suite — Phenotype Classification Pipeline
==============================================

This document describes the test suite in ``tests/``. The tests cover the
three critical layers of the pipeline:

1. **Signal engine** — :func:`growth_signals.classify_growth`
2. **Segment filtering** — :func:`growth_signals.filter_primary_segments`
3. **Single/Multiple classification** — :func:`03_group_genes.classify_phenotype_count`

All 77 tests pass:

.. code-block:: text

    $ mamba run -n bioinformatics python -m pytest tests/ -v
    77 passed in 0.21s

---

Test file: ``tests/test_growth_signals.py``
--------------------------------------------

Tests for :func:`growth_signals.classify_growth` as a pure function —
given a description text, what ``(category, tier)`` pair is returned?

**No signal → WT-like (tier 5)**

+---------------------------------------------------------------------+--------------+--------+
| Description                                                         | Category     | Tier   |
+=====================================================================+==============+========+
| ``VIABLE WT cells``                                                 | WT-like      | 5      |
+---------------------------------------------------------------------+--------------+--------+
| ``" "`` (empty string)                                              | WT-like      | 5      |
+---------------------------------------------------------------------+--------------+--------+
| ``VIABLE misshapen cells at 25,32``                                 | WT-like      | 5      |
+---------------------------------------------------------------------+--------------+--------+
| ``VIABLE slightly long cells``                                      | WT-like      | 5      |
+---------------------------------------------------------------------+--------------+--------+
| ``VIABLE cells at 25,32``                                           | WT-like      | 5      |
+---------------------------------------------------------------------+--------------+--------+

**Single signals**

+---------------------------------------------------------------------+--------------------+--------+
| Description                                                         | Category           | Tier   |
+=====================================================================+====================+========+
| ``ESSENTIAL spores``                                                | spores             | 1      |
+---------------------------------------------------------------------+--------------------+--------+
| ``ESSENTIAL spores misshapen cells``                                | spores             | 1      |
+---------------------------------------------------------------------+--------------------+--------+
| ``ESSENTIAL germinated spores``                                     | germinated         | 2      |
+---------------------------------------------------------------------+--------------------+--------+
| ``ESSENTIAL germinated spores long``                                | germinated         | 2      |
+---------------------------------------------------------------------+--------------------+--------+
| ``VIABLE germination at 25``                                        | germinated         | 2      |
+---------------------------------------------------------------------+--------------------+--------+
| ``VIABLE WT cells, some germination long``                          | germinated         | 2      |
+---------------------------------------------------------------------+--------------------+--------+
| ``ESSENTIAL germinated spores divide``                              | germinated and     | 2      |
|                                                                     | divided            |        |
+---------------------------------------------------------------------+--------------------+--------+
| ``ESSENTIAL germinated spores division``                            | germinated and     | 2      |
|                                                                     | divided            |        |
+---------------------------------------------------------------------+--------------------+--------+
| ``germinated spores divides once``                                  | germinated and     | 2      |
|                                                                     | divided            |        |
+---------------------------------------------------------------------+--------------------+--------+
| ``ESSENTIAL microcolonies``                                         | microcolonies      | 3      |
+---------------------------------------------------------------------+--------------------+--------+
| ``ESSENTIAL microcolonies misshapen cells``                         | microcolonies      | 3      |
+---------------------------------------------------------------------+--------------------+--------+
| ``VIABLE small colonies``                                           | small colonies     | 4      |
+---------------------------------------------------------------------+--------------------+--------+
| ``VIABLE small colonies long cells``                                | small colonies     | 4      |
+---------------------------------------------------------------------+--------------------+--------+
| ``VIABLE very small colonies``                                      | very small         | 4      |
|                                                                     | colonies           |        |
+---------------------------------------------------------------------+--------------------+--------+

**Combined / mixed signals**

+---------------------------------------------------------------------+----------------------------+--------+
| Description                                                         | Category                   | Tier   |
+=====================================================================+============================+========+
| ``spores, germinated spores``                                       | germinated, spores         | 2      |
+---------------------------------------------------------------------+----------------------------+--------+
| ``spores, germinated spores, microcolonies``                        | germinated, microcolonies, | 3      |
|                                                                     | spores                     |        |
+---------------------------------------------------------------------+----------------------------+--------+
| ``spores, germinated spores, microcolonies long cells,              | germinated, microcolonies, | 3      |
| occasionally misshapen branched``                                   | spores                     |        |
+---------------------------------------------------------------------+----------------------------+--------+
| ``ESSENTIAL spores, germinated spores``                             | germinated, spores         | 2      |
+---------------------------------------------------------------------+----------------------------+--------+
| ``germinated spores, microcolonies``                                | germinated, microcolonies  | 3      |
+---------------------------------------------------------------------+----------------------------+--------+

**Standalone-spores logic (germinated spores → spores is implied)**

+---------------------------------------------------------------------+---------------------+--------+
| Description                                                         | Category            | Tier   |
+=====================================================================+=====================+========+
| ``germinated spores long`` → spores appears only as part of         | germinated          | 2      |
| "germinated spores" → no standalone spores                          |                     |        |
+---------------------------------------------------------------------+---------------------+--------+
| ``spores, germinated spores`` → standalone "spores" exists →        | germinated, spores  | 2      |
| mixed population                                                    |                     |        |
+---------------------------------------------------------------------+---------------------+--------+

**Case / whitespace robustness**

+---------------------------------------------------------------------+---------+--------+
| Description                                                         | Categ.  | Tier   |
+=====================================================================+=========+========+
| ``Essential Spores`` (mixed case)                                   | spores  | 1      |
+---------------------------------------------------------------------+---------+--------+
| ``  ESSENTIAL spores  `` (leading/trailing space)                   | spores  | 1      |
+---------------------------------------------------------------------+---------+--------+

---

Test file: ``tests/test_classification.py``
--------------------------------------------

Tests for three functions: ``filter_primary_segments()``,
``count_growth_segments()``, and ``classify_phenotype_count()``.

**Filter — which comma segments are kept?**

The rule: a comma segment that **starts with** a modifier word is a
secondary description and is dropped.  All other segments are kept.

+---------------------------------------------------------------------+------------------------------------------------+--------+
| Input                                                               | After ``filter_primary_segments()``            | Reason |
+=====================================================================+================================================+========+
| ``VIABLE WT cells, some germination long``                          | ``VIABLE WT cells``                            | ``some`` is a modifier → |
|                                                                     |                                                | ``some germination long`` dropped |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``ESSENTIAL spores, germinated spores, some division``              | ``ESSENTIAL spores, germinated spores``        | ``some division`` dropped |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``spores, germinated spores, microcolonies long cells,              | ``spores, germinated spores, microcolonies     | ``occasionally misshapen`` |
| occasionally misshapen branched``                                   | long cells``                                   | dropped |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``ESSENTIAL germinated spores, slightly misshapen cells``           | ``ESSENTIAL germinated spores``                | ``slightly misshapen`` dropped |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``ESSENTIAL germinated spores very large swollen branched,          | ``ESSENTIAL germinated spores very large       | ``initially germinated`` |
| initially germinated spores long``                                  | swollen branched``                             | dropped |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``ESSENTIAL spores, germinated spores, microcolonies slightly       | ``ESSENTIAL spores, germinated spores,         | ``microcolonies slightly`` |
| misshapen cells``                                                   | microcolonies slightly misshapen cells``       | is KEPT (``microcolonies`` |
|                                                                     |                                                | not a modifier → only     |
|                                                                     |                                                | ``slightly misshapen cells``|
|                                                                     |                                                | starting with ``slightly`` |
|                                                                     |                                                | would be dropped)          |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``ESSENTIAL germinated spores, divide once``                        | ``ESSENTIAL germinated spores, divide once``   | ``divide`` is not a       |
|                                                                     |                                                | modifier → kept           |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``VIABLE cells, occasionally long, slightly misshapen``             | ``VIABLE cells``                                | All other segments start  |
|                                                                     |                                                | with modifiers → kept     |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``VIABLE small colonies long cells, possibly diploidising``         | ``VIABLE small colonies long cells``           | ``possibly diploidising`` |
|                                                                     |                                                | dropped                   |
+---------------------------------------------------------------------+------------------------------------------------+--------+
| ``ESSENTIAL germinated spores, often divide once``                  | ``ESSENTIAL germinated spores``                | ``often divide once``     |
|                                                                     |                                                | dropped                   |
+---------------------------------------------------------------------+------------------------------------------------+--------+

**Count growth segments — how many comma segments contain a growth signal?**

+---------------------------------------------------------------------+-------+
| Description                                                         | Count |
+=====================================================================+=======+
| ``ESSENTIAL germinated spores`` (no comma, 1 signal)                | 1     |
+---------------------------------------------------------------------+-------+
| ``VIABLE WT cells`` (no comma, 0 signals)                           | 0     |
+---------------------------------------------------------------------+-------+
| ``spores, germinated spores``                                       | 2     |
+---------------------------------------------------------------------+-------+
| ``spores, germinated spores, microcolonies, occasionally            | 3     |
| misshapen``                                                         |       |
+---------------------------------------------------------------------+-------+
| ``germinated spores slightly misshapen, some division``             | 1     |
| (``some division`` is modifier-led → not counted)                   |       |
+---------------------------------------------------------------------+-------+
| ``germinated spores, divide once`` (``divide`` not a modifier)      | 2     |
+---------------------------------------------------------------------+-------+
| ``np.nan`` (non-string input)                                       | 0     |
+---------------------------------------------------------------------+-------+
| ``ESSENTIAL spores, germinated spores, often divide once``          | 2     |
| (``often divide once`` is modifier-led → not counted)               |       |
+---------------------------------------------------------------------+-------+

**Single/Multiple classification — ``classify_phenotype_count()``**

+---------------------------------------------------------------------+----------------------+----------+
| Description                                                         | Result               | Reason   |
+=====================================================================+======================+==========+
| ``ESSENTIAL germinated spores`` (no comma)                          | Single               | 1 signal |
+---------------------------------------------------------------------+----------------------+----------+
| ``spores, germinated spores``                                       | **Multiple**         | 2 signals|
+---------------------------------------------------------------------+----------------------+----------+
| ``spores, germinated spores, microcolonies, occasionally            | **Multiple**         | 3 signals|
| misshapen``                                                         |                      |          |
+---------------------------------------------------------------------+----------------------+----------+
| ``ESSENTIAL germinated spores slightly misshapen, some division``   | Single               | 1 signal |
| (``some division`` is modifier-led → skipped)                       |                      |          |
+---------------------------------------------------------------------+----------------------+----------+
| ``ESSENTIAL spores, germinated spores slightly misshapen, some      | **Multiple**         | 2 signals|
| division`` (spores + germinated = 2, some division → skipped)       |                      |          |
+---------------------------------------------------------------------+----------------------+----------+
| ``ESSENTIAL microcolonies skittle cells, small colonies WT cells``  | **Multiple**         | 2 signals|
+---------------------------------------------------------------------+----------------------+----------+
| ``VIABLE WT cells, some germination long``                          | Single               | 0 signals|
+---------------------------------------------------------------------+----------------------+----------+
| ``ESSENTIAL germinated spores, divide once``                        | **Multiple**         | 2 signals|
+---------------------------------------------------------------------+----------------------+----------+
| ``VIABLE WT cells`` (no comma)                                      | Single               | 0 signals|
+---------------------------------------------------------------------+----------------------+----------+
| ``VIABLE misshapen cells`` (no comma, only morphology)              | Single               | 0 signals|
+---------------------------------------------------------------------+----------------------+----------+
| ``ESSENTIAL spores`` (no comma)                                     | Single               | 1 signal |
+---------------------------------------------------------------------+----------------------+----------+
| ``ESSENTIAL microcolonies`` (no comma)                              | Single               | 1 signal |
+---------------------------------------------------------------------+----------------------+----------+

**Constant integrity checks**

+---------------------------------------------------------------------------------+
| Check                                                                           |
+=================================================================================+
| MODIFIER_WORDS has no duplicates                                                |
+---------------------------------------------------------------------------------+
| GROWTH_KEYWORDS has no duplicates                                               |
+---------------------------------------------------------------------------------+
| All MODIFIER_WORDS entries are lowercase                                        |
+---------------------------------------------------------------------------------+
| All GROWTH_KEYWORDS entries are lowercase                                       |
+---------------------------------------------------------------------------------+

---

Running the tests
-----------------

.. code-block:: bash

    # Full suite
    mamba run -n bioinformatics python -m pytest tests/ -v

    # Single file
    mamba run -n bioinformatics python -m pytest tests/test_growth_signals.py -v

    # Single test class
    mamba run -n bioinformatics python -m pytest tests/ -k "TestClassifyPhenotypeCount" -v

    # Stop on first failure
    mamba run -n bioinformatics python -m pytest tests/ -x