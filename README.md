** Description of app goes here **


DEVELOPMENT NOTES

As of August 14, 2025:

The following layout changes were made:
- The layout is now 3 overarching tabs with the same sidear. In an earlier version, a gene would have to be selected in the "Data Visualization" tab's sidebar before the "Gene Expression Heatmaps" tab could be accessed.
- The dropdown menus in the "Gene Expression Heatmaps" tab have been rearranged

The following glitches were found and fixed:
- The "Gene Expression Heatmaps" tab crashes if no gene is selected.
- If a second gene is selected before the first selected gene fully loads, the app flickers, freezes and then crashes.
- If a gene is deselected from the list, the associated boxplot does not always get deleted.
- When the non-deletion glitch was fixed, another glitch arose: selected a second dataset does not always produce a sub-tab under "Data Visualization" for that dataset.

The following glitches remain:
- If a dataset is selected, then deselected, the app crashes.
