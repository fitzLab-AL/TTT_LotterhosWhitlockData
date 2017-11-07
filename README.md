# fitzLab-AL_TTT_LotterhosWhitlockData

For each simulation, the true selective environment is included in the `simfiles` folder.

In the file `results_AdaptreeEnviFor_R90.txt` the following columns mean:
"PopID" 
"X_Pops" x-location
"Y_Pops" y-location
"R90" "R60" "R30" "NSRangeTrans_30" These columns are for different sampling designs not included and can be ignored.
"Elevation" 
"MAT" Mean annual temperature (overlaid)
"MWMT" Mean warmest month temperature (overlaid)
"MCMT" Mean coldest month temperature (overlaid)
"TD" Temperature difference (overlaid)
"log.MAP." log-Mean annual precipitation (overlaid)
"log.MSP." log-Mean spring precipitation (overlaid)
"AHM" Annual heat:moisture index (overlaid)
"SHM" Summer heat:moisture index (overlaid)
"DD0" Degree-days below 0c
"DD5" Degree-days above 5c
"NFFD" 
"bFFP" begin frost-free period
"eFFP" end frost-free period
"FFP" frost-free period
"PAS" precipitation as snow
"EMT" extreme minimum temperature
"EXT" extreme maximum temperature
"Eref" 
"CMD" climate moisture deficit

The overlaid environmental variables are from real data (Adaptree project, British Colubmia and Alberta), which has been overlaid onto the simulations and interpolated.
KEL has confirmed that the interpolation preserves the correlation structure among the environmental variables.
