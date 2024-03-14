=====================================================================
				MC_coincidences.C
=====================================================================

Macro for simulating number of accidental events using Monte-Carlo 
method.


Macro takes following arguments:

	1) activity - variable of type double; defines activity of the 
			source of radiation, given in [Bq]. Defaut value is 1 MBq.
			
	2) time - defines time length of measurement in [s].
			Default value is 1 ms.

	3) timeFrame - parameter defining maximum time difference in [s] 
			between	mached creation and decay of ortho-positronium. 
			If time difference exceeds this parameter, the match is 
			discarded. When set to values <= 0.0, time difference is 
			ignored when matching.
			Default value is 0.0 
