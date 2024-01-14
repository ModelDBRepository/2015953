COMMENT
Based on diffusion equation from Nmda.mod from 
"LTP regulates burst initiation and frequency at mossy fiber-granule
cell synapses of rat cerebellum: experimental observations and
theoretical predictions." Nieus T, Sola E, Mapelli J, Saftenku E,
Rossi P, D'Angelo E.  J Neurophysiol. 2006 Feb;95(2):686-99"

Originally published in Balmer TS, Borges-Merjane C, Trussell LO (2021) 
Incomplete removal of extracellular glutamate controls synaptic transmission and 
integration at a cerebellar synapse. eLife 10:e63819.
  
ENDCOMMENT

NEURON {
    POINT_PROCESS diff3D
    RANGE M, Diff, R, lambd, T, alpha, amb
}

UNITS {
    PI	= (pi)		(1)
    (mM) = (milli/liter)
}

PARAMETER {
	: with the following values a single spike should produce a transient (T) that matches Fig 2C in Barbour & Hausser 1997 TINS
	M		= 4700				 
	R		= 1100 (nm)
	Diff		= 0.76 (um2/ms)
	lambd		= 1.55 
    alpha = 0.21
    amb = 0.005 (mM)
}

ASSIGNED {
    T (mM)
	tspike[5000]	(ms) 
	tsyn		(ms)
	NTdiffusion	(mM)
	numpulses
}

INITIAL {
    T=0 (mM)
    tspike[0]=1e12	(ms)
	numpulses = 0
}

FUNCTION NTdiffWave(){
	LOCAL ijk,t0
	: sums up diffusion contributes
	NTdiffusion=0
	FROM ijk=1 TO numpulses{
		t0=tspike[ijk-1]
		if(t>t0){		
 		NTdiffusion = NTdiffusion+((exp(-((1e-8*R)^2)/(4*Diff*1e-10*(t - t0)/(lambd^2))))*(M/6.022e23)/(8*alpha*((PI*Diff*1e-10*(t - t0)/(lambd^2))^1.5) ) )*1e3 
        }
	}					
	NTdiffWave=NTdiffusion
}



BREAKPOINT {
    T = NTdiffWave() + amb
}

:::::::::::::::::::::::
NET_RECEIVE(weight, on, t0 (ms)) {
		if (!on) {
			t0 = t
			on = 1				

			tspike[numpulses] = t
			numpulses = numpulses + 1
		}
		on = 0
}
