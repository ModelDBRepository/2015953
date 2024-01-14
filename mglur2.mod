TITLE detailed model of mGluR2 receptors

COMMENT
-----------------------------------------------------------------------------

  Model of mGluR2 receptors modified from GABA-B receptor model 
  described in:

  Destexhe, A. and Sejnowski, T.J.  G-protein activation kinetics and
  spill-over of GABA may account for differences between inhibitory responses
  in the hippocampus and thalamus.  Proc. Natl. Acad. Sci. USA  92:
  9515-9519, 1995.

  See also: 

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1998, pp. 1-25.

  Written by Alain Destexhe, Laval University, 1995

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS mglur2
	POINTER C
	RANGE R, D, G, g, gmax
	NONSPECIFIC_CURRENT i
	GLOBAL K1, K2, K3, K4, KD, d1, d2, Erev
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

:
:	From simplex fitting to experimental data
:	(Destexhe and Sejnowski, 1995)
:
	K1	= 0.66	(/ms mM)	: forward binding rate to receptor
	K2	= 0.020 (/ms)		: backward (unbinding) rate of receptor
	K3	= 0.083 (/ms)		: rate of G-protein production
	K4	= 0.0079 (/ms)		: rate of G-protein decay
	d1	= 0.017 (/ms)		: rate of desensitization
	d2	= 0.0053 (/ms)		: rate of re-sensitization
	KD	= 100			: dissociation constant of K+ channel
	n	= 4  : nb of binding sites of G-protein on K+
	Erev   = -90 (mV)		: reversal potential (E_K)
	gmax   (umho)		: maximum conductance
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: pointer to transmitter concentration
	Gn
}


STATE {
	R				: fraction of activated receptor
	D				: fraction of desensitized receptor
	G				: fraction of activated G-protein
}


INITIAL {
	R = 0
	D = 0
	G = 0
}

BREAKPOINT {
	SOLVE bindkin METHOD cnexp
	Gn = G^n
	g = gmax * Gn / (Gn+KD)
	i = g*(v - Erev)
}


DERIVATIVE bindkin {

	R' = K1 * C * (1-R-D) - K2 * R + d2 * D
	D' = d1 * R - d2 * D
	G' = K3 * R - K4 * G

}


