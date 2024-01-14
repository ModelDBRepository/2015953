TITLE I-h channel for Thalamic neurons from McCormick and Pape (1990) 
: modified to match values of UBC in Subramaniyam et al 2014

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v  (mV)
    eh = -30    (mV)
	celsius    (degC)
	ghbar=.00005   (mho/cm2)
    vhalft = -91.5  (mV) :this is based on 
    a0t=0.0005  (/ms)
    zetat=0.2   (1)
    gmt=.65 (1)
	q10=4.5
}

NEURON {
	SUFFIX htc
	NONSPECIFIC_CURRENT i
        RANGE ghbar, eh
        GLOBAL linf,taul
}

STATE {
    l
}

ASSIGNED {
	i  (mA/cm2)
    linf      
    taul
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*l*(v-eh)
}

FUNCTION alpt(v(mV)) {
  alpt = exp(zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(zetat*gmt*(v-vhalft)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v)
        l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL qt
        qt=q10^((celsius-36)/10)
        linf = 1/(1+exp((v+91.5)/7.7))
        taul = bett(v)/(qt*a0t*(1+alpt(v)))
}




















