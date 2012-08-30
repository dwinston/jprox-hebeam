package sim;

import java.util.Random;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.Math.asin;
import static java.lang.Math.acos;
import static sim.Standard.PI;
import static sim.Standard.ESU;

public class Electron extends Particle {

	@Override
	double getStepLength(double effectiveElectronDensity, 
			double[] atomDensities, double[] atomicNumbers) {
		return (-lambda(effectiveElectronDensity, cutoffE, 
				atomDensities, atomicNumbers)
				* log(new Random().nextDouble()));
	}


	/**
	 * mean free path, in angstroms
	 * @param effectiveElectronDensity
	 * @param cutoffE
	 * @param atomDensities
	 * @param atomicNumbers
	 * @return
	 */
	private double lambda(double effectiveElectronDensity, double cutoffE, 
			double[] atomDensities, double[] atomicNumbers) {
	    double sum = 1/lambda_el(atomDensities, atomicNumbers);
	    if (maxRank > 0) {
	    	sum += 1/lambda_inel(effectiveElectronDensity, cutoffE);
	    }
	    return 1/sum;
	}

	/**
	 * mean free path, in angstroms, for inelastic scattering
	 * @param effectiveElectronDensity
	 * @param cutoffE
	 * @return
	 */
	private double lambda_inel(double effectiveElectronDensity, 
			double cutoffE) {
	    return 1/(effectiveElectronDensity
	    		* sigma_inel(cutoffE));
	}

	/**
	 * Based on Murata (JAP '81)
	 * 
	 * Based on classical cross-section
	 * 
	 *                PI e^4    1-2e_c
	 * sigma_inel = -------- ---------- (E in ergs, e in ESU, sigma in cm^2)
	 *                E^2    e_c(1-e_c)
	 * @param cutoffE cutoff energy for lower bound of integration.               
	 * @return
	 */
	private double sigma_inel(double cutoffE) {
	    double result, sigma_inel_prefactor;
	    double e_c = cutoffE/energy;
	    double Eerg=energy*1.60219e-12;
	    sigma_inel_prefactor = (PI*ESU*ESU*ESU*ESU 
	    		*(1.0 - 2.0*e_c)/(e_c*(1.0 - e_c)));
	    result = sigma_inel_prefactor/(Eerg*Eerg);
	    return(result*1.0e16); /* convert to square Angstroms */
	}


	/**
	 * mean free path, in angstroms, for elastic scattering
	 * @param atomDensities
	 * @param atomicNumbers
	 * @return
	 */
	private double lambda_el(double[] atomDensities,
			double[] atomicNumbers) {
		double sum = 0;
		for (int i=0; i < atomDensities.length; i++) {
	        sum += (atomDensities[i] * sigma_el(atomicNumbers[i]));
	    }
	    return 1/sum;
	}

	/**
	 * from RJH thesis pp.5-6 (corrected to include e^4 by RG 3/19/92)
     *
     *                      Z(Z+1)PI e^4
     * sigma_elastic = -------------------------- (E in ergs, e in ESU, sigma in cm^2)
     *                   4E^2alpha^2(alpha^2+1)
     *                   
	 * @param Z atomic number
	 * @return
	 */
	private double sigma_el(double Z) {
	    double asqr = alphasqr(Z);
	    double result, sigma_el_prefactor;
	    double Eerg=energy*1.60219e-12;

	    sigma_el_prefactor = PI*ESU*ESU*ESU*ESU/4.0;
	    result = sigma_el_prefactor*Z*(Z+1.0)/(Eerg*Eerg*asqr*(asqr + 1.0));

	    return(result*1.0e16); /* convert to square Angstroms */
	}

	/**
	 * alpha^2 = (2.33 Z^(1/3) / sqrt(E))^2 (E in eV,alpha unitless)
	 * 
	 * @param Z atomic number
	 * @return
	 */
	private double alphasqr(double Z) {
	    return 2.33*2.33*pow(Z,2.0/3.0)/energy;
	}


	@Override
	boolean willSpawn(double randUnif, double effectiveElectronDensity,
			double[] atomDensities, double[] atomicNumbers) {
		double lamb = lambda(effectiveElectronDensity, cutoffE, 
				atomDensities, atomicNumbers);
		double lamb_el = lambda_el(atomDensities, atomicNumbers);
		return (randUnif > (lamb/lamb_el));
	}


	@Override
	/**
	 * Returns a secondary with a link to a scattered this as its partner.
	 */
	Particle spawn() {
		double ep=epsilon();
		/* energy normalized to rest energy of electron (eV) */
        double tau=energy/511000.0; 

        /* Assign direction and energy for this */
        theta = asin(sqrt(2.0*ep/(2.0 + tau*(1.0-ep))));
        energy = energy*(1.0-ep);
        phi = 2.0 * PI * Simulation.rand.nextDouble();

        /* Assign direction and energy for secondary */
        Particle sec = new Electron();
        sec.theta = asin(sqrt((2.0*(1.0-ep))/(2.0+tau*ep)));
        sec.energy = ep*energy;
        sec.phi = phi - PI;
        sec.rank = rank + 1;
        sec.partner = this;
        sec.r = r; sec.rNew = rNew;
        
		return sec;
	}

	/**
	 * From Joy '83
	 * @return fraction of energy transfer to secondaries (0<ep<0.5)
	 */
	private double epsilon() {
	    double e_c = cutoffE/energy;
	    double A = ((Simulation.rand.nextDouble() - 1.0)
	    		* (1.0-2.0*e_c)/(e_c*(1.0-e_c)));
	    double result = (-2.0 + A + sqrt(4.0 + A*A))/(2.0*A);
	    if (result<=0.0 || result>0.5) {
	        System.err.println("Error: ep messed up:  ep="+result);
	    }
	    return (result);
	}


	@Override
	void scatter(double[] atomDensities, double[] atomicNumbers) {
		double Z = getScatteringAtom(atomDensities, atomicNumbers);
		phi = 2.0 * PI * Simulation.rand.nextDouble();
        double R = Simulation.rand.nextDouble();
        double alphai = alphasqr(Z);
        theta = acos((R*(1.0+2.0*alphai) - alphai)/(R + alphai));
	}


	private double getScatteringAtom(double[] atomDensities,
			double[] atomicNumbers) {
	    double R =  Simulation.rand.nextDouble();
	    double sum = 0.0;
	    double running = 0.0;
	    int numSpecies = atomDensities.length;
	    double[] P = new double[numSpecies];
	    double result = -1;
	    for (int i = 0 ; i < numSpecies ; i++) {
	        P[i] = atomDensities[i] * sigma_el(atomicNumbers[i]);
	        sum += P[i];
	    }
	    for (int i = 0; i < numSpecies; i++) {
	        running += P[i];
	        if (R < (running/sum)) {
	            result = atomicNumbers[i];
	            break;
	        }
	    }
	    if (result < 0) {
	        System.err.println("Error in getScatteringAtom");
	    }
	    return (result);
	}


	@Override
	public double dE(double effectiveElectronDensity,
			double meanExcitationEnergy) {
		/*
	     * from RJH thesis pp.7-8
	     *
	     *   dE             2PIe^4n_e       E
	     * (----)      = - ----------  ln(---- sqrt(e(base of ln)/2) )
	     *   ds  Bethe         E           <I>
	     *
	     * CHANGE--JF 11/19/08: Use Joy and Luo's modification (Scanning-1989)
	     * 
	     *   dE             2PIe^4n_e             E+0.85*<I>
	     * (----)      = - ----------  ln(1.166*------------- )
	     *   ds  Bethe         E                     <I>
	     * 
	     * 
	     * from Murata (JAP '81)
	     *
	     *   dE                 /0.5      dsigma      
	     * (----)       = - n_e |    Eep(--------)dep
	     *   ds  Single         /e_c       dep        
	     *
	     *                   Pie^4n_e          1
	     *              = - ---------- (2 - ------- - ln(4e_c(1-e_c)))
	     *                      E            1-e_c
	     *
	     *  dE returned as eV/angstrom
	     */
	    double J = meanExcitationEnergy;
	    double result = -(1302.83*effectiveElectronDensity/energy)
	                    * log(1.166 * (energy+0.85*J)/J );
	    if (rank != maxRank) {
	        /* this spawns secondaries */
	    	double e_c = cutoffE/energy;
	        result += (651.415*effectiveElectronDensity/energy)
	                  *(2.0 - 1.0/(1.0-e_c) - log(4.0*e_c*(1.0-e_c)));
	    }
	    return (result);
	}
}
