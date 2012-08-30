package sim;

public abstract class Particle {
	
	Coord r, rNew;
	Particle partner = null;
	int rank = 0;
	// How many levels of secondaries? 0 -> no secondaries
	static int maxRank;
	static double cutoffE;
	double phi, theta;
	double energy;
	
	/**
	 * Step length between scattering events
	 * @param effectiveElectronDensity
	 * @param cutoffE
	 * @param atomDensities
	 * @param atomicNumbers
	 * @return
	 */
	abstract double getStepLength(double effectiveElectronDensity, 
			double[] atomDensities, double[] atomicNumbers);
	

	/**
	 * Will one or more secondaries be generated? This is random
	 * @param randUnif a uniformly distributed random number between 0 and 1
	 * @param effectiveElectronDensity
	 * @param cutoffE
	 * @param atomDensities
	 * @param atomicNumbers
	 * @return
	 */
	abstract boolean willSpawn(double randUnif, 
			double effectiveElectronDensity,
			double[] atomDensities, double[] atomicNumbers);
	
	/**
	 * Makes a collection of particles, each linking to another
	 * except for a "root", which is returned.
	 * @return
	 */
	abstract Particle spawn();

	/**
	 * Elastically scatter this, resulting in a change of theta and phi
	 * @param atomDensities
	 * @param atomicNumbers
	 */
	abstract void scatter(double[] atomDensities, double[] atomicNumbers);

	/**
	 * Stopping power, in eV/Angstrom
	 * @param effectiveElectronDensity
	 * @param meanExcitationEnergy
	 * @return
	 */
	abstract public double dE(double effectiveElectronDensity,
			double meanExcitationEnergy);

}
